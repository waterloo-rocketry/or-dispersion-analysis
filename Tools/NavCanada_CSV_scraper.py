"""
NAV Canada WXRecall wind-profile scraper.

Workflow
--------
1. Scrape the raw upper-wind text block for a given aerodrome from
   plan.navcanada.ca/wxrecall.
2. Parse it into WindProfile objects (one per forecast window).
3. Optionally attach speed/direction standard-deviation values looked
   up from OpenMeteo-derived histogram tables.
4. Write one CSV per WindProfile to a user-selected folder.
"""

import os
from csv import writer
from tkinter import filedialog

import numpy as np
import pandas as pd
from playwright.sync_api import Playwright, sync_playwright

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

SITE_URL = "https://plan.navcanada.ca/wxrecall/"

TOGGLES_TO_UNCHECK = [
    "sigmet-toggle",
    "airmet-toggle",
    "notam-toggle",
    "metar-toggle",
    "taf-toggle",
    "pirep-toggle",
    "space_weather-toggle",
]

ALTITUDES = [3000, 6000, 9000, 12000, 18000, 24000, 30000, 34000, 39000, 45000, 53000]

HEATMAP_DATA_PATH = "/Users/ricechessmaster/NavCanada_Scraper/data/heatmap_historical_2021-2025.csv"

# altitude -> (low_weight, high_weight, low_openmeteo_alt, high_openmeteo_alt)
ALT_LOOKUPS = {
    3000: (0.4280001463039954, 0.5719998536960046, 800, 1000),
    6000: (0.17800014630399508, 0.8219998536960049, 1500, 1900),
    9000: (0.3513846829095363, 0.6486153170904637, 1900, 3200),
    12000: (0.5424001170431962, 0.4575998829568037, 3200, 4200),
    18000: (0.08114298254628205, 0.918857017453718, 4200, 5600),
    24000: (0.942400117043196, 0.05759988295680392, 7200, 9200),
    30000: (0.028000146303995166, 0.9719998536960048, 7200, 9200),
    34000: (0.030666943018657133, 0.9693330569813429, 9200, 10400),
    39000: (0.9487061061119922, 0.051293893888007765, 11800, 13500),
    45000: (0.906087147353037, 0.09391285264696309, 13500, 15800),
    53000: (0.8134739562846226, 0.18652604371537743, 15800, 17700),
}


# ---------------------------------------------------------------------------
# WindProfile model
# ---------------------------------------------------------------------------

class WindProfile:
    def __init__(self, release_date: str, weather_window: str, hours: int, start:str, end:str):
        self.release_date = release_date
        self.weather_window = weather_window
        self.hours = hours
        self.start_hour = start
        self.end_hour = end
        self.profile: dict[int, tuple[int, int]] = {}
        self.spd_std_dev: dict[int, float] = {}
        self.dir_std_dev: dict[int, float] = {}

    def __str__(self):
        lines = [
            f"Release Date   : {self.release_date}",
            f"Weather Window : {self.weather_window}",
            f"Duration       : {self.hours} hours",
            f'From           : {self.start_hour}-{self.end_hour}',
            "Wind Profile:",
        ]
        for altitude in sorted(self.profile):
            direction, speed = self.profile[altitude]
            lines.append(f"  {altitude:>5} ft : {direction}\u00b0 @ {speed} kt")
        return "\n".join(lines)


# ---------------------------------------------------------------------------
# Scraping
# ---------------------------------------------------------------------------

def uncheck_buttons(page, toggles_to_uncheck: list[str]) -> None:
    """
    div.product-inputs
        └── div.product-input
            └── label.{toggle_name}-toggle
                └── input#{toggle_name}
    """
    for toggle in toggles_to_uncheck:
        page.locator(f"#{toggle}").uncheck()
        print(f"{toggle} is unchecked")


def scrape_wind_profile(playwright: Playwright, site: str, toggles: list[str]) -> str:
    browser = playwright.chromium.launch()
    try:
        page = browser.new_page()
        print("Page made")

        page.goto(site)
        print("Page entered")

        uncheck_buttons(page, toggles)
        print("Buttons unchecked")

        page.get_by_placeholder("Enter Aerodrome, FIR, Navaid, etc.").fill("CYYU")
        print("CYYU entered")

        # First click focuses/opens results, second click runs the search.
        page.click("div.btn.btn-primary.search-button")
        print("Click 1")
        page.click("div.btn.btn-primary.search-button")
        print("Click 2")

        return page.locator("pre").inner_text()
    finally:
        browser.close()


def run_scraper() -> list[str]:
    with sync_playwright() as playwright:
        raw_weather = scrape_wind_profile(playwright, SITE_URL, TOGGLES_TO_UNCHECK)
    return raw_weather.split("VALID")


# ---------------------------------------------------------------------------
# Parsing
# ---------------------------------------------------------------------------

def clean_raw_weather(raw_weather: list[str]) -> list[WindProfile]:
    wind_profiles = []
    scenarios = [case.split("\n")[:3] for case in raw_weather if case]

    for window in scenarios:
        profile = None

        for section in window:
            section = section.strip()
            if not section:
                continue

            if "FOR USE" in section:
                parts = section.split()
                release_date = parts[0].rstrip("Z")       # e.g. "251800"
                weather_window = parts[-1]                # e.g. "14-21"
                start_hour, end_hour = (int(h) for h in weather_window.split("-"))

                hours = end_hour - start_hour
                if hours <= 0:
                    hours += 24  # window wraps past midnight

                profile = WindProfile(release_date, weather_window, hours, str(start_hour), str(end_hour))

            elif "YYU" in section and profile is not None:
                data = section[3:].strip()
                parts = data.split("|")

                for altitude, part in zip(ALTITUDES, parts):
                    part = part.strip()
                    direction = int(part[:3])
                    speed = int(part[4:7])
                    profile.profile[altitude] = (direction, speed)

        if profile is not None:
            wind_profiles.append(profile)

    return wind_profiles


# ---------------------------------------------------------------------------
# Standard-deviation lookups (optional enrichment step)
# ---------------------------------------------------------------------------

def generate_lookup_tables(filename: str = HEATMAP_DATA_PATH) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Builds intraday direction/speed standard-deviation heatmaps from
    historical data, binned by altitude (rows) and daily mean speed
    (columns, e.g. "0-10", "10-20", ...).
    """
    data = pd.read_csv(filename)

    speed_cols = ['110', '320', '500', '800', '1000', '1500', '1900', '3200', '4200', '5600', '7200', '9200',
                  '10400', '11800', '13500', '15800', '17700', '19300', '22000']
    dir_cols = [f'direction [{col}]' for col in speed_cols]
    altitudes = [int(col) for col in speed_cols]

    # Extract the day to separate macro (daily) from micro (intraday) variance
    data['day'] = pd.to_datetime(data['date'], utc=True).dt.date

    # Precompute unit vector components for accurate circular math
    for s_col, d_col in zip(speed_cols, dir_cols):
        rad = np.radians(data[d_col])
        data[f'_sin_{s_col}'] = np.sin(rad)
        data[f'_cos_{s_col}'] = np.cos(rad)

    g = data.groupby('day')

    # 1. Macro: daily means (used only to bin days by speed regime below)
    daily_speed_mean = g[speed_cols].mean()

    # 2. Micro: intraday standard deviations (hour-to-hour fluctuation)
    intraday_speed_std = g[speed_cols].std()
    intraday_dir_std = pd.DataFrame(index=daily_speed_mean.index)

    for s_col in speed_cols:
        sin_mean = g[f'_sin_{s_col}'].mean()
        cos_mean = g[f'_cos_{s_col}'].mean()
        R = np.hypot(sin_mean, cos_mean)
        R = np.where(R > 0, R, np.nan)
        intraday_dir_std[s_col] = np.degrees(np.sqrt(-2 * np.log(R)))

    # 3. Categorize intraday variability by daily speed regime (heatmaps)
    max_speed = daily_speed_mean.max().max()
    bins = np.arange(0, max_speed + 10, 10)
    bin_labels = [f"{int(bins[i])}-{int(bins[i + 1])}" for i in range(len(bins) - 1)]

    heatmap_dir = np.full((len(speed_cols), len(bin_labels)), np.nan)
    heatmap_spd = np.full((len(speed_cols), len(bin_labels)), np.nan)

    for i, s_col in enumerate(speed_cols):
        speed_bins = pd.cut(daily_speed_mean[s_col], bins=bins, labels=bin_labels, include_lowest=True)

        mean_intraday_dir = intraday_dir_std[s_col].groupby(speed_bins, observed=False).mean()
        mean_intraday_spd = intraday_speed_std[s_col].groupby(speed_bins, observed=False).mean()

        for j, label in enumerate(bin_labels):
            heatmap_dir[i, j] = mean_intraday_dir.get(label, np.nan)
            heatmap_spd[i, j] = mean_intraday_spd.get(label, np.nan)

    df_hm_dir = pd.DataFrame(heatmap_dir, index=altitudes, columns=bin_labels)
    df_hm_spd = pd.DataFrame(heatmap_spd, index=altitudes, columns=bin_labels)

    return df_hm_dir, df_hm_spd


def _speed_bin_label(speed: int, bin_labels: list[str]) -> str:
    idx = min(speed // 10, len(bin_labels) - 1)
    return bin_labels[idx]


def _perform_lookups(df: pd.DataFrame, low: int, high: int, speed: int) -> tuple[float, float]:
    label = _speed_bin_label(speed, df.columns.tolist())
    return df.loc[low, label], df.loc[high, label]


def lookup_values(
    df_hm_dir: pd.DataFrame,
    df_hm_spd: pd.DataFrame,
    weather_profiles: list[WindProfile],
) -> list[WindProfile]:
    """Attach speed/direction std-dev estimates to each WindProfile in place."""
    for profile in weather_profiles:
        for altitude, (direction, speed) in profile.profile.items():
            low_weight, high_weight, low_alt, high_alt = ALT_LOOKUPS[altitude]

            dir_low, dir_high = _perform_lookups(df_hm_dir, low_alt, high_alt, speed)
            spd_low, spd_high = _perform_lookups(df_hm_spd, low_alt, high_alt, speed)

            profile.dir_std_dev[altitude] = round(dir_low * low_weight + dir_high * high_weight, 1)
            profile.spd_std_dev[altitude] = round(spd_low * low_weight + spd_high * high_weight, 1)

    return weather_profiles


# ---------------------------------------------------------------------------
# CSV output
# ---------------------------------------------------------------------------

def _write_csv(profile: WindProfile, destination: str, scenario: str, bit: str, columns: list[str]) -> None:
    has_std_dev = len(columns) == 5
    name = f"{scenario}_{bit}_SimParams_released_{profile.release_date}_{profile.start_hour}-{profile.end_hour}_{profile.hours}h.csv"
    full_path = os.path.join(destination, name)

    with open(full_path, "w", newline="") as file:
        csv_writer = writer(file)
        csv_writer.writerow(columns)

        for altitude in sorted(ALTITUDES):
            direction, speed = profile.profile[altitude]
            row = [altitude, speed, direction]
            if has_std_dev:
                row += [profile.spd_std_dev[altitude], profile.dir_std_dev[altitude]]
            csv_writer.writerow(row)


def create_csvs(weather_profiles: list[WindProfile], destination: str, scenario: str) -> None:
    """
    Writes CSVs for each WindProfile into `destination`.

    Always writes the OPR file (altitude/speed/direction). Additionally
    writes a PLG file (with spd/dir std-dev columns) for any profile
    that has been enriched via lookup_values().
    """
    os.makedirs(destination, exist_ok=True)

    for profile in weather_profiles:
        _write_csv(profile, destination, scenario, "OPR", ["altitude", "speed", "direction"])

        has_std_dev = bool(profile.spd_std_dev) and bool(profile.dir_std_dev)
        if has_std_dev:
            _write_csv(
                profile,
                destination,
                scenario,
                "PLG",
                ["altitude", "speed", "direction", "spd_std_dev", "dir_std_dev"],
            )


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main():
    raw_weather_list = run_scraper()
    weather_profiles = clean_raw_weather(raw_weather_list)

    print("Building std-dev lookup tables from historical data...")
    df_hm_dir, df_hm_spd = generate_lookup_tables()
    weather_profiles = lookup_values(df_hm_dir, df_hm_spd, weather_profiles)

    print("Opening folder picker window...")
    destination = filedialog.askdirectory(title="Select Destination Folder")

    if not destination:
        print("Selection cancelled. Files were not saved.")
        return

    # NOTE: naming convention for manual delivery is:
    #   L-X_Polaris_Cycle6_mc_[recovery situation]_[mmddhhhh-hhhh]
    scenario = "L-X_Polaris_Cycle6"

    try:
        create_csvs(weather_profiles, destination, scenario)
        print(f"Success! Files saved to: {destination}")
    except Exception as e:
        print("There was an error:", e)

main()