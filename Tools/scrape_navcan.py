from playwright.sync_api import sync_playwright, Playwright
from tkinter import filedialog
import os
from csv import writer
from collections import defaultdict


def uncheck_buttons(page, toggles_to_uncheck):
    """
    div.product-inputs
        └── div.product-input
            └──label.{toggle_name}-toggle
                └──input.id="upperwind-toggle"
    """

    for toggle in toggles_to_uncheck:
        checkbox = page.locator(f"#{toggle}")
        checkbox.uncheck()
        print(f"{toggle} is unchecked")



def scraper(playwright: Playwright, site: str, toggles:list[str]):

    # Default browser due to chrome
    chromium = playwright.chromium
    browser = chromium.launch()
    page = browser.new_page()
    print("Page made")

    # Current instance is on the correct browser page
    page.goto(site)
    print("Page entered")

    # Make sure all buttons are unchecked
    uncheck_buttons(page,toggles)
    print("Buttons unchecked")

    # Make sure the CYYU airport is in the search bar
    page.get_by_placeholder("Enter Aerodrome, FIR, Navaid, etc.").fill("CYYU")
    print("CYYU entered")

    # Click on the search button to initialize the search
    page.click("div.btn.btn-primary.search-button")
    print("Click 1")

    # Click on the search button to find the weather information
    page.click("div.btn.btn-primary.search-button")
    print("Click 2")

    # Parse the website for the specific "pre" container
    # Inner text extracts the text as a string
    raw_weather = page.locator("pre").inner_text()

    browser.close()
    return raw_weather


def run_scraper():
    site = "https://plan.navcanada.ca/wxrecall/"
    toggles_to_uncheck = [
            "sigmet-toggle", "airmet-toggle",
            "notam-toggle", "metar-toggle",
            "taf-toggle", "pirep-toggle",
            "space_weather-toggle",
        ]
    with sync_playwright() as playwright:
        raw_weather_string = scraper(playwright, site, toggles_to_uncheck)

    return raw_weather_string.split("VALID")

class WindProfile:
    def __init__(self, release_date: str, weather_window: str, hours: int):
        self.release_date = release_date
        self.weather_window = weather_window
        self.hours = hours
        self.profile = defaultdict(tuple)

    # For testing purposes
    def __str__(self):
        lines = [
            f"Release Date   : {self.release_date}",
            f"Weather Window : {self.weather_window}",
            f"Duration       : {self.hours} hours",
            "Wind Profile:",
        ]

        for altitude in sorted(self.profile):
            direction, speed = self.profile[altitude]
            lines.append(f"  {altitude:>5} ft : {direction}° @ {speed} kt")

        return "\n".join(lines)


def clean_raw_weather(raw_weather: list[str]):
    wind_profiles = []
    scenarios = []
    ALTITUDES = [
        3000, 6000, 9000, 12000,
        18000, 24000, 30000,
        34000, 39000, 45000, 53000
    ]

    for case in raw_weather:
        # The case of the empty string
        if not case:
            continue
        else:
            scenarios.append(case.split("\n")[:3])

    # Each window is a list
    for window in scenarios:
        for section in window:
            if not section:
                continue
            else:
                section = section.strip()

                # If this is a header row
                if "FOR USE" in section:
                    parts = section.split()

                    release_date = parts[0].rstrip("Z")  # "251800"
                    weather_window = parts[-1]  # "14-21"
                    start_hour, end_hour = weather_window.split("-")  # "14", "21"
                    hours = int(start_hour) - int(end_hour)
                    if end_hour < start_hour:
                        hours -= 24

                    altitude_profile = WindProfile(release_date, weather_window, -hours)

                # If this is a wind data row
                elif "YYU" in section:
                    section = section.strip()[3:]

                    parts = section.split("|")

                    for altitude, data in zip(ALTITUDES, parts):
                        part = data.strip()
                        degree = int(part[:3])
                        speed = int(part[4:7])

                        altitude_profile.profile[altitude] = (degree, speed)

        wind_profiles.append(altitude_profile)

    return wind_profiles

def create_csvs(weather_objects: list[WindProfile], root_path):
    BASE_DIR = root_path
    file_path = os.path.join(BASE_DIR, "csv_files")
    os.makedirs(file_path, exist_ok=True)

    for weather_object in weather_objects:
        name = f"OpenRocket_Polaris_Cycle6_Params_{weather_object.release_date}_{weather_object.hours}hourwindow.csv"
        full_path = os.path.join(file_path, name)

        with open(full_path, "w", newline="") as file:
            csv_writer = writer(file)
            csv_writer.writerow(["altitude", "speed", "direction"])

            for altitude in sorted(weather_object.profile):
                direction, speed = weather_object.profile[altitude]
                csv_writer.writerow([altitude, speed, direction])


def main():
    raw_weather_list = run_scraper()
    weather_profile_objects = clean_raw_weather(raw_weather_list)
    print("Opening folder picker window...")


    # Open a native OS window to visually browse and pick a folder
    file_path = filedialog.askdirectory(title="Select Destination Folder")

    # Check if the user closed the window without selecting anything
    if not file_path:
        print("Selection cancelled. Files were not saved.")
        return

    try:
        create_csvs(weather_profile_objects, file_path)
        print(f"Success! Files saved to: {file_path}")
    except Exception as e:
        print("There was an error: ", e)


main()
