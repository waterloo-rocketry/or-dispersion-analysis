import streamlit as st
from datetime import date, timedelta
from weather_router import run_many
import re
import pandas as pd
import numpy as np
from io import BytesIO
import zipfile
import config

# Launch site used when the app starts and when "Reset Coords" is pressed
default_launch_latitude = 47.965378
default_launch_longitude = -81.873536

st.title("Meteo Data to CSV")

unit_label = st.selectbox(
    "Wind speed unit",
    options=["Default (km/h)", "kn", "ms", "mph"],
    index=0,
)

# Empty string means "don't send the parameter", so Open-Meteo uses its km/h default
wind_speed_unit = "" if unit_label == "Default (km/h)" else unit_label

round_wind_direction = st.checkbox("Round wind direction")

st.markdown(
    """
    <p style="font-size:14px;">
        <em>
            Note: Checking this box rounds the wind direction for all queries to the nearest 10 degrees
        </em>
    </p>
    """,
    unsafe_allow_html=True
)

st.space()

DAILY_VARS = ["wind_gusts_10m_mean"]

# Request-only options passed straight through to Open-Meteo
EXTRA_PARAMS = {"models": "best_match",}
if wind_speed_unit:
    EXTRA_PARAMS["wind_speed_unit"] = wind_speed_unit


def sanitize_filename(name: str) -> str:
    """
    Make a string safe to use as a filename on any OS.

    Replaces every run of characters other than letters, numbers, dash and
    underscore with a single underscore. Falls back to "meteo" if nothing is left.
    """
    name = name.strip()
    name = re.sub(r"[^A-Za-z0-9_\-]+", "_", name)
    return name or "meteo"


def ensure_sorted_hourly(df: pd.DataFrame):
    """
    Return a copy of `df` with its time column parsed as UTC datetimes and sorted.

    The time column is whichever of "time" or "date" exists. Rows whose time
    can't be parsed are dropped.

    Returns:
        (sorted_df, time_column_name)
    """
    tcol = "time" if "time" in df.columns else ("date" if "date" in df.columns else None)
    if tcol is None:
        raise KeyError(f"Hourly df missing time column. Columns: {list(df.columns)}")

    out = df.copy()
    out[tcol] = pd.to_datetime(out[tcol], utc=True, errors="coerce")
    out = out.dropna(subset=[tcol]).sort_values(tcol).reset_index(drop=True)
    return out, tcol


def make_filtered_sorted(df_hourly: pd.DataFrame) -> pd.DataFrame:
    """
    Keep only the hours that are plausible launch conditions, sorted by time.

    Filters:
      - daylight only (is_day == 1)
      - hour of day between 08 and 19 inclusive
      - low cloud cover at or below 50 %
    """
    df, tcol = ensure_sorted_hourly(df_hourly)

    required = ["is_day", "cloud_cover_low", tcol]
    missing = [c for c in required if c not in df.columns]

    if missing:
        raise KeyError(
            f"Missing required columns for filtering: {missing}. "
            f"Have: {list(df.columns)}"
        )

    df = df[
        (df["is_day"] == 1) &
        (df["cloud_cover_low"] <= 50) &
        (df[tcol].dt.hour >= 8) &
        (df[tcol].dt.hour <= 19)
    ]

    return df.reset_index(drop=True)


# ----------------------------
# Coords UI
# ----------------------------
def coords():
    """
    Draw the latitude/longitude inputs and the "Reset Coords" button.

    Values live in st.session_state.latitude / .longitude so they survive reruns.
    """
    # First run: seed the inputs with the default launch site
    if "latitude" not in st.session_state:
        st.session_state.latitude = default_launch_latitude
    if "longitude" not in st.session_state:
        st.session_state.longitude = default_launch_longitude

    if st.button("Reset Coords"):
        st.session_state.latitude = default_launch_latitude
        st.session_state.longitude = default_launch_longitude
        st.rerun()

    col1, col2 = st.columns(2)
    with col1:
        st.number_input("Latitude", -90.0, 90.0, key="latitude", format="%.6f")
    with col2:
        st.number_input("Longitude", -180.0, 180.0, key="longitude", format="%.6f")


# ----------------------------
# Date routing logic
# ----------------------------
def decide_open_meteo_endpoint(
    start_date: date,
    end_date: date,
    today: date,
    historical_days: int = 90,
    forecast_horizon_days: int = 16,
):
    """
    Decide which Open-Meteo endpoint(s) a date range should be fetched from.

    Dates older than `historical_days` ago go to the historical endpoint; anything
    newer (up to `forecast_horizon_days` ahead) goes to the forecast endpoint. A range
    that crosses the cutoff is split into one historical and one forecast piece.

    Returns:
        dict with keys
          ok     - False if the range is invalid
          error  - message to show when ok is False
          mode   - "historical", "forecast" or "split"
          ranges - list of (kind, start, end) tuples for weather_router
    """
    if end_date < start_date:
        return {"ok": False, "error": "End date must be on/after start date.", "mode": None, "ranges": []}

    cutoff = today - timedelta(days=historical_days)
    forecast_max = today + timedelta(days=forecast_horizon_days)

    if end_date > forecast_max:
        return {"ok": False, "error": f"End date cannot be more than {forecast_horizon_days} days in the future.", "mode": None, "ranges": []}

    # Entirely historical
    if end_date < cutoff:
        return {"ok": True, "error": None, "mode": "historical", "ranges": [("historical", start_date, end_date)]}

    # Entirely forecast
    if start_date >= cutoff:
        return {"ok": True, "error": None, "mode": "forecast", "ranges": [("forecast", start_date, end_date)]}

    # Split: historical up to the day before the cutoff, forecast from the cutoff on
    hist_end = cutoff - timedelta(days=1)
    fcst_start = cutoff
    ranges = []
    if start_date <= hist_end:
        ranges.append(("historical", start_date, hist_end))
    if fcst_start <= end_date:
        ranges.append(("forecast", fcst_start, end_date))

    return {"ok": True, "error": None, "mode": "split", "ranges": ranges}


def format_query(item: dict) -> list[str]:
    """
    Describe a saved query's ranges for display, e.g. ["forecast: 2026-11-01 → 2026-11-05"].
    """
    ranges = item.get("ranges", [])
    parts = [f"{kind}: {s.isoformat()} → {e.isoformat()}" for kind, s, e in ranges]
    return parts


# ----------------------------
# Date input + Save/Remove UI
# ----------------------------
def date_ui(today: date):
    """
    Draw the date-range picker and the list of saved queries.

    The picked range is checked with decide_open_meteo_endpoint(); "Add Date" stores
    it in st.session_state.saved_queries, where each query gets an editable file
    name and a "Remove" button.
    """
    if "saved_queries" not in st.session_state:
        st.session_state.saved_queries = []

    picked = st.date_input("Select date range", value=(), key="date_range")

    # Normalize picked -> (start, end); the widget returns a 1-tuple mid-selection
    start_date, end_date = (None, None)
    if isinstance(picked, tuple) and len(picked) == 2:
        start_date, end_date = picked
    elif not isinstance(picked, tuple):
        start_date, end_date = picked, None  # single date picked

    # Only a valid, complete range can be added
    st.session_state.current_decision = None

    if start_date and end_date:
        decision = decide_open_meteo_endpoint(start_date, end_date, today=today)
        if decision["ok"]:
            st.session_state.current_decision = decision
        else:
            st.error(decision["error"])

    colA, colB = st.columns([0.25, 1.5])
    with colA:
        if st.button("Add Date", disabled=(st.session_state.current_decision is None)):
            decision = dict(st.session_state.current_decision)  # copy so later edits don't touch current_decision
            decision["label"] = ""  # file name, filled in by the text box below
            st.session_state.saved_queries.append(decision)
            st.rerun()

    with colB:
        if st.button("Clear all", disabled=(len(st.session_state.saved_queries) == 0)):
            st.session_state.saved_queries = []
            st.rerun()

    if st.session_state.saved_queries:
        st.subheader(f"Saved date queries ({len(st.session_state.saved_queries)})")

        for i, item in enumerate(st.session_state.saved_queries):

            with st.container(border=True):

                left, mid, right = st.columns([0.55, 0.30, 0.15])

                with left:
                    st.write(f"**Query {i+1}**")
                    st.caption(format_query(item))

                with mid:
                    key = f"label_{i}"
                    default_val = item.get("label", f"query_{i+1:02d}")
                    new_label = st.text_input(
                        "File name",
                        value=default_val,
                        key=key
                    )
                    item["label"] = sanitize_filename(new_label)

                with right:
                    # Two blank lines push the button down to line up with the text box
                    st.write("")
                    st.write("")
                    if st.button("Remove", key=f"rm_{i}"):
                        st.session_state.saved_queries.pop(i)
                        st.rerun()


# ----------------------------
# Run app
# ----------------------------
today = date.today()
coords()
date_ui(today)

st.divider()

st.subheader("Export settings")
st.markdown(
    """
    <p style="font-size:14px;">
        Note: Name of zipfile is "open-meteo-download_{current date}"
    </p>
    """,
    unsafe_allow_html=True
)


if st.button("Fetch weather"):
    saved = st.session_state.get("saved_queries", [])
    if not saved:
        st.warning("Add at least one date query first.")
        st.stop()

    lat = st.session_state.latitude
    lon = st.session_state.longitude

    with st.spinner("Fetching data..."):
        try:
            # One combined DataFrame for all queries, tagged with a 1-based query_id
            res = run_many(
                lat,
                lon,
                saved,
                hourly_vars=config.HOURLY_VARS,
                daily_vars=DAILY_VARS,
                extra_params=EXTRA_PARAMS,
            )
        except Exception as e:
            st.error(f"Open-Meteo request failed: {e}")
            st.stop()

    if "hourly" not in res or res["hourly"] is None or len(res["hourly"]) == 0:
        st.info("No hourly data returned.")
        st.stop()

    hourly_all = res["hourly"].copy()

    if "query_id" not in hourly_all.columns:
        st.error("Expected 'query_id' column in hourly output. Check weather_router.run_many().")
        st.stop()

    # Build the ZIP in memory so it can be handed straight to the download button
    zip_buffer = BytesIO()

    with zipfile.ZipFile(zip_buffer, "w", compression=zipfile.ZIP_DEFLATED) as zf:
        # Each query produces four CSVs: full, filtered, sim-params full, sim-params filtered
        for qid, decision in enumerate(saved, start=1):
            label = sanitize_filename(decision.get("label", f"query_{qid:02d}"))

            # This query's rows only
            df_q = hourly_all[hourly_all["query_id"] == qid].copy()
            df_q, tcol = ensure_sorted_hourly(df_q)

            if round_wind_direction:
                for col in config.WIND_DIRECTIONS:
                    df_q[col] = df_q[col].round(-1)  # round to the nearest 10 degrees

            # FULL + FILTERED
            df_q = df_q.rename(columns={tcol: 'date'})
            raw_filtered = make_filtered_sorted(df_q)

            # Format timestamps as plain text for the CSVs
            df_q["date"] = df_q["date"].dt.strftime("%Y-%m-%d %H:%M")
            raw_filtered["date"] = raw_filtered["date"].dt.strftime("%Y-%m-%d %H:%M")

            # Internal bookkeeping columns, not wanted in the exports
            del df_q["source"]
            del df_q["query_id"]

            zf.writestr(
                f"{label}_full-data.csv",
                df_q.to_csv(index=False)
            )

            # The filter columns are always the same after filtering, so drop them too
            del raw_filtered["query_id"]
            del raw_filtered["source"]
            del raw_filtered["is_day"]
            del raw_filtered["cloud_cover_low"]

            zf.writestr(
                f"{label}_filtered-data.csv",
                raw_filtered.to_csv(index=False)
            )

            # SIM_PARAMETER_FULL_DATA
            # Reshape into the simulation-parameter layout: one speed, stdev and
            # direction column per altitude
            df_q2 = pd.DataFrame()

            df_q2["date"] = df_q["date"]
            df_q2["temperature"] = df_q["temperature_2m"]
            df_q2["pressure"] = df_q["pressure_msl"]

            csv_size = len(df_q["date"])

            # Relies on ALTITUDES, WIND_SPEEDS and WIND_DIRECTIONS being in the same order
            for altitude, wind_speed, wind_direction in zip(config.ALTITUDES, config.WIND_SPEEDS, config.WIND_DIRECTIONS):
                df_q2[f"{altitude}"] = df_q[wind_speed]
                df_q2[f"stdev [{altitude}]"] = np.zeros(csv_size)  # placeholder; no std-dev data yet
                df_q2[f"direction [{altitude}]"] = df_q[wind_direction]

            # Kept only so make_filtered_sorted can filter on them below
            df_q2["cloud_cover_low"] = df_q["cloud_cover_low"]
            df_q2["is_day"] = df_q["is_day"]

            zf.writestr(
                f"sim_parameters_{label}_full-data.csv",
                df_q2.to_csv(index=False)
            )

            # SIM_PARAMETER_FILTERED_DATA
            sim_filtered = make_filtered_sorted(df_q2)

            del sim_filtered["cloud_cover_low"]
            del sim_filtered["is_day"]

            zf.writestr(
                f"sim_parameters_{label}_filtered-data.csv",
                sim_filtered.to_csv(index=False)
            )

            # Copy for the on-screen preview, with the timezone stripped for display
            last_sim_parameters_filtered = sim_filtered.copy()

            last_sim_parameters_filtered["date"] = (
                last_sim_parameters_filtered["date"]
                .dt.tz_localize(None)
                .dt.strftime("%Y-%m-%d %H:%M")
            )

    # Rewind so getvalue() returns the whole archive
    zip_buffer.seek(0)

    st.download_button(
        "Download all files (ZIP)",
        data=zip_buffer.getvalue(),
        file_name=f"open-meteo-download_{today}.zip",
        mime="application/zip",
    )

    # raw_filtered and last_sim_parameters_filtered still hold the last loop iteration
    st.subheader("Raw Filtered Preview (last query)")
    st.dataframe(
        raw_filtered.head(200),
        use_container_width=True
    )

    st.subheader("Sim Parameters Filtered Preview (last query)")
    st.dataframe(
        last_sim_parameters_filtered.head(200),
        use_container_width=True
    )

    st.markdown(
    """
    <p style="font-size:14px;">
        Note: Displaying two (filtered) of the four files (raw + filtered)
    </p>
    """,
    unsafe_allow_html=True
    )
