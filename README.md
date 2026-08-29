# Project Atlas

Project Atlas is Waterloo Rocketry's toolkit for visualizing and analyzing sounding rocket dispersion data — essentially, a way of answering "where is the rocket going to land, and how confident is that answer." It started as a way to plot OpenRocket Monte Carlo exports on a map, and has since grown into a small suite of tools covering dispersion analysis, flight dynamics (damping/coning), and a handful of geography and weather utilities that feed into the main app.

The core deliverable is a desktop app that takes Monte Carlo simulation CSVs from OpenRocket and plots the landing zone dispersion on a real map, complete with confidence ellipses, water landing risk, and a stats panel. Everything else in the repo either feeds that app or performs standalone analysis on the side.

---

## Table of Contents

- [Project Breakdown](#project-breakdown)
- [The Main App](#the-main-app)
  - [Main App Functionality](#main-app-functionality)
  - [Main App Usage](#main-app-usage)
  - [Expected Input Format](#expected-input-format)
- [Supporting Scripts](#supporting-scripts)
  - [Flight Dynamics Analysis](#flight-dynamics-analysis)
  - [Geography & Weather Utilities](#geography--weather-utilities)
- [Setup](#setup)
- [Repo Structure](#repo-structure)
- [Known Quirks](#known-quirks)
- [Contributing](#contributing)

---

## Project Breakdown

| File | What It Does |
|---|---|
| `main.py` | Entry point — launches the dispersion analysis GUI |
| `app.py` | The PySide6 (Qt) GUI itself: buttons, file lists, stats panel, plot window |
| `data_engine.py` | Backend data analysis: reading CSVs, computing stats, outlier detection, water landing checks |
| `plotting.py` | Matplotlib rendering: the actual dispersion plot, ellipses, basemap, styling |
| `calculate_dampingratio.py` | Standalone script for pitch/yaw damping ratio over a flight |
| `coning_analysis.py` | Standalone script for natural frequency vs. roll rate (coning risk) |
| `download_lakes.py` | One-off script to pull lake/water body geometry for water-landing detection |
| `download_map.py` | One-off script to pull a satellite basemap tile for offline use |
| `extract_roads.py` | Converts a Google Earth KML of nearby highways into a CSV for plotting |
| `generate_montecarlo.py` | Interactive CLI tool to build randomized wind-dispersion input files for OpenRocket |
| `scrape_navcan.py` | Scrapes NavCanada's upper wind forecasts for use in the above |

The four files that do the heavy lifting for the actual dispersion tool are `main.py`, `app.py`, `data_engine.py`, and `plotting.py`. Everything else is a smaller, mostly standalone utility — some are meant to be run once to generate a static file the main app depends on (lakes, basemap, roads), others are independent analysis scripts for other parts of vehicle design.

---

## The Main App

### Main App Functionality

When `main.py` is run, a window opens that:

1. Allows one or more OpenRocket Monte Carlo export CSVs (landing dispersion data) to be selected
2. Plots every simulated landing point on a map centered on the launch site, with:
   - An optional 10 NM "waiver" radius ellipse (Launch Canada Advanced Pad accuracy requirement)
   - Optional 1σ/2σ/3σ dispersion ellipses
   - An optional statistical confidence ellipse (confidence level is user-set, e.g. 0.95)
   - Water landings highlighted in cyan (checked against a local lakes dataset)
   - The top N outlier landings highlighted in yellow
   - Local geography for context: highways, a nearby lodge, launch site marker
3. Produces a "Flight Statistics" panel per file, showing mean apogee, landing distance, accuracy %, min stability, lateral velocity, wind speed, etc.
4. Supports drilling into wind-shear outlier analysis for a given simulation, provided a matching sim-parameters file has been uploaded
5. Allows the plot to be saved (PNG/JPG/MP4) or the stats panel to be exported to CSV

CSV handling is run on a background thread so the UI isn't frozen on large files, and files are cached in memory so re-plotting or exporting doesn't trigger a re-parse from disk.

### Main App Usage

```bash
python main.py
```

No CLI arguments are required. The last folder files were selected from is remembered automatically (via `QSettings`), so it doesn't need to be re-navigated to on the next run.

### Expected Input Format

OpenRocket Monte Carlo export CSVs are expected. Column-matching is done by substring, so exact column order doesn't matter, but the following are searched for:

- `Apogee`
- `Landing Latitude` / `Landing Longitude` (both `(°N)`/`(°E)` and `(deg N)`/`(deg E)` naming conventions are handled)
- `Min Stability`
- `Lateral Velocity at Apogee`
- `Max Windspeed` (assumed to already be in knots)

If a column can't be found, a clear error is raised describing what was expected and what columns actually exist in the file, rather than a cryptic `StopIteration` or `KeyError`.

---

## Supporting Scripts

These aren't wired into the GUI — they're standalone scripts, each with its own config block near the top (usually marked `[MODIFY ACCORDINGLY]`) and a usage docstring describing what's needed. The general pattern: open the script, adjust the paths/values, run it.

### Flight Dynamics Analysis

**`calculate_dampingratio.py`** — Computes pitch/yaw damping ratio across the flight, combining jet damping (from motor mass flow rate) and aerodynamic damping (from per-component CNa/CP data). An OpenRocket flight export plus two Component Analysis exports (CP vs. Mach and CNa vs. Mach, per-component) are required as inputs. A static PNG or an animated MP4 is produced depending on the config flags.

**`coning_analysis.py`** — Computes natural frequency (ωₙ) across the flight and compares it to roll rate, flagging roll-pitch coupling ("coning") risk. A 300-sample Monte Carlo dispersion is run on mass properties, aero coefficients, and atmospheric conditions to produce a 95% confidence band, with roll rate curves at different fin cant angles overlaid. Static PNG or animated MP4 output, same as above.

Both scripts are self-contained — the main app doesn't need to be running, just the right OpenRocket export CSVs sitting in a `sim files/` folder.

### Geography & Weather Utilities

These generate static assets that `plotting.py` (and therefore the main app) depends on. They're generally only needed once, or when the launch site/area of interest changes — pulling the repo from GitHub should already include the generated files, so these are best treated as "regenerate if missing or outdated" tools rather than routine steps.

**`download_lakes.py`** — Pulls lake/water body geometry from OpenStreetMap (via `osmnx`) within a 20 km radius of the launch site, saved as `LC Geography/lakes.geojson`. This is what powers the water-landing detection in the main app.

**`download_map.py`** — Downloads a satellite basemap tile (ESRI World Imagery, via `contextily`) for the area around the launch site, saved as `LC Geography/lc_basemap.tif`. This exists because GitHub doesn't handle large binary files gracefully, so the basemap can be regenerated locally instead of relying on it being tracked in the repo.

**`extract_roads.py`** — Takes a KML file (exportable from Google Earth) containing nearby highways and converts it to a CSV of lat/lon points, which `plotting.py` overlays on the dispersion map for context.

**`generate_montecarlo.py`** — Interactive CLI tool that takes a NavCanada wind export (see below) and generates a randomized set of wind/temperature/pressure conditions, formatted for OpenRocket's Monte Carlo plugin. Standard deviations per altitude layer (or a single blanket value), temperature, and pressure are configured interactively, followed by the number of simulations to generate.

**`scrape_navcan.py`** — Scrapes NavCanada's upper wind forecast tool for Kapuskasing Airport (CYYU), the closest station to the launch site, and exports the 7/9/12-hour wind prediction windows to CSV. Output feeds directly into `generate_montecarlo.py`. Playwright is used to drive a headless browser, since the data isn't available through a public API.

---

## Setup

Python 3.10+ is required, and a virtual environment is recommended. Not every script needs every dependency — installation can be scoped to whatever's actually being run — but for the full toolkit:

```bash
pip install pandas numpy scipy matplotlib PySide6 geopandas shapely osmnx contextily playwright
playwright install chromium   # only needed for scrape_navcan.py
```

Notes:
- The core dispersion tool (`main.py` / `app.py` / `data_engine.py` / `plotting.py`) requires `pandas`, `numpy`, `scipy`, `matplotlib`, `PySide6`, `geopandas`, and `shapely`.
- `download_lakes.py` requires `osmnx`.
- `download_map.py` requires `contextily`.
- `scrape_navcan.py` requires `playwright` (plus the browser binary — see above).
- `generate_montecarlo.py` uses `tkinter`, which ships with most Python installs but occasionally needs a separate OS package on Linux (`python3-tk`).

---

## Repo Structure

The project is roughly organized as such:

```
project-root/
├── Dispersion Analysis/        <- main.py, app.py, data_engine.py, plotting.py live here
├── LC Geography/                <- lakes.geojson, lc_basemap.tif, highway_144+101.csv
├── Road Coordinates/            <- highway_144+101.kml (source for extract_roads.py)
└── sim files/                   <- OpenRocket exports for damping/coning scripts
```

`plotting.py` locates `LC Geography/` relative to its own file location (two directories up), so as long as the folder structure above is preserved, paths should resolve correctly without edits. The flight-dynamics scripts and `extract_roads.py` use relative paths from the working directory at runtime, so care should be taken there. These paths should be identified with the `[MODIFY ACCORDINGLY]` paths at the top of each file can simply be edited directly.

---

## Known Quirks

- **Launch site coordinates are hardcoded** in `data_engine.py` (`LAUNCH_LAT`, `LAUNCH_LON`) and duplicated in `download_lakes.py`/`download_map.py`. If the launch site ever changes, all of these need to be updated.
- **Wind speed units in knots** throughout (main app stats panel, outlier analysis, `generate_montecarlo.py`). No conversion is applied, so input files pulled from a source using different units should be checked carefully.
- **`haversine_nm` returns nautical miles**, not statute miles; any lateral distances are treated in nautical miles.
- The core app functionality is designed to accommodate two-stage vehicles as well; however some development work is still necessary both in terms of data analysis and visual display settings. 
- **`scrape_navcan.py`** depends on NavCanada's site structure (selectors, layout) remaining unchanged. If it breaks, a site update is the most likely cause — the CSS selectors and toggle IDs are the first things worth checking.
- The animated MP4 outputs in `calculate_dampingratio.py` and `coning_analysis.py` require `ffmpeg` to be installed and available on the system PATH.

---

## Contributing

- Since `app.py` imports specific names from both `data_engine.py` and `plotting.py`, those import lists are worth checking before anything is renamed or removed in either file.
- Docstrings at the top of each file describe expected usage and should be kept up to date if input formats change, since these scripts are mostly run manually, as needed, by whoever on the team requires them.
