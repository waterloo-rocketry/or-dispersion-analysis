import os
import re
import numpy as np
import pandas as pd
import tkinter as tk
import contextily as cx
import matplotlib as mpl
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

from pathlib import Path
from scipy.stats import chi2
from tkinter import filedialog


def _get_col(data, *names):
    """
    Returns the first matching column from a DataFrame by trying multiple possible column names.
    Avoids the ambiguous truth-value error that occurs when using `or` with pandas Series.
    :param data:    pandas DataFrame
    :param names:   column names to try in order
    :return:        first matching Series, or None if none found
    """
    for name in names:
        col = data.get(name)
        if col is not None:
            return col
    return None


def print_stats(name, series):
    """
    Function to print statistics to console
    :param name:    .csv filename
    :param series:  pandas series of data
    :return:        console output; average, standard deviation, min and max values of data series
    """
    if series is not None:
        print(f"{name}")
        print(f">> Average: {round(series.mean(), 3)}")
        print(f">> Std Dev: {round(series.std(), 3)}")
        print(f">> Min: {round(series.min(), 3)}")
        print(f">> Max: {round(series.max(), 3)}")
    else:
        print(f"*** Series {name} not found: TypeError ***")


def get_file_paths(initial_directory):
    """
    Function enabling user-selected files using filedialog.
    :param initial_directory:   custom initial directory so that filedialog does not open at machine root directory
    :return:                    array of selected .csv filenames
    """
    root = tk.Tk()
    root.withdraw()

    file_paths = filedialog.askopenfilenames(
        title="Select files",
        initialdir=initial_directory,
        filetypes=[("CSV files", "*.csv"), ("All files", "*.*")]
    )

    file_paths = list(file_paths)
    for f in file_paths:
        print(f)

    root.destroy()

    if not file_paths:
        print(f"No file selected")
        exit(0)

    return file_paths


def generate_labels(file_paths):
    """
    Auto-generates label names for each .csv file.
    :param file_paths:  array of .csv file_paths
    :return:            array of label names
    """
    labels = []

    for path in file_paths:
        filename = os.path.basename(path)
        name, _ = os.path.splitext(filename)
        name = name.replace("_", " ").replace("-", " ")
        label = name.title()
        labels.append(label)

    print(labels)

    return labels


def get_colors(file_paths):
    """
    Chooses evenly spaced colors for plotting several different .csv files; ensures sufficient contrast between white
    background and data points.
    :param file_paths:  array of .csv file_paths
    :return:            array of evenly spaced colors; array of defined colors for plotting sigma ellipses
    """
    colors = cm.viridis(np.linspace(0, 1, len(file_paths) + 1))  # add +1 to avoid pale colors

    sigma_colors = [
        "xkcd:red wine",
        "xkcd:cherry red",
        "xkcd:orangish"
    ]

    return colors, sigma_colors


def make_safe_filename(title, ext="png"):
    """
    Removes any special characters from plot title to prevent errors when saving file.
    :param title:   user-provided plot title
    :param ext:     specifies file extension
    :return:        file path as type Path
    """
    if not title:
        title = "plot"
    safe = re.sub(r'[^\w\s-]', '', title)
    safe = re.sub(r'\s+', '_', safe)
    return Path(f"{safe}.{ext}")


def set_default_style():
    """
    Sets default plotting style.
    :return:
    """
    mpl.rcParams.update({
        "text.usetex": False,
        "mathtext.fontset": "cm",
        "font.family": "serif",
        "font.size": 5,
        "axes.labelsize": 9,
        "axes.titlesize": 10,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
        "legend.fontsize": 8,
        "axes.linewidth": 0.8,
        "grid.linewidth": 0.5,
        "lines.linewidth": 1.0,
    })


# ── Ellipse math ──────────────────────────────────────────────────────────────

def ellipse_math(x, y):
    """
    Function to compute properties of the ellipse required to plot it.
    :param x:   x values of population
    :param y:   y values of population
    :return:    center (x,y) coordinates of ellipse; semi-major and semi-minor axis length; rotation of ellipse
    """
    mean_x, mean_y =    np.mean(x), np.mean(y)
    cov =               np.cov(x, y)
    vals, vecs =        np.linalg.eigh(cov)
    order =             vals.argsort()[::-1]
    vals =              vals[order]
    vecs =              vecs[:, order]
    theta =             np.degrees(np.arctan2(*vecs[:, 0][::-1]))

    return mean_x, mean_y, vals, theta


def plot_ellipse(x, y, major, minor, tilt, edge_color, style, name, ax):
    """
    Plots ellipses with provided arguments.
    """
    ellipse = mpatches.Ellipse(
        (x, y),
        width=major,
        height=minor,
        angle=tilt,
        edgecolor=edge_color,
        facecolor="none",
        linestyle=style,
        linewidth=1.5,
        label=name,
        alpha=1,
    )
    ax.add_patch(ellipse)


# ── Known locations ───────────────────────────────────────────────────────────

def plot_known_locations(launch_lon, launch_lat, ax):
    """
    Function to plot known launch and obstacle locations around Launch Canada.
    """
    lodge_lat       = 47.9017898
    lodge_lon       = -81.6497071

    highway_coords  = pd.read_csv("../LC Geography/highway_144+101.csv")
    highway_lat     = highway_coords["Latitude"]
    highway_lon     = highway_coords["Longitude"]

    ax.plot(highway_lon, highway_lat, color="grey", linewidth=2, label="Hwy 144 & Hwy 101")
    ax.scatter(launch_lon, launch_lat, edgecolors="black", color="yellow", marker="*", s=150, label="Launch Site")
    ax.scatter(lodge_lon, lodge_lat, edgecolors="black", color="xkcd:lightblue", marker="^", s=80, label="Tata Chika Pika Lake Lodge")


# ── Haversine / outliers ──────────────────────────────────────────────────────

def haversine_nm(ref_latitude, ref_longitude, latitude, longitude):
    """
    Compute the great-circle distance between two points on Earth (in nautical miles).
    """
    R_nm = 3443.92

    ref_latitude, ref_longitude, latitude, longitude = map(np.radians, [ref_latitude, ref_longitude, latitude, longitude])

    dlat = latitude - ref_latitude
    dlon = longitude - ref_longitude

    a = np.sin(dlat / 2) ** 2 + np.cos(ref_latitude) * np.cos(latitude) * np.sin(dlon / 2) ** 2
    c = 2 * np.arctan2(np.sqrt(a), np.sqrt(1 - a))

    return R_nm * c


def get_outliers(radius_nm, ref_latitude, ref_longitude, latitude, longitude):
    """
    Returns indices of points outside the given radius.
    """
    distance_nm = haversine_nm(ref_latitude, ref_longitude, latitude, longitude)
    print(distance_nm)

    outside_pts = distance_nm > radius_nm

    if hasattr(latitude, "index"):
        indices = latitude.index[outside_pts]
    else:
        indices = np.where(outside_pts)[0]

    return indices


# ── Shared extract helper ─────────────────────────────────────────────────────

def _extract_columns(data):
    """
    Extracts rocket and payload lat/lon columns from a DataFrame, trying
    both Polaris and Aurora column name conventions.
    :param data:    pandas DataFrame read from a .csv file
    :return:        tuple of (rocket_latitude, rocket_longitude, payload_latitude, payload_longitude)
                    any of which may be None if not found
    """

    lat_suffixes = ("Landing Latitude (°N)", "Landing Latitude (deg N)")
    lon_suffixes = ("Landing Longitude (°E)", "Landing Longitude (deg E)")

    lat_col = next((col for col in data.columns if col.endswith(lat_suffixes)), None)
    lon_col = next((col for col in data.columns if col.endswith(lon_suffixes)), None)

    rocket_latitude = data[lat_col] if lat_col else None
    rocket_longitude = data[lon_col] if lon_col else None

    payload_latitude = None
    payload_longitude = None

    # rocket_latitude = _get_col(
    #     data,
    #     "Polaris Rocket Landing Latitude (°N)",
    #     "Aurora Rocket Landing Latitude (°N)"
    # )
    # rocket_longitude = _get_col(
    #     data,
    #     "Polaris Rocket Landing Longitude (°E)",
    #     "Aurora Rocket Landing Longitude (°E)"
    # )
    # payload_latitude = _get_col(
    #     data,
    #     "Polaris Deployable Payload Landing Latitude (°N)",
    #     "Deployable Payload Landing Latitude (°N)"
    # )
    # payload_longitude = _get_col(
    #     data,
    #     "Polaris Deployable Payload Landing Longitude (°E)",
    #     "Deployable Payload Landing Longitude (°E)"
    # )

    return rocket_latitude, rocket_longitude, payload_latitude, payload_longitude


# ── Shared plotting function ──────────────────────────────────────────────

def draw_plot_elements(ax, file_paths, plot_title, plot_LC_ellipse, plot_sigma_ellipses, plot_confidence_ellipse, confidence, plot_top_outliers=True):
    """
    Core logic to populate a matplotlib Axes object with the dispersion data.
    Plots all files: rocket and payload scatter points plus optional ellipses onto ax.
    :param ax:                      matplotlib axes object where data is plotted
    :param file_paths:              array of .csv file file_path names with longitude/latitude data
    :param plot_title:              user-provided title of plot
    :param plot_LC_ellipse:         bool, user-option to plot arbitrary 10nm radius around advanced pad launch site
    :param plot_sigma_ellipses:     bool, user-option to plot ellipses enclosing n * standard deviation
    :param plot_confidence_ellipse: bool, user-option to plot ellipse enclosing n% of data points
    :param confidence:              float, user-option specifying perecentage of data points to enclose
    :return:                        plots data in dispersion_zone app
    """
    scatter_labels = generate_labels(file_paths)
    colors, sigma_colors = get_colors(file_paths)

    launch_lat = 47.965378
    launch_lon = -81.873536

    # Loop to read files
    for i, (file_path, label, color) in enumerate(zip(file_paths, scatter_labels, colors)):
        if not os.path.exists(file_path):
            print(f"File not found: {file_path}")
            return

        data = pd.read_csv(file_path, encoding="utf-8", sep=None, engine="python", index_col=False)

        rocket_latitude, rocket_longitude, payload_latitude, payload_longitude = _extract_columns(data)

        if rocket_longitude is not None:
            ax.scatter(
                rocket_longitude,
                rocket_latitude,
                facecolor=color,
                edgecolor=color,
                marker=".",
                s=10,
                alpha=0.7,
                label=f"{label}"
            )

            if plot_top_outliers:
                # Calculate distances of all points from the launch pad
                distances = haversine_nm(launch_lat, launch_lon, rocket_latitude, rocket_longitude)
                top_50 = distances.nlargest(100)

                # Hardcode the exact string name of the column
                sim_col = 'Simulation'

                # Dictionary to track the frequency of specific dates
                outlier_dates = {}

                for idx in top_50.index:
                    # np.atleast_1d safeguards against duplicate indices returning a Series
                    pt_lat = float(np.atleast_1d(rocket_latitude.loc[idx])[0])
                    pt_lon = float(np.atleast_1d(rocket_longitude.loc[idx])[0])

                    # Fetch sim name safely using the string column name
                    raw_val = data[sim_col].loc[idx]
                    sim_name = str(np.atleast_1d(raw_val)[0])

                    # Extract just the date by splitting the string at the space
                    # "2021-08-12 17:00:00+00:00" -> "2021-08-12"
                    sim_date = sim_name.split(" ")[0]

                    # Track duplicates by incrementing the count for this date
                    if sim_date in outlier_dates:
                        outlier_dates[sim_date] += 1
                    else:
                        outlier_dates[sim_date] = 1

                    # Draw highlight circle
                    ax.scatter(pt_lon, pt_lat, facecolors="none", edgecolors="xkcd:dandelion", s=60, linewidths=1.5,
                               zorder=5)

                # Print the summary to the console after the loop finishes
                print(f"\n--- Outlier Dates Summary ({label}) ---")
                for date, count in sorted(outlier_dates.items()):
                    print(f"Date: {date} | Outliers: {count}")
                print("-" * 40 + "\n")

        if payload_longitude is not None:
            ax.scatter(
                payload_longitude,
                payload_latitude,
                color=color,
                marker="x",
                s=10,
                alpha=0.7,
                label=f"{label} - Payload"
            )

        if plot_sigma_ellipses and rocket_longitude is not None:
            mean_x, mean_y, vals, theta = ellipse_math(rocket_longitude.values, rocket_latitude.values)

            sigma_levels =   [1, 2]
            ellipse_params = [(np.sqrt(vals[0]) * level * 2, np.sqrt(vals[1]) * level * 2) for level in sigma_levels]

            for (width, height), level, sigma_color in zip(ellipse_params, sigma_levels, sigma_colors):
                plot_ellipse(
                    x=mean_x,
                    y=mean_y,
                    major=width,
                    minor=height,
                    tilt=theta,
                    edge_color=sigma_color,
                    style="-",
                    # name=f"{level}σ ellipse for {label}",
                    name=f"",
                    ax=ax
                )

        if plot_confidence_ellipse and rocket_longitude is not None:
            mean_x, mean_y, vals, theta = ellipse_math(rocket_longitude.values, rocket_latitude.values)

            chi2_val   = chi2.ppf(confidence, 2)
            major_axis = 2 * np.sqrt(vals[0] * chi2_val)
            minor_axis = 2 * np.sqrt(vals[1] * chi2_val)

            plot_ellipse(
                x=mean_x,
                y=mean_y,
                major=major_axis,
                minor=minor_axis,
                tilt=theta,
                edge_color="xkcd:bluish purple",
                style="-",
                name=f"{confidence * 100}% Confidence Ellipse",
                ax=ax
            )

    if plot_LC_ellipse:
        diameter_deg_lat = 2 * (18520 / 111320)
        diameter_deg_lon = 2 * (18520 / (111320 * np.cos(np.radians(launch_lat))))

        plot_ellipse(
            x=launch_lon,
            y=launch_lat,
            major=diameter_deg_lon,
            minor=diameter_deg_lat,
            tilt=0,
            edge_color="red",
            style="--",
            name="10 NM Radius",
            ax=ax
        )

    # 1x1 Mile Grid Math
    lat_ref = 47.704847
    long_ref = -82.510814
    width, length = 47, 35
    deg_lat_per_mile = 1 / 69.0
    deg_lon_per_mile = 1 / (69.17 * np.cos(np.radians(lat_ref)))
    width_deg, height_deg = width * deg_lon_per_mile, length * deg_lat_per_mile

    # Waiver Rectangle
    rectangle = mpl.patches.Rectangle(
        xy=(long_ref, lat_ref),
        width=width_deg,
        height=height_deg,
        edgecolor="red",
        linestyle="-",
        linewidth=2,
        facecolor="none",
        label="LC Full Waiver Boundary",
        zorder=3
    )
    ax.add_patch(rectangle)

    # Grid Lines
    x_grid = [long_ref + (i * deg_lon_per_mile) for i in range(1, width)]
    y_grid = [lat_ref + (i * deg_lat_per_mile) for i in range(1, length)]
    ax.vlines(x=x_grid, ymin=lat_ref, ymax=lat_ref + height_deg, colors='gray', linestyles=':', linewidth=0.5, zorder=2)
    ax.hlines(y=y_grid, xmin=long_ref, xmax=long_ref + width_deg, colors='gray', linestyles=':', linewidth=0.5, zorder=2)

    buffer = 0.1
    ax.set_xlim(long_ref - buffer, long_ref + width_deg + buffer)
    ax.set_ylim(lat_ref - buffer, lat_ref + height_deg + buffer)

    ax.margins(0.03)
    ax.set_aspect('equal', adjustable='datalim')
    ax.figure.canvas.draw()

    plot_known_locations(launch_lon, launch_lat, ax)

    try:
        cx.add_basemap(
            ax,
            crs="EPSG:4326",
            source="LC Geography/lc_basemap.tif",
            zorder=0,
            alpha=0.9,
            attribution=False
        )
    except Exception as e:
        print(f"Could not load local basemap: {e}")
        cx.add_basemap(
            ax,
            crs="EPSG:4326",
            source=cx.providers.Esri.WorldImagery,
            zorder=0,
            alpha=0.9,
            attribution=False
        )

    handles, labels = ax.get_legend_handles_labels()
    if handles:
        sorted_pairs = sorted(zip(labels, handles), key=lambda x: 0 if x[0] in scatter_labels else 1)
        sorted_labels, sorted_handles = zip(*sorted_pairs)
        ax.legend(sorted_handles, sorted_labels)

    ax.set_xlabel("Longitude (°)")
    ax.set_ylabel("Latitude (°)")
    ax.set_title(plot_title)

# ── Public plotting functions ─────────────────────────────────────────────────

def plot_data(file_paths, plot_title, fig, ax, plot_LC_ellipse=True, plot_sigma_ellipses=False, plot_confidence_ellipse=False, confidence=0.95, plot_top_outliers=True):
    """Updates the Tkinter UI plot."""
    ax.clear()
    draw_plot_elements(ax, file_paths, plot_title, plot_LC_ellipse, plot_sigma_ellipses, plot_confidence_ellipse, confidence, plot_top_outliers)


def save_plot(file_paths, plot_title, output_path: str | None = None, plot_LC_ellipse=True, plot_sigma_ellipses=False, plot_confidence_ellipse=False, confidence=0.95, plot_top_outliers=True):
    """Saves a high-resolution plot bypassing Tkinter display quirks."""
    set_default_style()
    dispersion_plot, axes = plt.subplots(figsize=(10, 8))
    draw_plot_elements(axes, file_paths, plot_title, plot_LC_ellipse, plot_sigma_ellipses, plot_confidence_ellipse, confidence, plot_top_outliers)
    dispersion_plot.tight_layout()

    # Choose filename
    if output_path:
        filename = Path(output_path)
    else:
        filename = make_safe_filename(plot_title)

    filename.parent.mkdir(parents=True, exist_ok=True)
    dispersion_plot.savefig(filename, transparent=False, dpi=300)
    plt.close(dispersion_plot)
