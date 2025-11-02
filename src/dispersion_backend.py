import os
import re
import numpy as np
import pandas as pd
import tkinter as tk
import matplotlib as mpl
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

from pathlib import Path
from scipy.stats import chi2
from tkinter import filedialog


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

    # Clean up file name
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


def ellipse_math(x, y):
    """
    Function to compute properties of the ellipse required to plot it:
        1.  Computes center (x,y) coordinates of ellipse
        2.  Computes covariance matrix of population
            --> Gives the variance of x values w.r.t. y values and vice versa
        3.  Computes eigenvalues and eigenvectors of population
            --> Eigenvalues correspond to semi-major and semi-minor axes
            --> Eigenvectors define rotation of the ellipse
        4.  Computes rotation of the ellipse using eigenvectors
    :param x:   x values of population
    :param y:   y values of population
    :return:    center (x,y) coordinates of ellipse; semi-major and semi-minor axis length; rotation of ellipse
    """
    mean_x, mean_y =    np.mean(x), np.mean(y)
    cov =               np.cov(x, y)
    vals, vecs =        np.linalg.eigh(cov)     # eigen decomposition
    order =             vals.argsort()[::-1]    # sort largest eigenvalue first
    vals =              vals[order]
    vecs =              vecs[:, order]
    theta =             np.degrees(np.arctan2(*vecs[:, 0][::-1]))  # Angle of ellipse in degrees

    return mean_x, mean_y, vals, theta


def plot_ellipse(x, y, major, minor, tilt, edge_color, style, name, ax):
    """
    Plots ellipses with provided arguments.
    :param x:
    :param y:
    :param major:       2 * semi-major axis
    :param minor:       2 * semi-minor axis
    :param tilt:        rotation of the ellipse relative to horizontal; in radians
    :param edge_color:
    :param style:
    :param name:
    :param ax:
    :return:            plots ellipse
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


def plot_known_locations(launch_lon, launch_lat, ax):
    """
    Function to plot known launch and obstacle locations around Launch Canada
    :param launch_lon:  advanced pad launch longitude
    :param launch_lat:  advanced pad launch latitude
    :param ax:          matplotlib axes object where data is plotted
    :return:            plots known locations on ax object
    """
    # Known Coordinates
    lodge_lat =         47.9017898
    lodge_lon =         -81.6497071

    highway_144 =       pd.read_csv("../Road Coordinates/highway_144.csv")
    highway_144_lat =   highway_144["Latitude"]
    highway_144_lon =   highway_144["Longitude"]

    ax.plot(highway_144_lon, highway_144_lat, color="grey", linewidth=2, label="Highway 144")
    ax.scatter(launch_lon, launch_lat, edgecolors="black", color="yellow", marker="*", s=150, label="Launch Site")
    ax.scatter(lodge_lon, lodge_lat, edgecolors="black", color="xkcd:lightblue", marker="^", s=80, label="Tata Chika Pika Lake Lodge")


def haversine_nm(ref_latitude, ref_longitude, latitude, longitude):
    """
    Compute the great-circle distance between two points on Earth (in nautical miles).
    All inputs and outputs are in degrees and nautical miles respectively.
    :param ref_latitude:
    :param ref_longitude:
    :param latitude:
    :param longitude:
    :return:
    """
    R_nm = 3443.92  # Earth's radius in nautical miles

    # Convert all angles to radians
    ref_latitude, ref_longitude, latitude, longitude = map(np.radians, [ref_latitude, ref_longitude, latitude, longitude])

    dlat = latitude - ref_latitude
    dlon = longitude - ref_longitude

    a = np.sin(dlat / 2) ** 2 + np.cos(ref_latitude) * np.cos(latitude) * np.sin(dlon / 2) ** 2
    c = 2 * np.arctan2(np.sqrt(a), np.sqrt(1 - a))

    return R_nm * c


def get_outliers(radius_nm, ref_latitude, ref_longitude, latitude, longitude):
    """
    Incomplete.
    :param radius_nm:
    :param ref_latitude:
    :param ref_longitude:
    :param latitude:
    :param longitude:
    :return:
    """
    distance_nm = haversine_nm(ref_latitude, ref_longitude, latitude, longitude)
    print(distance_nm)

    # Boolean mask for points outside the radius
    outside_pts = distance_nm > radius_nm

    # If latitude is a Series, use its index; otherwise, use NumPy where()
    if hasattr(latitude, "index"):
        indices = latitude.index[outside_pts]
    else:
        indices = np.where(outside_pts)[0]

    return indices


def plot_data(file_paths, plot_title, ax, plot_LC_ellipse=True, plot_sigma_ellipses=False, plot_confidence_ellipse=False, confidence=0.95):
    """
    Function to plot data from dispersion analysis plugin including options for statistical analysis of data points.
    :param file_paths:              array of .csv file file_path names with longitude/latitude data
    :param plot_title:              user-provided title of plot
    :param ax:                      matplotlib axes object where data is plotted
    :param plot_LC_ellipse:         bool, user-option to plot arbitrary 10nm radius around advanced pad launch site
    :param plot_sigma_ellipses:     bool, user-option to plot ellipses enclosing n * standard deviation
    :param plot_confidence_ellipse: bool, user-option to plot ellipse enclosing n% of data points
    :param confidence:              float, user-option specifying perecentage of data points to enclose
    :return:                        plots data in dispersion_zone app
    """
    scatter_labels =            generate_labels(file_paths)
    colors, sigma_colors =      get_colors(file_paths)

    ax.clear()

    launch_lat = 47.965378
    launch_lon = -81.873536

    # Loop to read files
    for i, (file_path, label, color) in enumerate(zip(file_paths, scatter_labels, colors)):

        if not os.path.exists(file_path):
            print(f"File not found: {file_path}")
            continue

        # Extract data
        data = pd.read_csv(file_path, encoding="utf-8", sep=None, engine="python", index_col=False)

        rocket_latitude         = data.get("Aurora Rocket Landing Latitude (°N)")
        rocket_longitude        = data.get("Aurora Rocket Landing Longitude (°E)")
        payload_latitude        = data.get("Deployable Payload Landing Latitude (°N)")
        payload_longitude       = data.get("Deployable Payload Landing Longitude (°E)")

        # Plot rocket
        if rocket_longitude is not None:
            ax.scatter(
                rocket_longitude,
                rocket_latitude,
                facecolor=color,
                edgecolor=color,
                marker=".",
                alpha=0.7,
                label=f"{label}"
            )

        # Plot payload
        if payload_longitude is not None:
            ax.scatter(
                payload_longitude,
                payload_latitude,
                color=color,
                marker="x",
                alpha=0.7,
                label=f"{label} - Payload"
            )

        if plot_sigma_ellipses:
            # Ellipse math: compute covariance matrix to get size and rotation of ellipse
            mean_x, mean_y, vals, theta = ellipse_math(rocket_longitude.values, rocket_latitude.values)

            # Ellipse radii for 1σ, 2σ, 3σ
            sigma_levels    = [1, 2, 3]
            ellipse_params  = [(np.sqrt(vals[0]) * level * 2, np.sqrt(vals[1]) * level * 2) for level in sigma_levels]

            # For each sigma level, plot the ellipse
            for (width, height), level, sigma_color in zip(ellipse_params, sigma_levels, sigma_colors):
                name = f"{level}σ ellipse for {label}"
                plot_ellipse(
                    x=mean_x,
                    y=mean_y,
                    major=width,
                    minor=height,
                    tilt=theta,
                    edge_color=sigma_color,
                    style="-",
                    name=name,
                    ax=ax
                )

        if plot_confidence_ellipse:
            # Ellipse math: compute covariance matrix to get size and rotation of ellipse
            mean_x, mean_y, vals, theta = ellipse_math(rocket_longitude.values, rocket_latitude.values)
            name        = f"{confidence * 100}% Confidence Ellipse"
            edge_color  = "xkcd:bluish purple"

            # Define size of ellipse based on confidence level
            chi2_val    = chi2.ppf(confidence, 2)
            major_axis  = (2 * np.sqrt(vals[0] * chi2_val))
            minor_axis  = (2 * np.sqrt(vals[1] * chi2_val))

            # 95% confidence ellipse
            plot_ellipse(
                x=mean_x,
                y=mean_y,
                major=major_axis,
                minor=minor_axis,
                tilt=theta,
                edge_color=edge_color,
                style="-",
                name=name,
                ax=ax)

    if plot_LC_ellipse:
        # Arbitrary Launch Canada limit on dispersion to prevent rockets from landing near highway 144
        # 1 NM = 1852 meters, so 10 NM = 18520 m; Convert meters to degrees approx (1 degree latitude ≈ 111.32 km)
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
            ax=ax)

    plot_known_locations(launch_lon, launch_lat, ax)

    # Tidy the labels and display in the legend
    handles, labels = ax.get_legend_handles_labels()
    if handles:
        sorted_pairs = sorted(zip(labels, handles), key=lambda x: 0 if x[0] in scatter_labels else 1)
        sorted_labels, sorted_handles = zip(*sorted_pairs)
        ax.legend(sorted_handles, sorted_labels)

    ax.set_xlabel("Longitude (°)")
    ax.set_ylabel("Latitude (°)")
    ax.set_title(plot_title)
    ax.grid(True)
    # ax.relim()
    # ax.autoscale_view()
    ax.margins(0.03)  # Add small margin so points/ellipses aren't on the plot edge
    ax.set_aspect('equal', adjustable='datalim')


def save_plot(file_paths, plot_title, output_path: str | None = None, plot_LC_ellipse=True, plot_sigma_ellipses=False, plot_confidence_ellipse=False, confidence=0.95):
    """
    Function to save user-generated plot. Addresses issue with high-def displays (specifically Mac machines) where plots
    generated in tkinter are blurry; gets around this issue by re-plotting in matplotlib.pyplot.
    :param file_paths:              array of .csv file file_path names with longitude/latitude data
    :param plot_title:              user-provided title of plot
    :param output_path:             file path and file name of plot to be saved
    :param plot_LC_ellipse:         bool, user-option to plot arbitrary 10nm radius around advanced pad launch site
    :param plot_sigma_ellipses:     bool, user-option to plot ellipses enclosing n * standard deviation
    :param plot_confidence_ellipse: bool, user-option to plot ellipse enclosing n% of data points
    :param confidence:              float, user-option specifying perecentage of data points to enclose
    :return:                        saves plot to output_path
    """
    scatter_labels =            generate_labels(file_paths)
    colors, sigma_colors =      get_colors(file_paths)

    set_default_style()
    dispersion_plot, axes = plt.subplots(figsize=(10, 8))

    launch_lat = 47.965378
    launch_lon = -81.873536

    # Loop to read files
    for i, (file_path, label, color) in enumerate(zip(file_paths, scatter_labels, colors)):

        if not os.path.exists(file_path):
            print(f"File not found: {file_path}")
            continue

        # Extract data
        data = pd.read_csv(file_path, encoding="utf-8", sep=None, engine="python", index_col=False)

        rocket_latitude     = data.get("Aurora Rocket Landing Latitude (°N)")
        rocket_longitude    = data.get("Aurora Rocket Landing Longitude (°E)")
        payload_latitude    = data.get("Deployable Payload Landing Latitude (°N)")
        payload_longitude   = data.get("Deployable Payload Landing Longitude (°E)")

        # Plot rocket
        if rocket_longitude is not None:
            axes.scatter(
                rocket_longitude,
                rocket_latitude,
                facecolor=color,
                edgecolor=color,
                marker=".",
                alpha=0.7,
                label=f"{label}"
            )

        # Plot payload
        if payload_longitude is not None:
            axes.scatter(
                payload_longitude,
                payload_latitude,
                color=color,
                marker="x",
                alpha=0.7,
                label=f"{label} - Payload"
            )

        if plot_sigma_ellipses:
            # Ellipse math: compute covariance matrix to get size and rotation of ellipse
            mean_x, mean_y, vals, theta = ellipse_math(rocket_longitude.values, rocket_latitude.values)

            # Ellipse radii for 1σ, 2σ, 3σ
            sigma_levels    = [1, 2, 3]
            ellipse_params  = [(np.sqrt(vals[0]) * level * 2, np.sqrt(vals[1]) * level * 2) for level in sigma_levels]

            for (width, height), level, sigma_color in zip(ellipse_params, sigma_levels, sigma_colors):
                name = f"{level}σ ellipse for {label}"
                plot_ellipse(
                    x=mean_x,
                    y=mean_y,
                    major=width,
                    minor=height,
                    tilt=theta,
                    edge_color=sigma_color,
                    style="-",
                    name=name,
                    ax=axes
                )

        if plot_confidence_ellipse:
            # Ellipse math: compute covariance matrix to get size and rotation of ellipse
            mean_x, mean_y, vals, theta = ellipse_math(rocket_longitude.values, rocket_latitude.values)
            name        = f"{confidence * 100}% Confidence Ellipse"
            edge_color  = "xkcd:bluish purple"

            # Define size of ellipse based on confidence level
            chi2_val    = chi2.ppf(confidence, 2)
            major_axis  = (2 * np.sqrt(vals[0] * chi2_val))
            minor_axis  = (2 * np.sqrt(vals[1] * chi2_val))

            # 95% confidence ellipse
            plot_ellipse(
                x=mean_x,
                y=mean_y,
                major=major_axis,
                minor=minor_axis,
                tilt=theta,
                edge_color=edge_color,
                style="-",
                name=name,
                ax=axes
            )

    if plot_LC_ellipse:
        # 1 NM = 1852 meters, so 10 NM = 18520 m; Convert meters to degrees approx (1 degree latitude ≈ 111.32 km)
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
            ax=axes
        )

    plot_known_locations(launch_lon, launch_lat, axes)

    # Tidy the labels and display in the legend
    handles, labels = axes.get_legend_handles_labels()
    if handles:
        sorted_pairs = sorted(zip(labels, handles), key=lambda x: 0 if x[0] in scatter_labels else 1)
        sorted_labels, sorted_handles = zip(*sorted_pairs)
        axes.legend(sorted_handles, sorted_labels)

    axes.set_xlabel("Longitude (°)")
    axes.set_ylabel("Latitude (°)")
    axes.set_title(plot_title)
    axes.grid(True)
    axes.margins(0.03)  # Add small margin so points/ellipses aren't on the plot edge
    axes.set_aspect('equal', adjustable='datalim')

    dispersion_plot.tight_layout()
    # plt.show()

    # Choose filename: user-provided output_path preferred; else crate a safe name
    if output_path:
        filename = Path(output_path)
    else:
        filename = make_safe_filename(plot_title)

    # Ensure parent dir exists (optional)
    filename.parent.mkdir(parents=True, exist_ok=True)

    # Save plot
    dispersion_plot.savefig(filename, transparent=False, dpi=300)


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
