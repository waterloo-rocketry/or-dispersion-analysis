import os
import re
import numpy as np
import pandas as pd
import contextily as cx
import matplotlib as mpl
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

from pathlib import Path
from scipy.stats import chi2

# Imports from data_engine.py
from data_engine import generate_labels, _extract_columns, ellipse_math, haversine_nm, ALTITUDES, _read_csv

TOP_OUTLIERS_COUNT  = 20
LC_GEOGRAPHY_DIR    = Path("../LC Geography")


def get_colors(file_paths):
    """
    Chooses evenly spaced colors for plotting several different .csv files; ensures sufficient contrast between white
    background and data points.
    :param file_paths:  array of .csv file_paths
    :return:            array of evenly spaced colors; array of defined colors for plotting sigma ellipses
    """
    n_files = len(file_paths)
    if n_files == 1:
        colors = [cm.viridis(0.65)]
    else:
        colors = cm.viridis(np.linspace(0.25, 1.0, n_files))

    sigma_colors = ["xkcd:red wine", "xkcd:sunflower", "xkcd:cherry red"]
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
    return Path(f"{re.sub(r'\s+', '_', safe)}.{ext}")


def set_default_style():
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


def plot_ellipse(x, y, major, minor, tilt, edge_color, style, name, ax):
    """
    Plots ellipses with provided arguments.
    """
    ellipse = mpatches.Ellipse(
        xy=(x, y),
        width=major,
        height=minor,
        angle=tilt,
        edgecolor=edge_color,
        facecolor="none",
        linestyle=style,
        linewidth=1.5,
        label=name,
        alpha=1
    )
    ax.add_patch(ellipse)


def plot_known_locations(launch_lon, launch_lat, ax):
    """
    Function to plot known launch and obstacle locations around Launch Canada.
    """
    lodge_lat = 47.9017898
    lodge_lon = -81.6497071
    try:
        highway_coords = pd.read_csv(LC_GEOGRAPHY_DIR / "highway_144+101.csv")
        ax.plot(
            highway_coords["Longitude"],
            highway_coords["Latitude"],
            color="grey",
            linewidth=2,
            label="Hwy 144 & Hwy 101"
        )
    except Exception as e:
        print(f"Highway layout missing: {e}")

    ax.scatter(
        launch_lon,
        launch_lat,
        edgecolors="black",
        color="yellow",
        marker="*",
        s=150,
        label="Launch Site"
    )
    ax.scatter(
        lodge_lon,
        lodge_lat,
        edgecolors="black",
        color="xkcd:lightblue",
        marker="^",
        s=80,
        label="Tata Chika Pika Lake Lodge"
    )


def draw_plot_elements(ax, file_paths, plot_title, plot_LC_ellipse, plot_sigma_ellipses, plot_confidence_ellipse,
                       confidence, plot_top_outliers=True):
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
    :param plot_top_outliers:       bool, user-option to plot/circle top outliers on the map
    :return:                        dict mapping each file_path to a {date: count} summary of its
                                    top-outlier landing dates (empty dict for files with no
                                    'Simulation' column, or where outlier highlighting was skipped)
    """
    scatter_labels = generate_labels(file_paths)
    colors, sigma_colors = get_colors(file_paths)
    launch_lat, launch_lon = 47.965378, -81.873536
    outlier_date_summary = {}

    for i, (file_path, label, color) in enumerate(zip(file_paths, scatter_labels, colors)):
        if not os.path.exists(file_path):
            continue
        data = _read_csv(file_path)
        rocket_lat, rocket_lon, payload_lat, payload_lon = _extract_columns(data)

        if rocket_lon is not None:
            ax.scatter(
                rocket_lon,
                rocket_lat,
                facecolor=color,
                edgecolor=color,
                marker=".",
                s=10,
                alpha=0.7,
                label=label
            )

            if plot_top_outliers:
                distances = haversine_nm(launch_lat, launch_lon, rocket_lat, rocket_lon)
                uh_ohs = distances.nlargest(TOP_OUTLIERS_COUNT)

                # Single vectorized scatter call instead of one ax.scatter() per point
                pts_lat = rocket_lat.loc[uh_ohs.index].values
                pts_lon = rocket_lon.loc[uh_ohs.index].values
                ax.scatter(
                    pts_lon,
                    pts_lat,
                    facecolors="none",
                    edgecolors="xkcd:dandelion",
                    s=60,
                    linewidths=1.5,
                    zorder=5
                )

                # Tally outlier counts by date, if a 'Simulation' column with date-prefixed
                # names is present (guarded so a missing column doesn't crash the plot).
                # Handed back to the caller instead of printed, so the UI can display it
                # (e.g. in the Outlier Analytics sidebar) rather than the console.
                outlier_dates = {}
                if 'Simulation' in data.columns:
                    for idx in uh_ohs.index:
                        sim_name = str(np.atleast_1d(data['Simulation'].loc[idx])[0])
                        sim_date = sim_name.split(" ")[0]
                        outlier_dates[sim_date] = outlier_dates.get(sim_date, 0) + 1

                outlier_date_summary[file_path] = outlier_dates

        if payload_lon is not None:
            ax.scatter(
                payload_lon,
                payload_lat,
                color=color,
                marker="x",
                s=10,
                alpha=0.7,
                label=f"{label} - Payload"
            )

        # Compute the ellipse geometry once and reuse it for both ellipse types,
        # rather than recomputing identical eigen-decomposition math twice.
        ellipse_stats = None
        if (plot_sigma_ellipses or plot_confidence_ellipse) and rocket_lon is not None:
            ellipse_stats = ellipse_math(rocket_lon.values, rocket_lat.values)

        if plot_sigma_ellipses and ellipse_stats is not None:
            mean_x, mean_y, vals, theta = ellipse_stats
            ax.scatter(
                mean_x,
                mean_y,
                edgecolors="black",
                color="xkcd:salmon",
                marker="P",
                s=80
            )
            for level, sigma_color in zip([1, 2], sigma_colors):
                plot_ellipse(
                    x=mean_x,
                    y=mean_y,
                    major=np.sqrt(vals[0]) * level * 2,
                    minor=np.sqrt(vals[1]) * level * 2,
                    tilt=theta,
                    edge_color=sigma_color,
                    style="-",
                    name="",
                    ax=ax
                )

        if plot_confidence_ellipse and ellipse_stats is not None:
            mean_x, mean_y, vals, theta = ellipse_stats
            chi2_val = chi2.ppf(confidence, 2)
            ax.scatter(
                mean_x,
                mean_y,
                edgecolors="black",
                color="xkcd:salmon",
                marker="P",
                s=80
            )
            plot_ellipse(
                x=mean_x,
                y=mean_y,
                major=2 * np.sqrt(vals[0] * chi2_val),
                minor=2 * np.sqrt(vals[1] * chi2_val),
                tilt=theta,
                edge_color="xkcd:bluish purple",
                style="-",
                name=f"{confidence * 100}% Confidence Ellipse",
                ax=ax
            )

    if plot_LC_ellipse:
        plot_ellipse(
            x=launch_lon,
            y=launch_lat,
            major=2 * (18520 / (111320 * np.cos(np.radians(launch_lat)))),
            minor=2 * (18520 / 111320),
            tilt=0,
            edge_color="red",
            style="--",
            name="10 NM Radius",
            ax=ax
        )

    lat_ref, long_ref = 47.704847, -82.510814
    width, length = 47, 35
    deg_lat_per_mile = 1 / 69.0
    deg_lon_per_mile = 1 / (69.17 * np.cos(np.radians(lat_ref)))
    width_deg, height_deg = width * deg_lon_per_mile, length * deg_lat_per_mile

    ax.add_patch(
        mpl.patches.Rectangle(
            xy=(long_ref, lat_ref),
            width=width_deg,
            height=height_deg,
            edgecolor="red",
            linestyle="-",
            linewidth=2,
            facecolor="none",
            label="LC Full Waiver",
            zorder=3
        )
    )

    ax.vlines(
        x=[long_ref + (i * deg_lon_per_mile) for i in range(1, width)],
        ymin=lat_ref,
        ymax=lat_ref + height_deg,
        colors='gray',
        linestyles=':',
        linewidth=0.5,
        zorder=2
    )
    ax.hlines(
        y=[lat_ref + (i * deg_lat_per_mile) for i in range(1, length)],
        xmin=long_ref,
        xmax=long_ref + width_deg,
        colors='gray',
        linestyles=':',
        linewidth=0.5,
        zorder=2
    )

    ax.set_aspect('equal', adjustable='datalim')
    plot_known_locations(launch_lon, launch_lat, ax)

    # `adjustable='datalim'` normally only gets applied lazily, at draw time -
    # matplotlib stretches xlim/ylim then to force a 1:1 aspect ratio. But
    # add_basemap() below reads ax.get_xlim()/get_ylim() *immediately* to decide
    # how much of the satellite image to crop. If that read happens before the
    # aspect ratio has actually been applied, it crops a smaller region than the
    # axes ends up occupying once rendered - leaving white space in the gap
    # between the (smaller) fetched image and the (larger) final view extent.
    # Calling apply_aspect() here forces that adjustment to happen now, so the
    # basemap is cropped to the true final extent instead.
    ax.apply_aspect()

    try:
        cx.add_basemap(
            ax,
            crs="EPSG:4326",
            source=str(LC_GEOGRAPHY_DIR / "lc_basemap.tif"),
            zorder=0,
            alpha=0.85,
            attribution=False
        )
    except Exception:
        # cx.add_basemap(
        #     ax,
        #     crs="EPSG:4326",
        #     source=cx.providers.Esri.WorldImagery,
        #     zorder=0,
        #     alpha=0.9,
        #     attribution=False
        # )
        print("Could not find basemap.")

    handles, labels = ax.get_legend_handles_labels()
    if handles:
        ax.legend(handles, labels)
    ax.set_xlabel("Longitude (°)")
    ax.set_ylabel("Latitude (°)")
    ax.set_title(plot_title)

    return outlier_date_summary


def plot_data(file_paths, plot_title, fig, ax, **kwargs):
    """Updates the Tkinter UI plot. Returns the outlier-date summary dict from draw_plot_elements."""
    ax.clear()
    return draw_plot_elements(ax, file_paths, plot_title, **kwargs)


def save_plot(file_paths, plot_title, output_path=None, **kwargs):
    """Saves a high-resolution plot bypassing Tkinter display quirks."""
    set_default_style()
    fig, axes = plt.subplots(figsize=(10, 8))
    draw_plot_elements(axes, file_paths, plot_title, **kwargs)
    fig.tight_layout()
    filename = Path(output_path) if output_path else make_safe_filename(plot_title)
    filename.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(filename, transparent=False, dpi=300)
    plt.close(fig)


def plot_outlier_analysis(ax, outliers, summary_stats):
    """Plots all outlier wind profiles with a highlighted average overlay."""
    if not outliers:
        ax.text(0.5, 0.5, "No valid outlier profiles found.", ha='center', va='center')
        return

    # Overlay individual outliers
    for outlier in outliers:
        speeds = np.array([outlier["WindProfile"][str(a)]["speed"] for a in ALTITUDES])
        dirs = np.array([outlier["WindProfile"][str(a)]["direction"] + 180 for a in ALTITUDES])
        dirs = np.where(dirs > 360, dirs - 360, dirs)
        theta = np.deg2rad(dirs)

        ax.plot(
            speeds,
            ALTITUDES,
            linestyle='-',
            color='grey',
            linewidth=1,
            alpha=0.2
        )
        ax.quiver(
            speeds,
            ALTITUDES,
            0.4 * np.sin(theta),
            0.4 * np.cos(theta),
            angles='uv',
            scale_units='width',
            scale=40,
            width=0.002,
            alpha=0.15,
            color='grey'
        )

    # Overlay the Aggregate Mean Line for Outliers
    mean_speeds = np.array(summary_stats["mean_speeds"])
    mean_dirs = np.array(summary_stats["mean_dirs"]) + 180
    mean_dirs = np.where(mean_dirs > 360, mean_dirs - 360, mean_dirs)
    mean_theta = np.deg2rad(mean_dirs)

    ax.plot(
        mean_speeds,
        ALTITUDES,
        linestyle='-',
        color='xkcd:cherry red',
        linewidth=3,
        label="Mean Outlier Trend"
    )
    ax.quiver(
        mean_speeds,
        ALTITUDES, 0.4 * np.sin(mean_theta), 0.4 * np.cos(mean_theta),
        angles='uv',
        scale_units='width',
        width=0.005,
        color='xkcd:dark purple',
        alpha=1
    )

    # Overlay the Aggregate Mean Line for the Entire Population
    if "pop_mean_speeds" in summary_stats and "pop_mean_dirs" in summary_stats:
        pop_speeds = np.array(summary_stats["pop_mean_speeds"])
        pop_dirs = np.array(summary_stats["pop_mean_dirs"]) + 180
        pop_dirs = np.where(pop_dirs > 360, pop_dirs - 360, pop_dirs)
        pop_theta = np.deg2rad(pop_dirs)

        ax.plot(
            pop_speeds,
            ALTITUDES,
            linestyle='--',
            color='xkcd:azure',
            linewidth=3,
            label="Population Mean"
        )
        ax.quiver(
            pop_speeds,
            ALTITUDES, 0.4 * np.sin(pop_theta), 0.4 * np.cos(pop_theta),
            angles='uv',
            scale_units='width',
            width=0.005,
            color='xkcd:sea blue',
            alpha=1
        )

    ax.legend(loc="upper right")
    ax.set_xlabel("Speed (kn)")
    ax.set_ylabel("Altitude (m)")
    ax.set_title("Overlaid Outlier Wind Profiles (>10 NM)")
    ax.grid(True)
