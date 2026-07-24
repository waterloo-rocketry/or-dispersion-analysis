import os
import numpy as np
import pandas as pd
from scipy.stats import chi2

ALTITUDES = [
    110, 320, 500, 800, 1000, 1500, 1900, 3200, 4200, 5600, 7200, 9200, 10400, 11800, 13500, 15800, 17700, 19300, 22000
]


def _read_csv(file_path):
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"The file {file_path} was not found.")

    data = pd.read_csv(
        file_path,
        encoding="utf-8",
        sep=None,
        engine="python",
        index_col=False
    )
    return data


def _find_col(data, substring, file_label=""):
    """
    Returns the first column in a DataFrame whose name contains the given substring.
    Raises a clear, actionable error instead of letting a bare `next()` throw a StopIteration
    when the expected column is missing from a CSV.
    :param data:        pandas DataFrame
    :param substring:   text to search for within column names
    :param file_label:  optional filename/context to include in the error message
    :return:            first matching Series
    """
    col_name = next((col for col in data.columns if substring in col), None)
    if col_name is None:
        context = f" in '{file_label}'" if file_label else ""
        raise ValueError(
            f"Could not find a column containing '{substring}'{context}. "
            f"Available columns: {list(data.columns)}"
        )
    return data[col_name]


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


def ellipse_math(x, y):
    """
    Function to compute properties of the ellipse required to plot it.
    :param x:   x values of population
    :param y:   y values of population
    :return:    center (x,y) coordinates of ellipse; semi-major and semi-minor axis length; rotation of ellipse
    """
    mean_x, mean_y = np.mean(x), np.mean(y)
    cov = np.cov(x, y)
    vals, vecs = np.linalg.eigh(cov)
    order = vals.argsort()[::-1]
    vals = vals[order]
    vecs = vecs[:, order]
    theta = np.degrees(np.arctan2(*vecs[:, 0][::-1]))
    return mean_x, mean_y, vals, theta


def haversine_nm(ref_latitude, ref_longitude, latitude, longitude):
    """
    Compute the great-circle distance between two points on Earth (in nautical miles).
    """
    R_nm = 3443.92
    ref_latitude, ref_longitude, latitude, longitude = map(np.radians,
                                                           [ref_latitude, ref_longitude, latitude, longitude])
    dlat = latitude - ref_latitude
    dlon = longitude - ref_longitude
    a = np.sin(dlat / 2) ** 2 + np.cos(ref_latitude) * np.cos(latitude) * np.sin(dlon / 2) ** 2
    c = 2 * np.arctan2(np.sqrt(a), np.sqrt(1 - a))
    return R_nm * c


def get_outliers(radius_nm, data, ref_latitude, ref_longitude):
    """
    Returns a list of 'Simulation' names for points outside the given radius.
    Accepts a Pandas DataFrame directly to avoid redundant file I/O.
    """
    lat_col = next((col for col in data.columns if "Landing Latitude" in col), None)
    lon_col = next((col for col in data.columns if "Landing Longitude" in col), None)

    if not lat_col or not lon_col:
        raise ValueError("Latitude or Longitude columns could not be found in the DataFrame.")

    distances = haversine_nm(ref_latitude, ref_longitude, data[lat_col], data[lon_col])

    outliers_mask = distances >= radius_nm

    if 'Simulation' in data.columns:
        outlier_sim_names = data.loc[outliers_mask, 'Simulation'].astype(str).tolist()
    else:
        outlier_sim_names = data[outliers_mask].index.astype(str).tolist()

    return outlier_sim_names


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
    # NOTE: adjust these to match your actual payload column headers -
    # this is a best guess (prefixed with "Payload ") since no payload
    # columns were previously being searched for at all (they were hardcoded to None).
    payload_lat_suffixes = tuple(f"Payload {s}" for s in lat_suffixes)
    payload_lon_suffixes = tuple(f"Payload {s}" for s in lon_suffixes)

    lat_col = next((col for col in data.columns if col.endswith(lat_suffixes)), None)
    lon_col = next((col for col in data.columns if col.endswith(lon_suffixes)), None)
    payload_lat_col = next((col for col in data.columns if col.endswith(payload_lat_suffixes)), None)
    payload_lon_col = next((col for col in data.columns if col.endswith(payload_lon_suffixes)), None)

    rocket_latitude = data[lat_col] if lat_col else None
    rocket_longitude = data[lon_col] if lon_col else None
    payload_latitude = data[payload_lat_col] if payload_lat_col else None
    payload_longitude = data[payload_lon_col] if payload_lon_col else None
    return rocket_latitude, rocket_longitude, payload_latitude, payload_longitude


def all_names(list_of_files):
    names = set()
    for file in list_of_files:
        filename = os.path.basename(file)
        components = filename.strip().split("_")
        names.add(f"{components[0]} , {components[1]}")
    return names


class RocketStats:
    def __init__(self, total_sims, mean_apogee, std_apogee, mean_landing_distance, std_landing_distance,
                 max_landing_distance, avg_lat, avg_lon, accuracy_launches, mean_min_stability, mean_lateral_velocity, mean_wind_speed):
        self.total_simulations = total_sims
        self.mean_apogee = mean_apogee
        self.std_apogee = std_apogee
        self.mean_landing_distance = mean_landing_distance
        self.std_landing_distance = std_landing_distance
        self.max_landing_distance = max_landing_distance
        self.avg_lat = avg_lat
        self.avg_lon = avg_lon
        self.accuracy_launches = accuracy_launches
        self.mean_min_stability = mean_min_stability
        self.mean_lateral_velocity = mean_lateral_velocity
        self.mean_wind_speed = mean_wind_speed


def coordinate_stats(historical_file_path):
    file_label = os.path.basename(historical_file_path)
    launch_data_frame = _read_csv(historical_file_path)

    total_sims = len(launch_data_frame)

    # Apogee
    apogee_series = _find_col(launch_data_frame, 'Apogee', file_label)
    mean_apogee = round(apogee_series.mean(), 3)
    std_apogee = round(apogee_series.std(), 3)

    # Lat and Lon for Haversine distances
    lat_series = _find_col(launch_data_frame, 'Landing Latitude', file_label)
    lon_series = _find_col(launch_data_frame, 'Landing Longitude', file_label)

    # Hardcoded Advanced Pad coordinates
    launch_lat, launch_lon = 47.965378, -81.873536
    landing_distances = haversine_nm(launch_lat, launch_lon, lat_series, lon_series)

    # NOTE: haversine_nm() returns NAUTICAL miles (R_nm = 3443.92 NM). These are
    # NOT converted to statute miles - the UI label should say "NM", not "miles"
    # (fixed in app.py's stats display). If you actually want statute miles here,
    # multiply by 1.15078 before returning.
    mean_landing_distance = round(landing_distances.mean(), 3)
    std_landing_distance = round(landing_distances.std(), 3)
    max_landing_distance = round(landing_distances.max(), 3)

    avg_lat = round(lat_series.mean(), 6)
    avg_lon = round(lon_series.mean(), 6)

    successes = (landing_distances <= 10).sum()
    accuracy_launches = successes / total_sims if total_sims > 0 else 0

    # Stability
    stability_series = _find_col(launch_data_frame, 'Min Stability', file_label)
    mean_min_stability = round(stability_series.mean(), 3)

    # Lateral Velocity
    lateral_vel_series = _find_col(launch_data_frame, 'Lateral Velocity at Apogee', file_label)
    mean_lateral_velocity = round(lateral_vel_series.mean(), 3)

    # Wind Speed
    wind_series = _find_col(launch_data_frame, 'Max Windspeed', file_label)
    mean_wind_speed = round(wind_series.mean(), 3)

    return RocketStats(
        total_sims,
        mean_apogee,
        std_apogee,
        mean_landing_distance,
        std_landing_distance,
        max_landing_distance,
        avg_lat,
        avg_lon,
        accuracy_launches,
        mean_min_stability,
        mean_lateral_velocity,
        mean_wind_speed
    )


def compute_circular_mean(degrees):
    """Computes the true average of angular data using vector math."""
    rads = np.deg2rad(degrees)
    mean_sin = np.mean(np.sin(rads))
    mean_cos = np.mean(np.cos(rads))
    mean_deg = np.rad2deg(np.arctan2(mean_sin, mean_cos))
    return (mean_deg + 360) % 360


def circular_diff(a, b):
    """Computes the shortest angular difference between two degrees."""
    return 180 - abs(abs(a - b) - 180)


def analyze_outlier_winds(historical_file, sim_parameters_file):
    """Extracts all >10 NM outliers and computes aggregate atmospheric statistics."""
    launch_lat, launch_lon = 47.965378, -81.873536

    sim_params = _read_csv(sim_parameters_file)
    sim_results = _read_csv(historical_file)

    if 'date' in sim_params.columns:
        sim_params = sim_params.set_index('date')
    if 'Simulation' in sim_results.columns:
        sim_results = sim_results.set_index('Simulation')

    outlier_names = get_outliers(
        radius_nm=10,
        data=sim_results,
        ref_latitude=launch_lat,
        ref_longitude=launch_lon
    )

    if not outlier_names:
        return [], {}  # Early exit if no outliers exist

    outlier_params = sim_params[sim_params.index.isin(outlier_names)].groupby(level=0).first()
    outlier_results = sim_results[sim_results.index.isin(outlier_names)].groupby(level=0).first()
    valid_names = outlier_params.index.intersection(outlier_results.index)

    # Compute Nominal Population Baselines
    nominal_params = sim_params[~sim_params.index.isin(outlier_names)]

    # Track per-altitude population averages for the graph
    pop_profile_speeds = []
    pop_profile_dirs = []

    for alt in ALTITUDES:
        str_a = str(alt)
        # For per-altitude profile plotting
        if str_a in nominal_params.columns:
            alt_speeds = nominal_params[str_a].dropna().values
            pop_profile_speeds.append(np.mean(alt_speeds) if len(alt_speeds) > 0 else 0)
        else:
            pop_profile_speeds.append(0)

        dir_col = f"direction [{alt}]"
        if dir_col in nominal_params.columns:
            alt_dirs = nominal_params[dir_col].dropna().values
            pop_profile_dirs.append(compute_circular_mean(alt_dirs) if len(alt_dirs) > 0 else 0)
        else:
            pop_profile_dirs.append(0)

    outliers = []
    speed_matrix = {alt: [] for alt in ALTITUDES}
    dir_matrix = {alt: [] for alt in ALTITUDES}

    for name in valid_names:
        wind_row = outlier_params.loc[name]

        wind_profile = {}
        for alt in ALTITUDES:
            str_alt = str(alt)
            speed_val = wind_row.get(str_alt, wind_row.get(alt, 0))
            dir_val = wind_row.get(f"direction [{alt}]", 0)

            wind_profile[str_alt] = {"speed": speed_val, "direction": dir_val}
            speed_matrix[alt].append(speed_val)
            dir_matrix[alt].append(dir_val)

        outliers.append({
            "WindProfile": wind_profile,
            "avg_speed": np.mean([wind_profile[str(alt)]["speed"] for alt in ALTITUDES]),
            "avg_dir": compute_circular_mean([wind_profile[str(alt)]["direction"] for alt in ALTITUDES])
        })

    summary_stats = {
        "total_outliers": len(outliers),
        "mean_speeds": [],
        "mean_dirs": [],
        "overall_max_speed": 0,
        "overall_max_speed_alt": 0,
        "pop_mean_speeds": pop_profile_speeds,  # Passed to plotting.py
        "pop_mean_dirs": pop_profile_dirs  # Passed to plotting.py
    }

    if outliers:
        for alt in ALTITUDES:
            mean_spd = np.mean(speed_matrix[alt])
            mean_dir = compute_circular_mean(dir_matrix[alt])

            summary_stats["mean_speeds"].append(mean_spd)
            summary_stats["mean_dirs"].append(mean_dir)

            if mean_spd > summary_stats["overall_max_speed"]:
                summary_stats["overall_max_speed"] = mean_spd
                summary_stats["overall_max_speed_alt"] = alt

    return outliers, summary_stats
