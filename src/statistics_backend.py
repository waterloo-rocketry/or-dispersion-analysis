from dispersion_backend import get_outliers
import pandas as pd
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import os


def all_names(list_of_files):
    names = set()
    for file in list_of_files:
        filename = os.path.basename(file)
        components = filename.strip().split("_")
        #Split and clean into format of a file name
        names.add(f"{components[0]} , {components[1]}")

    return names


def get_outlier_indices(historical_file):
    try:
        general_df = pd.read_csv(historical_file, index_col=0)
    except Exception as e:
        print(e)
        return None

    print(general_df.columns.tolist())

    launch_lat = 47.965378
    launch_lon = -81.873536

    # I might have to change the name of these depending on the names of the columns
    latitudes = general_df[f"Polaris Rocket Landing Latitude (°N)"]
    longitudes = general_df[f"Polaris Rocket Landing Longitude (°E)"]

    outlier_indices = get_outliers(
        radius_nm=10,  # or whatever limit you want
        ref_latitude=launch_lat,
        ref_longitude=launch_lon,
        latitude=latitudes,
        longitude=longitudes
    )

    # List of all indices (on the table) where the launch was outside the 10 NM radius
    return outlier_indices


class RocketStats:
    """Storage of all the generic flight simulation statistics"""

    def __init__(self, total_sims, mean_apogee, std_apogee, mean_landing_distance, std_landing_distance,
                 max_landing_distance, accuracy_launches, mean_min_stability, mean_lateral_velocity, mean_wind_speed):
        self.total_simulations = total_sims
        self.mean_apogee = mean_apogee
        self.std_apogee = std_apogee
        self.mean_landing_distance = mean_landing_distance
        self.std_landing_distance = std_landing_distance
        self.max_landing_distance = max_landing_distance
        self.accuracy_launches = accuracy_launches
        self.mean_min_stability = mean_min_stability
        self.mean_lateral_velocity = mean_lateral_velocity
        self.mean_wind_speed = mean_wind_speed


def coordinate_stats(historical_file_path):
    try:
        launch_data_frame = pd.read_csv(historical_file_path, index_col=False)  # All landing data to sort outliers
    except Exception as e:
        print(e)
        return None

    print(launch_data_frame)

    # 1/10 | total_simulations
    total_sims = len(launch_data_frame)

    # 2/10 | mean_apogee
    mean_apogee = round(launch_data_frame['Apogee (ft)'].mean(), 2)

    # 4/10 | mean_landing_distance
    landing_distances_feet = np.sqrt(
        launch_data_frame['Polaris Rocket Position East of Launch (ft)'] ** 2 +
        launch_data_frame['Polaris Rocket Position North of Launch (ft)'] ** 2
    )
    landing_distances = landing_distances_feet * 0.000189394
    mean_landing_distance = landing_distances.mean()

    # 3/10 | std_apogee
    std_apogee = launch_data_frame['Apogee (ft)'].std()

    # 5/10 | std_landing_distance
    std_landing_distance = landing_distances.std()

    # 6/10 | max_landing_distance
    max_landing_distance = landing_distances.max()

    # 7/10 | accuracy_launches
    perimeter = 60761.2  # 10 Nautical Miles in ft.
    successes = (landing_distances <= perimeter).sum()

    try:
        accuracy_launches = successes / total_sims
    except ZeroDivisionError:
        accuracy_launches = 0

    # 8/10 | mean_min_stability
    mean_min_stability = launch_data_frame['Polaris Rocket Min Stability'].mean()

    # 9/10 | mean_lateral_velocity
    mean_lateral_velocity = launch_data_frame['Polaris Rocket Lateral Velocity at Apogee (m/s)'].mean()

    # 10/10 | mean_wind_speed
    mean_wind_speed = launch_data_frame['Max Windspeed (mph)'].mean()

    launch_stats = RocketStats(total_sims, mean_apogee, std_apogee, mean_landing_distance, std_landing_distance,
                               max_landing_distance, accuracy_launches, mean_min_stability, mean_lateral_velocity,
                               mean_wind_speed, )

    return launch_stats

class Outlying_Launch:
    def __init__(self, row, date, max_speed: dict, WindProfile: dict, AverageWindSpeed):
        self.row = row
        self.date = date
        self.max_speed = max_speed
        self.WindProfile = WindProfile
        self.AverageWindSpeed = AverageWindSpeed


def store_launch_information(outlier_indices,
                             sim_parameters_file="sim_parameters_historical_combined.csv",
                             historical_file="historical-main-at-apogee_2021-2025_35ftps.csv"):

    wind_df = pd.read_csv(sim_parameters_file, index_col=0)  # All the wind speeds at different altitudes
    general_df = pd.read_csv(historical_file)

    print("outlier_indices:", outlier_indices[:5])
    print("wind_df index:", wind_df.index.tolist()[:5])
    print("general_df index:", general_df.index.tolist()[:5])

    altitudes = [
        "110", "320", "500", "800", "1000", "1500", "1900",
        "3200", "4200", "5600", "7200", "9200", "10400",
        "11800", "13500", "15800", "17700", "19300", "22000"
    ]

    All_Outlying_launches_data = []

    for idx in outlier_indices:
        wind_profile = {}
        max_speed = {}
        max_speed["speed"] = general_df.loc[idx, "Max Windspeed (mph)"]
        max_speed["direction"] = general_df.loc[idx, "Wind Direction (deg)"]

        #Names of the columns
        speed_cols = altitudes
        #Names of the direction columns
        dir_cols = [f"direction [{a}]" for a in altitudes]
        #All speeds of all altitudes of a single launch
        row_speeds = wind_df.loc[idx, speed_cols]
        #All directions of all altitudes of a single launch
        row_dirs = wind_df.loc[idx, dir_cols]
        wind_profile = {
            alt: {"speed": row_speeds[alt], "direction": row_dirs[f"direction [{alt}]"]}
            for alt in altitudes
        }

        launch = Outlying_Launch(
            row=idx,
            date=idx,
            max_speed=max_speed,
            WindProfile=wind_profile,
            AverageWindSpeed=np.mean([wind_profile[alt]["speed"] for alt in altitudes])
        )

        All_Outlying_launches_data.append(launch)

    return All_Outlying_launches_data


def sort_outliers_winds(launch_data):
    """Returns lists ranking from least --> most in wind speeds, direction changes"""

    def get_max_speed(launch):
        """Return the max speed of a launch"""
        return launch.max_speed["speed"]

    def get_avg_speed(launch):
        """Return the average wind speed of a launch"""
        return launch.AverageWindSpeed

    sorted_by_max_speed = sorted(launch_data, key=get_max_speed)
    sorted_by_avg_speed = sorted(launch_data, key=get_avg_speed)

    return sorted_by_max_speed, sorted_by_avg_speed

def generate_plots(ax, list_indices, sorted_speeds, altitudes, name):
    for idx in list_indices:
        speeds = np.array([sorted_speeds[idx].WindProfile[str(a)]["speed"] for a in altitudes])
        directions = np.array([sorted_speeds[idx].WindProfile[str(a)]["direction"] + 180 for a in altitudes])
        directions = np.where(directions > 360, directions - 360, directions)
        name = sorted_speeds[idx].date

        # Convert to radians
        theta = np.deg2rad(directions)

        # Unit vectors (direction only)
        U = 0.4 * np.sin(theta)
        V = 0.4 * np.cos(theta)

        # Plot line connecting points
        ax.plot(speeds, altitudes, '-', linewidth=3, alpha=0.6, label=name)

        # Overlay arrows
        ax.quiver(
            speeds, altitudes, U, V,
            angles='uv',
            scale_units='width',
            scale=30,
            width=0.004,
            headwidth=3,
            headlength=4,
            alpha=1
        )

    ax.legend(loc="upper right")
    ax.set_xlabel("Speed in mph")
    ax.set_ylabel("Altitude in ft")
    ax.set_title(f"Top/Bottom {name} Speeds")
    ax.grid(True)

def plotting_top_outliers(sorted_by_max_speed, sorted_by_avg_speed, sets=3, save_path=None, ax1=None, ax2=None):
    # Code from Luca that initializes a certain format of the plots
    mpl.rcParams.update({
        "text.usetex": False,
        "mathtext.fontset": "cm",
        "font.family": "serif",
        "font.size": 11,
        "legend.fontsize": 9
    })

    altitudes = [
        110, 320, 500, 800, 1000, 1500, 1900,
        3200, 4200, 5600, 7200, 9200, 10400,
        11800, 13500, 15800, 17700, 19300, 22000
    ]

    # Prepare the list of indices: first few and last few
    list_indices = list(range(sets)) + list(range(-1, -sets - 1, -1))

    # Create a single figure with 2 subplots side by side
    if ax1 is None or ax2 is None:
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 7), sharey=True, num=f'...')
    else:
        fig = ax1.get_figure()

    # First subplot by maximum speeds
    generate_plots(ax1, list_indices, sorted_by_max_speed, altitudes, name = "Sorted by high/low Max. Speeds")
    # Second subplot by average speeds
    generate_plots(ax2, list_indices, sorted_by_avg_speed, altitudes, name = "Sorted by high/low Avg. Speeds")

    plt.tight_layout()

    # Save before show; show may clear the plot before saving it
    if save_path:
        plt.savefig(fname=save_path, dpi=300, bbox_inches='tight')

    plt.show()

    return fig