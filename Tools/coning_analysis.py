"""
Natural Frequency and Coning Analysis
Waterloo Rocketry
By Luca Scavone

This script was initially written by Luca. Gemini has modified it thoroughly since.

This script analyzes the pitch/yaw natural frequency (omega_n) of sounding rockets and
compares it against expected roll rates to evaluate the risk of roll-pitch coupling (coning).

It utilizes a Monte Carlo simulation approach to apply dispersions to mass properties,
aerodynamics, and atmospheric conditions, generating a 95% confidence interval for the
vehicle's natural frequency over the flight profile.

Various OpenRocket export files are needed:
    1. Complete OpenRocket export file of the vehicle through flight (make sure comment character is set to ' '.
    2. Under Tools > Component Analysis > Export:
        a. Select all vehicle components, export CNa vs Mach Number
    3. (optional) Complete OpenRocket flight exports with modified fin cant angles

Usage:
    1. Place the required simulation CSVs in the 'sim files/' directory.
    2. Adjust the configuration flags below to toggle between static PNG and MP4 animation.
    3. Run the script. Outputs are saved to the current working directory.
"""
import pandas as pd
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.animation as animation

# ==========================================
# CONFIGURATION
# ==========================================
ANIMATE         = False  # Set to True to output an MP4 animation, False for a static plot
SAVE_PLOT       = False  # Set to True to save the output (PNG or MP4) to disk
RKT_NAME        = 'Polaris'
CYCLE_NBR       = '4'

# File paths
FILE_EXPECTED   = 'sim files/Polaris_Cycle4_Expected.csv'
FILE_CNA        = 'sim files/Polaris_Cycle4_CNa-per-Component.csv'
FILE_CANT_01    = 'sim files/Polaris_Cycle4_Expected_010.csv'
FILE_CANT_02    = 'sim files/Polaris_Cycle4_Expected_020.csv'

# Plot styling
mpl.rcParams.update({
    "text.usetex": False,
    "mathtext.fontset": "cm",
    "font.family": "serif",
    "font.size": 11,
    "legend.fontsize": 9
})

# ==========================================
# DATA LOADING & PREPARATION
# ==========================================
data            = pd.read_csv(FILE_EXPECTED)
cna_data        = pd.read_csv(FILE_CNA)
cant_data_01    = pd.read_csv(FILE_CANT_01)
cant_data_02    = pd.read_csv(FILE_CANT_02)

cant_time_01    = cant_data_01[' Time (s)']
roll_rate_01    = cant_data_01['Roll rate (°/s)'] / 360  # Convert to Hz

cant_time_02    = cant_data_02[' Time (s)']
roll_rate_02    = cant_data_02['Roll rate (°/s)'] / 360  # Convert to Hz

time    = data[' Time (s)']
mach    = data['Mach number (​)']

# Convert imperial/custom units to SI metric
Izz     = data['Longitudinal moment of inertia (lb·ft²)'] * 0.0421401101  # kg·m²
Ixx     = data['Rotational moment of inertia (lb·ft²)'] * 0.0421401101  # kg·m²
vel     = data['Total velocity (ft/s)'] / 3.281  # m/s
rho     = data['Air density (g/cm³)'] * 1000  # kg/m³
cg      = data['CG location (in)'] / 39.37  # m
cp      = data['CP location (in)'] / 39.37  # m

q       = 0.5 * rho * vel ** 2
cal     = cp - cg
Aref    = np.pi * 0.1016 ** 2

# Interpolate CNa data against the main simulation Mach array
ref_mach    = cna_data[' Mach number']
ref_cna     = cna_data['CNα (Polaris Cycle 4) (​)']
cna_values  = np.interp(mach, ref_mach, ref_cna)
Cna         = pd.Series(cna_values, name='Cna')

# ==========================================
# MONTE CARLO SIMULATION
# ==========================================
num_sims = 300
wn_sims = []

for sim in range(num_sims):
    # Dispersion bounds
    m_Izz   = np.random.uniform(0.95, 1.05)
    m_Ixx   = np.random.uniform(0.95, 1.05)
    m_vel   = np.random.uniform(0.95, 1.05)
    m_rho   = np.random.uniform(0.99, 1.01)
    m_cg    = np.random.uniform(0.95, 1.05)
    m_cp    = np.random.uniform(0.9, 1.1)
    m_Cna   = np.random.uniform(0.95, 1.05)

    # Apply multipliers
    Izz_mc  = Izz * m_Izz
    Ixx_mc  = Ixx * m_Ixx
    vel_mc  = vel * m_vel
    rho_mc  = rho * m_rho
    cg_mc   = cg * m_cg
    cp_mc   = cp * m_cp
    Cna_mc  = Cna * m_Cna

    # Recalculate flight dynamics
    q_mc = 0.5 * rho_mc * vel_mc ** 2
    stability_mc = cp_mc - cg_mc

    C1_mc = q_mc * Aref * Cna_mc * stability_mc
    C1_mc = np.maximum(C1_mc, 0)  # Avoid NaNs if stability goes briefly negative

    wn_mc = np.sqrt(C1_mc / (Izz_mc + Ixx_mc)) / (2 * np.pi)
    wn_sims.append(wn_mc)

wn_sims     = np.array(wn_sims)
wn_lower    = np.percentile(wn_sims, 5, axis=0)
wn_median   = np.percentile(wn_sims, 50, axis=0)
wn_upper    = np.percentile(wn_sims, 95, axis=0)

# ==========================================
# PLOTTING
# ==========================================
fig, ax = plt.subplots(figsize=(12, 6))

# Static plot elements
ax.axvline(3.77, linestyle='-.', label='Mach 1 (3.77s)', color='blue', lw=0.9)
ax.axvline(11.28, linestyle='-.', label='MaxQ/Liquid-Phase Burnout (11.28s)', color='red', lw=0.9)
ax.axvline(14.3, linestyle='--', label='Gas-Phase Burnout (14.30s)', color='green', lw=0.9)

ax.set_title(f'{RKT_NAME} Cycle {CYCLE_NBR} Natural Frequency $\omega_n$ vs Time')
ax.set_xlabel('Time (s)')
ax.set_ylabel('Natural Frequency / Roll Rate (Hz)')
ax.set_xlim(0, 30)

if not ANIMATE:
    # Render static plot
    ax.fill_between(
        time,
        wn_lower,
        wn_upper,
        color='blue',
        alpha=0.15,
        label=r'95% Confidence Interval ($\pm 5\%$ inputs)'
    )
    ax.plot(time, wn_median, lw=1.1, color='black', label=r'Median $\omega_n$')
    ax.plot(cant_time_01, roll_rate_01, lw=0.7, color='orange', label='Roll Rate, fin cant = 0.1deg')
    ax.plot(cant_time_02, roll_rate_02, lw=0.7, color='red', label='Roll Rate, fin cant = 0.2deg')

    ax.legend(loc='upper right')
    fig.tight_layout()

    if SAVE_PLOT:
        plt.savefig(f'{RKT_NAME}_Cycle{CYCLE_NBR}_ConingAnalysis.png', dpi=300)
    plt.show()

else:
    # Render animated plot
    fill_original = ax.fill_between(
        time,
        wn_lower,
        wn_upper,
        color='blue',
        alpha=0.15,
        label=r'95% Confidence Interval ($\pm 5\%$ inputs)'
    )
    line_wn, = ax.plot(time, wn_median, lw=1.1, color='black', label=r'Median $\omega_n$')
    line_cant1, = ax.plot(cant_time_01, roll_rate_01, lw=0.7, color='orange', label='Roll Rate, fin cant = 0.1deg')
    line_cant2, = ax.plot(cant_time_02, roll_rate_02, lw=0.7, color='red', label='Roll Rate, fin cant = 0.2deg')

    ax.legend(loc='upper right')
    fig.tight_layout()

    # Tracker Dots
    dot_wn, = ax.plot([], [], 'ko', markersize=4, zorder=5)
    dot_cant1, = ax.plot([], [], 'o', color='orange', markersize=4, zorder=5)
    dot_cant2, = ax.plot([], [], 'ro', markersize=4, zorder=5)

    # Initialize empty data for frame 0
    full_time, full_wn = time.to_numpy(), wn_median
    full_cant1_time, full_cant1_roll = cant_time_01.to_numpy(), roll_rate_01.to_numpy()
    full_cant2_time, full_cant2_roll = cant_time_02.to_numpy(), roll_rate_02.to_numpy()

    line_wn.set_data([], [])
    line_cant1.set_data([], [])
    line_cant2.set_data([], [])
    fill_original.remove()
    fill_state = {'coll': None}

    fps = 30
    total_frames = fps * 5  # 5 seconds duration


    def update(frame):
        fraction = (frame + 1) / total_frames

        idx_time = max(1, int(len(full_time) * fraction))
        line_wn.set_data(full_time[:idx_time], full_wn[:idx_time])
        dot_wn.set_data([full_time[idx_time - 1]], [full_wn[idx_time - 1]])

        idx_cant1 = max(1, int(len(full_cant1_time) * fraction))
        line_cant1.set_data(full_cant1_time[:idx_cant1], full_cant1_roll[:idx_cant1])
        dot_cant1.set_data([full_cant1_time[idx_cant1 - 1]], [full_cant1_roll[idx_cant1 - 1]])

        idx_cant2 = max(1, int(len(full_cant2_time) * fraction))
        line_cant2.set_data(full_cant2_time[:idx_cant2], full_cant2_roll[:idx_cant2])
        dot_cant2.set_data([full_cant2_time[idx_cant2 - 1]], [full_cant2_roll[idx_cant2 - 1]])

        if fill_state['coll'] is not None:
            fill_state['coll'].remove()
        fill_state['coll'] = ax.fill_between(
            full_time[:idx_time],
            wn_lower[:idx_time],
            wn_upper[:idx_time],
            color='blue',
            alpha=0.15
        )
        return [line_wn, line_cant1, line_cant2, dot_wn, dot_cant1, dot_cant2, fill_state['coll']]


    ani = animation.FuncAnimation(fig, update, frames=total_frames, interval=1000 // fps, blit=False)

    if SAVE_PLOT:
        writer = animation.FFMpegWriter(
            fps=fps, codec='libx264', bitrate=5000, extra_args=['-crf', '18', '-pix_fmt', 'yuv420p']
        )
        ani.save(f'{RKT_NAME}_Cycle{CYCLE_NBR}_ConingAnalysis_Animated.mp4', writer=writer, dpi=200)

    plt.close(fig)
