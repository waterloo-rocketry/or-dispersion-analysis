"""
Sounding Rocket Damping Ratio Analysis
Waterloo Rocketry
By Luca Scavone

This script was written hastily and with not much care. Gemini cleaned it up.
This script calculates rocket pitch/yaw damping ratio over the course of its flight profile.

It accounts for both primary mechanisms of damping during atmospheric flight:
    1. Jet Damping: Restoring forces created by the motor mass flow rate.
    2. Aerodynamic Damping: Calculated using component-wise normal force coefficients (CNa)
       and Center of Pressure (CP) locations across the vehicle's geometry.

Various OpenRocket export files are needed:
    1. Complete OpenRocket export file of the vehicle through flight (make sure comment character is set to ' '.
    2. Under Tools > Component Analysis > Export:
        a. Select all vehicle components, export CP vs Mach Number
        b. Select all vehicle components, export CNa vs Mach Number

Usage:
    1. Provide path to required OpenRocket export files.
    2. Adjust the configuration flags below to toggle between static PNG and MP4 animation.
    3. Run the script. Outputs are saved to the current working directory (or you can provide an output path).
"""

import pandas as pd
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from scipy.interpolate import interp1d

# ==========================================
# CONFIGURATION [MODIFY ACCORDINGLY]
# ==========================================
ANIMATE     = False  # Set to True to output an MP4 animation, False for a static plot
SAVE_PLOT   = False  # Set to True to save the output (PNG or MP4) to disk
RKT_NAME    = 'Polaris'
CYCLE_NBR   = '4'

# File paths [MODIFY ACCORDINGLY]
FILE_SIM    = 'sim files/Polaris_Cycle4_Expected.csv'
FILE_CP     = 'sim files/Polaris_Cycle4_CP-per-Component.csv'
FILE_CNA    = 'sim files/Polaris_Cycle4_CNa-per-Component.csv'

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
sim_data    = pd.read_csv(FILE_SIM)
CP_data     = pd.read_csv(FILE_CP)
CNa_data    = pd.read_csv(FILE_CNA)

time        = sim_data[' Time (s)']
velocity    = sim_data['Total velocity (ft/s)']
mach        = sim_data['Mach number (​)']
aoa         = sim_data['Angle of attack (°)'] * np.pi / 180  # radians
prop_mass   = sim_data['Motor mass (lb)']
mdot        = np.abs(np.gradient(prop_mass, time))  # lb/s
thrust      = sim_data['Thrust (N)']
I_zz        = sim_data['Longitudinal moment of inertia (lb·ft²)']
CP_loc      = sim_data['CP location (in)'] / 12  # ft
CG_loc      = sim_data['CG location (in)'] / 12  # ft
normal_coef = sim_data['Normal force coefficient (​)']
CNa         = normal_coef / aoa
density     = sim_data['Air density (g/cm³)'] * 62.428  # lb/ft³
ref_area    = sim_data['Reference area (cm²)'] / 929  # ft²
rkt_length  = 220  # in [MODIFY ACCORDINGLY] – total rocket length
nzl_dist    = 5    # in [MODIFY ACCORDINGLY] – distance of nozzle throat from rocket aft
L_ne        = (rkt_length - nzl_dist) / 12  # nozzle distance from tip of rocket (ft)

# Interpolate CP and CNa arrays
ref_Mach    = CP_data[' Mach number']
ref_CP      = CP_data[CP_data.columns[1:]].values
ref_CNa     = CNa_data[CNa_data.columns[1:]].values

CP_interp   = interp1d(ref_Mach, ref_CP, axis=0, kind='linear', bounds_error=False, fill_value='extrapolate')
CP_values   = pd.DataFrame(CP_interp(mach), index=mach.index, columns=CP_data.columns[1:])

CNa_interp  = interp1d(ref_Mach, ref_CNa, axis=0, kind='linear', bounds_error=False, fill_value='extrapolate')
CNa_values  = pd.DataFrame(CNa_interp(mach), index=mach.index, columns=CNa_data.columns[1:])

# Component identifiers based on OpenRocket/Simulation export headers (the first column 'CP/CNa (Polaris Cycle 4) (in)'
# is ommitted because we sum up the contributions from the invidual components, not the entire rocket)
CNa_names = [
    'CNα (Nose Cone Tip Top) (​)', 'CNα (Nosecone Tip) (​)', 'CNα (Nosecone) (​)',
    'CNα (Avionics Bay) (​)', 'CNα (Canards Section) (​)', 'CNα (Canards) (​)',
    'CNα (Feedsystem) (​)', 'CNα (Fincan) (​)', 'CNα (Fins) (​)', 'CNα (Boattail) (​)'
]
CP_names = [
    'CP (Nose Cone Tip Top) (in)', 'CP (Nosecone Tip) (in)', 'CP (Nosecone) (in)',
    'CP (Avionics Bay) (in)', 'CP (Canards Section) (in)', 'CP (Canards) (in)',
    'CP (Feedsystem) (in)', 'CP (Fincan) (in)', 'CP (Fins) (in)', 'CP (Boattail) (in)'
]

# ==========================================
# DAMPING RATIO CALCULATION
# ==========================================
component_contribution = pd.Series(np.zeros(len(CG_loc)), index=CG_loc.index)

for CNa_vals, CP_vals in zip(CNa_names, CP_names):
    moment_arm = ((CP_values[CP_vals] / 12) - CG_loc) ** 2 * CNa_values[CNa_vals]
    component_contribution += moment_arm

C1          = 0.5 * density * velocity ** 2 * ref_area * CNa * (CP_loc - CG_loc)
C2_term1    = mdot * (L_ne - CG_loc) ** 2
C2_term2    = 0.5 * density * velocity * ref_area * component_contribution
Dr          = (C2_term1 + C2_term2) / (2 * np.sqrt(C1 * I_zz))

# Extract key flight events for plotting
mach1_mask = time[mach >= 1.0]
mach1_time = mach1_mask.iloc[0] if not mach1_mask.empty else None

maxq = 0.5 * density * velocity ** 2
maxq_time = time[maxq.idxmax()]

t_maxthrust = time[thrust.idxmax()]
t_burnout = time[(thrust <= 0.001) & (time > t_maxthrust)].min()

# ==========================================
# PLOTTING
# ==========================================
fig, ax = plt.subplots(figsize=(12, 6))

# Static plot elements
if mach1_time is not None:
    ax.axvline(mach1_time, linestyle='-.', label=f'Mach 1 ({mach1_time:.2f}s)', color='blue', lw=0.9)

ax.axvline(maxq_time, linestyle='-.', label=f'MaxQ / Liquid Phase Burnout ({maxq_time:.2f}s)', color='red', lw=0.9)
ax.axvline(t_burnout, linestyle='--', label=f'Vapor Phase Burnout ({t_burnout:.2f}s)', color='green', lw=0.9)

ax.set_title(f'{RKT_NAME} Cycle {CYCLE_NBR} Damping Ratio vs Time')
ax.set_xlabel('Time (s)')
ax.set_ylabel('Damping Ratio')
ax.set_xlim(0, 50)
ax.grid(True, alpha=0.3)

if not ANIMATE:
    # Render static plot
    ax.plot(time, Dr, lw=0.7, color='black', label='Damping Ratio')
    ax.legend(loc='upper right')
    fig.tight_layout()

    if SAVE_PLOT:
        plt.savefig(f'{RKT_NAME}_Cycle{CYCLE_NBR}_DampingRatio.png', dpi=300)
    plt.show()

else:
    # Render animated plot
    line_dr, = ax.plot(time, Dr, lw=0.7, color='black', label='Damping Ratio')
    ax.legend(loc='upper right')
    fig.tight_layout()

    full_time = time.to_numpy()
    full_dr = Dr.to_numpy()
    line_dr.set_data([], [])

    fps = 30
    total_frames = fps * 5  # 5 seconds duration


    def update(frame):
        fraction = (frame + 1) / total_frames
        idx = max(1, int(len(full_time) * fraction))
        line_dr.set_data(full_time[:idx], full_dr[:idx])
        return [line_dr]


    ani = animation.FuncAnimation(fig, update, frames=total_frames, interval=1000 // fps, blit=False)

    if SAVE_PLOT:
        writer = animation.FFMpegWriter(
            fps=fps, codec='libx264', bitrate=5000, extra_args=['-crf', '18', '-pix_fmt', 'yuv420p']
        )
        ani.save(f'{RKT_NAME}_Cycle{CYCLE_NBR}_DampingRatio.mp4', writer=writer, dpi=200)

    plt.close(fig)
