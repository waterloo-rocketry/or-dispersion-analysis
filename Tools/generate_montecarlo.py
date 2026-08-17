import pandas as pd
import numpy as np
from datetime import datetime
import sys
import tkinter as tk
from tkinter import filedialog


def main():
    # File Selection via Device File Explorer
    root = tk.Tk()
    root.withdraw()  # Hide the background tkinter window

    print("Opening file explorer to select input CSV...")
    input_file = filedialog.askopenfilename(
        title="Select the Wind Parameters CSV File",
        filetypes=[("CSV Files", "*.csv"), ("All Files", "*.*")]
    )

    if not input_file:
        print("\n[!] No input file selected. Exiting script.")
        sys.exit(1)

    try:
        df = pd.read_csv(input_file)
    except Exception as e:
        print(f"\n[!] Error reading the file: {e}")
        sys.exit(1)

    # Show Overview
    print("\n--- Input CSV Overview ---")
    print(df.to_string(index=False))
    print("--------------------------\n")

    # Standard Deviation Modifications
    edit_choice = input("Would you like to edit the wind standard deviations? (y/n): ").strip().lower()

    if edit_choice == 'y':
        mod_type = input(
            "Enter '1' to modify each altitude layer manually, or '2' for a blanket standard deviation: ").strip()

        if mod_type == '2':
            try:
                b_speed_std = float(input("  Enter blanket standard deviation for wind speed: "))
                b_dir_std = float(input("  Enter blanket standard deviation for wind direction: "))
                df['stddev'] = b_speed_std
                df['windDirStdDev'] = b_dir_std
                print("  [+] Blanket standard deviations applied.")
            except ValueError:
                print("  [!] Invalid input. Keeping original standard deviations.")
        elif mod_type == '1':
            for index, row in df.iterrows():
                print(f"\nAltitude: {row['altitude']} ft")

                # Speed Std Dev
                s_input = input(f"  New speed std dev (current: {row['stddev']}) [Press Enter to keep]: ").strip()
                if s_input:
                    try:
                        df.at[index, 'stddev'] = float(s_input)
                    except ValueError:
                        print("  [!] Invalid input, keeping current.")

                # Direction Std Dev
                d_input = input(
                    f"  New direction std dev (current: {row['windDirStdDev']}) [Press Enter to keep]: ").strip()
                if d_input:
                    try:
                        df.at[index, 'windDirStdDev'] = float(d_input)
                    except ValueError:
                        print("  [!] Invalid input, keeping current.")
        else:
            print("  [!] Invalid choice. Keeping original standard deviations.")

    # Temperature Setup
    print("\n--- Environmental Parameters ---")
    try:
        temp_input = input("Enter the temperature in °C (default: 15.0): ").strip()
        user_temp = float(temp_input) if temp_input else 15.0
    except ValueError:
        print("  [!] Invalid input. Defaulting to 15.0 °C.")
        user_temp = 15.0

    try:
        temp_std_input = input("Enter the temperature standard deviation (default: 0.0): ").strip()
        user_temp_std = float(temp_std_input) if temp_std_input else 0.0
    except ValueError:
        print("  [!] Invalid input. Defaulting to 0.0.")
        user_temp_std = 0.0

    # Pressure Setup
    try:
        press_input = input("\nEnter the pressure in hPa (default: 1013.2): ").strip()
        user_press = float(press_input) if press_input else 1013.2
    except ValueError:
        print("  [!] Invalid input. Defaulting to 1013.2 hPa.")
        user_press = 1013.2

    try:
        press_std_input = input("Enter the pressure standard deviation (default: 0.0): ").strip()
        user_press_std = float(press_std_input) if press_std_input else 0.0
    except ValueError:
        print("  [!] Invalid input. Defaulting to 0.0.")
        user_press_std = 0.0

    # Number of Simulations
    try:
        x_sims = int(input("\nHow many Monte Carlo simulations would you like to generate? "))
    except ValueError:
        print("  [!] Invalid number. Defaulting to 10 simulations.")
        x_sims = 10

    # Output File Location & Name
    date_str = datetime.now().strftime('%d%m%y')
    default_name = f"NavCan_MonteCarlo_{date_str}.csv"

    print("\nOpening file explorer to select save location...")
    out_path = filedialog.asksaveasfilename(
        title="Save Output Simulations As...",
        initialfile=default_name,
        defaultextension=".csv",
        filetypes=[("CSV Files", "*.csv"), ("All Files", "*.*")]
    )

    if not out_path:
        print("\n[!] No save location selected. Exiting script.")
        sys.exit(1)

    # Generate Simulations
    print(f"\nGenerating {x_sims} simulations...")

    cols = ["Simulation", "temperature", "pressure"]
    for alt_ft in df['altitude']:
        alt_m = int(round(alt_ft * 0.3048))
        cols.extend([str(alt_m), f"stdev [{alt_m}]", f"direction [{alt_m}]"])

    output_data = []

    for i in range(x_sims):
        # Generate randomized temp and pressure (Simulation 1 uses exact means)
        if i == 0:
            sim_temp = user_temp
            sim_press = user_press
        else:
            sim_temp = np.random.normal(user_temp, user_temp_std)
            sim_press = np.random.normal(user_press, user_press_std)

        row_data = [f"simulation {i + 1}", round(sim_temp, 2), round(sim_press, 2)]

        for idx, row in df.iterrows():
            mean_speed = row['speed']
            std_speed = row['stddev']
            mean_dir = row['direction']
            std_dir = row['windDirStdDev']

            if i == 0:
                # Simulation 1: No randomization, strictly use the mean values
                sim_speed = mean_speed
                sim_dir = mean_dir
            else:
                # Simulation 2 to X: Apply standard normal variation
                sim_speed = np.random.normal(mean_speed, std_speed)
                sim_speed = max(0.0, sim_speed)

                sim_dir = np.random.normal(mean_dir, std_dir)
                sim_dir = sim_dir % 360.0

            row_data.extend([round(sim_speed, 2), 0, round(sim_dir, 2)])

        output_data.append(row_data)

    # Export Data
    out_df = pd.DataFrame(output_data, columns=cols)
    out_df.to_csv(out_path, index=False)
    print(f"\n[SUCCESS] Saved {x_sims} simulations to '{out_path}'")


if __name__ == "__main__":
    main()
