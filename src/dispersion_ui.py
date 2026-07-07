import os
import tkinter as tk
import matplotlib.pyplot as plt
from matplotlib.colors import to_hex
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk

from dispersion_backend import plot_data
from dispersion_backend import save_plot
from dispersion_backend import set_default_style
from dispersion_backend import make_safe_filename
from dispersion_backend import get_colors

from statistics_backend import sort_outliers_winds
from statistics_backend import store_launch_information
from statistics_backend import get_outlier_indices
from statistics_backend import coordinate_stats
from statistics_backend import plotting_top_outliers
from statistics_backend import all_names

from tkinter import ttk, filedialog, messagebox


class FilePlotApp:
    def __init__(self, root, initial_dir=None):
        """
        Initialize class.
        :param root:        tkinter root
        :param initial_dir: user-specified initial directory, optional - returns root directory
        """
        self.root = root
        self.root.title("Waterloo Rocketry Dispersion Zone Analysis")
        self.initial_dir = initial_dir or "."
        self.file_paths = []
        self.launch_names = []
        self.file_vars = [] # New, helps track the selected / deselected file paths

        self._build_ui()
        self._setup_matplotlib()

    def _build_ui(self):
        """
        Creates GUI with the following configuration:
            > Left Panel:
                |--> 'Select files' button to select .csv files
                |--> 'Plot data' button to plot data from .csv files
                |--> 'Clear all' button to clear plot
                |--> 'Save plot' button to save generate plot as a .png image
                |--> File display window to view selected files
                |--> Plot title user input box
                |--> Optional checkbox plot options:
                    |--> Plot arbitrary 10nm radius around Launch Canada advanced pad
                    |--> Plot 1 sigma, 2 sigma, 3 sigma dispersion ellipses
                    |--> Plot confidence ellipse around data points
                        |--> Confidence level user input box
            > Right Panel:
                |--> Plot window
            > Bottom Panel:
                |--> Status bar to inform user of errors or updates
        :return:
        """
        # Main frame
        main_frame = ttk.Frame(self.root, padding=8)
        main_frame.grid(row=0, column=0, sticky="nsew")

        # Create status_var early so init-time callbacks can use it
        self.status_var = tk.StringVar(value="Ready")

        # Hi Adrien made this
        self.sim_param_files_selected = None

        # Allow main window to expand
        self.root.columnconfigure(0, weight=1)
        self.root.rowconfigure(0, weight=1)
        main_frame.columnconfigure(0, weight=1)
        main_frame.rowconfigure(0, weight=1)

        # Inset main frame
        middle_frame = ttk.Frame(main_frame)
        middle_frame.grid(row=0, column=0, sticky="nsew", pady=(0, 0))

        # Let left (buttons/files) and right (plot) share space dynamically
        middle_frame.columnconfigure(0, weight=1, minsize=200)
        middle_frame.columnconfigure(1, weight=4)
        middle_frame.rowconfigure(0, weight=0)
        middle_frame.rowconfigure(1, weight=1)
        middle_frame.rowconfigure(2, weight=0)

        # Left column stacked buttons
        button_stack = ttk.Frame(middle_frame)
        button_stack.grid(row=0, column=0, sticky="nsew", padx=(0, 8), pady=(0, 6))
        button_stack.columnconfigure(0, weight=1)
        button_stack.rowconfigure((0, 1, 2, 3), weight=0)

        self.select_btn = ttk.Button(button_stack, text="Select Files", command=self.select_files)
        self.select_btn.grid(row=0, column=0, padx=4, pady=(0, 4), sticky="ew")

        self.plot_btn = ttk.Button(button_stack, text="Plot", command=self.plot_selected)
        self.plot_btn.grid(row=1, column=0, padx=4, pady=(0, 4), sticky="ew")

        self.clear_btn = ttk.Button(button_stack, text="Clear", command=self.clear_all)
        self.clear_btn.grid(row=2, column=0, padx=4, pady=(0, 4), sticky="ew")

        self.save_file_btn = ttk.Button(button_stack, text="Save Plot", command=self.save_file)
        self.save_file_btn.grid(row=3, column=0, padx=4, pady=(0, 4), sticky="ew")

        # Left column file list
        list_frame = ttk.LabelFrame(middle_frame, text="Selected files", padding=6)
        list_frame.grid(row=1, column=0, sticky="nsew", padx=(0, 8))
        list_frame.rowconfigure(0, weight=1)
        list_frame.columnconfigure(0, weight=1)

        # A canvas lets the inner frame scroll
        self.file_canvas = tk.Canvas(list_frame, width=220)
        self.file_canvas.grid(row=0, column=0, sticky="nsew")

        scrollbar = ttk.Scrollbar(list_frame, orient="vertical", command=self.file_canvas.yview)
        scrollbar.grid(row=0, column=1, sticky="ns")
        self.file_canvas.config(yscrollcommand=scrollbar.set)

        # This inner frame is what actually holds the checkboxes
        self.file_check_frame = ttk.Frame(self.file_canvas)
        self.file_canvas.create_window((0, 0), window=self.file_check_frame, anchor="nw")

        # This line makes the scroll region update whenever checkboxes are added
        self.file_check_frame.bind("<Configure>", lambda e: self.file_canvas.configure(
            scrollregion=self.file_canvas.bbox("all")))

        # Set plot title
        set_title = ttk.Frame(middle_frame)
        set_title.grid(row=2, column=0, sticky="nsew", padx=(0,8), pady=(6,6))
        set_title.columnconfigure(0, weight=1)
        set_title.rowconfigure((0, 1), weight=0)

        self.plot_title = tk.StringVar(value="Blah blah, plot title, etc. etc.")
        self.title_label = ttk.Label(set_title, text="Set Plot Title:")
        self.title_entry = ttk.Entry(set_title, textvariable=self.plot_title, width=40)

        self.title_label.grid(row=0, column=0, sticky="w", padx=8, pady=(6, 0))
        self.title_entry.grid(row=1, column=0, sticky="w", padx=8, pady=(2, 0))

        # Configurable inputs for plots
        checkbox_stack = ttk.Frame(middle_frame)
        checkbox_stack.grid(row=3, column=0, sticky="nsew", padx=(0,8), pady=(0,6))
        checkbox_stack.columnconfigure(0, weight=1)
        checkbox_stack.rowconfigure((0, 1, 2, 3), weight=0)

        self.plot_LC_ellipse =           tk.BooleanVar(value=True)
        self.plot_sigma_ellipses =       tk.BooleanVar(value=False)
        self.plot_confidence_ellipse =   tk.BooleanVar(value=False)

        self.LC_ellipse_box = ttk.Checkbutton(checkbox_stack, text="Plot LC Ellipse", variable=self.plot_LC_ellipse)
        self.sigma_ellipse_box = ttk.Checkbutton(checkbox_stack, text="Plot Sigma Ellipses", variable=self.plot_sigma_ellipses)
        self.confidence_ellipse_box = ttk.Checkbutton(checkbox_stack, text="Plot Confidence Ellipse", variable=self.plot_confidence_ellipse)

        self.LC_ellipse_box.grid(row=0, column=0, sticky="w", padx=4, pady=(12,2))
        self.sigma_ellipse_box.grid(row=1, column=0, sticky="w", padx=4, pady=2)
        self.confidence_ellipse_box.grid(row=2, column=0, sticky="w", padx=4, pady=(2,0))

        self.confidence_level = tk.StringVar(value="0.95")
        self.confidence_label = ttk.Label(checkbox_stack, text="Confidence level (e.g. 0.95):")
        self.confidence_entry = ttk.Entry(checkbox_stack, textvariable=self.confidence_level, width=10)

        self.confidence_label.grid(row=3, column=0, sticky="w", padx=8, pady=(6,0))
        self.confidence_entry.grid(row=3, column=0, sticky="e", padx=8, pady=(6,0))
        self.confidence_label.grid_remove()
        self.confidence_entry.grid_remove()

        self._register_input_traces()

        # Right column plot area
        plot_frame = ttk.LabelFrame(middle_frame, text="Plot", padding=6)
        plot_frame.grid(row=0, column=1, rowspan=4, sticky="nsew")
        plot_frame.rowconfigure(0, weight=1)
        plot_frame.columnconfigure(0, weight=1)
        self.plot_container = plot_frame

        # Bottom status bar (using already-created self.status_var)
        status = ttk.Label(main_frame, textvariable=self.status_var, relief=tk.SUNKEN, anchor="w")
        status.grid(row=1, column=0, sticky="ew", pady=(6, 0))

        self.root.update_idletasks()
        self.root.minsize(700, 400)

    def _setup_matplotlib(self):
        """
        Creates matplotlib Figure + Axes and embeds into Tkinter.
        :return:
        """
        set_default_style()
        self.fig, self.ax = plt.subplots(figsize=(10, 8), dpi=50)
        self.canvas = FigureCanvasTkAgg(self.fig, master=self.plot_container)
        self.canvas_widget = self.canvas.get_tk_widget()
        self.canvas_widget.grid(row=0, column=0, sticky="nsew")

        # Navigation toolbar
        toolbar_frame = ttk.Frame(self.plot_container)
        toolbar_frame.grid(row=1, column=0, sticky="ew")
        nav = NavigationToolbar2Tk(self.canvas, toolbar_frame)
        nav.update()

    def _on_plot_option_changed(self):
        """
        Updates status bar message when one of the user plot option checkboxes is updated.
        :return:
        """
        self.status_var.set("Plot options changed — press 'Plot' to apply")

    def _toggle_confidence_entry(self):
        """
        Shows and hides confidence level user input textbox depending on if self.plot_confidence_ellipse is selected.
        :return:
        """
        if self.plot_confidence_ellipse.get():
            self.confidence_label.grid()
            self.confidence_entry.grid()
        else:
            self.confidence_label.grid_remove()
            self.confidence_entry.grid_remove()
        self._on_plot_option_changed()

    def _register_input_traces(self):
        """
        Adds a trace to provided variables to track changes from user.
        :return:
        """
        for var in (self.plot_LC_ellipse, self.plot_sigma_ellipses, self.plot_confidence_ellipse, self.plot_title):
            try:
                var.trace_add('write', lambda *a, v=var: self._on_plot_option_changed())
            except AttributeError:
                var.trace('w', lambda *a, v=var: self._on_plot_option_changed())

        try:
            self.plot_confidence_ellipse.trace_add('write', lambda *a: self._toggle_confidence_entry())
        except AttributeError:
            self.plot_confidence_ellipse.trace('w', lambda *a: self._toggle_confidence_entry())

        self._toggle_confidence_entry()

    def select_files(self):
        """
        Function enabling user-selected files using filedialog.
        :return:
        """
        # Open file dialog; allow multiple selection
        files = filedialog.askopenfilenames(
            parent=self.root,
            title="Select files to plot",
            initialdir=self.initial_dir,
            filetypes=[("CSV files", "*.csv"), ("All files", "*.*")]
        )
        if not files:
            self.status_var.set("No files selected")
            return

        # Convert to list and update internal state & UI
        self.file_paths = list(files)
        self.launch_names = list(all_names(self.file_paths))
        self._refresh_file_listbox()
        self.status_var.set(f"Selected {len(self.file_paths)} files")

    def _refresh_file_listbox(self):
        """
        Adds selected .csv filenames to the file display window so the user can see what files have been selected.
        :return:
        """
        # Wipe all existing checkboxes from the frame
        for widget in self.file_check_frame.winfo_children():
            widget.destroy()
        self.file_vars = []

        # Create one checkbox per file
        for file_path in self.file_paths:
            var = tk.BooleanVar(value=True)  # checked by default
            cb = ttk.Checkbutton(self.file_check_frame, text=os.path.basename(file_path), variable=var)
            cb.pack(anchor="w", padx=4, pady=2)
            self.file_vars.append(var)

    def plot_selected(self):
        """
        Function is called when 'Plot data' button is clicked; calls plot_data() to post-process data and display on
        the graph.
        :return:
        """
        self.active_paths = [p for p, v in zip(self.file_paths, self.file_vars) if v.get()]
        if not self.active_paths:
            messagebox.showwarning("No files checked", "Please check at least one file to plot.")
            return

        try:
            # clear axes, call user plotting routine, draw canvas
            self.ax.clear()

            LC_flag = self.plot_LC_ellipse.get()
            sigma_flag = self.plot_sigma_ellipses.get()
            confidence_flag = self.plot_confidence_ellipse.get()

            # Validate confidence level if checkbox is set
            confidence_level = None
            if confidence_flag:
                value = (self.confidence_level.get() or "").strip()
                if not value:
                    confidence_level = 0.95  # default value
                else:
                    try:
                        confidence_level = float(value)
                        if not (0.0 < confidence_level < 1.0):
                            raise ValueError("Confidence level must be between 0 and 1 (exclusive).")
                    except Exception as e:
                        messagebox.showerror("Invalid confidence level",
                                             f"Please enter a valid confidence value between 0 and 1 (e.g. 0.95).\n\nError: {e}")
                        self.status_var.set("Invalid confidence value")
                        return

            print("active_paths:", self.active_paths)
            print("file_vars:", [v.get() for v in self.file_vars])
            print("file_paths:", self.file_paths)

            plot_data(
                file_paths=self.active_paths,
                plot_title=self.plot_title.get(),
                fig=self.fig,
                ax=self.ax,
                plot_LC_ellipse=LC_flag,
                plot_sigma_ellipses=sigma_flag,
                plot_confidence_ellipse=confidence_flag,
                confidence=confidence_level
            )
            self.fig.tight_layout()
            self.canvas.draw()
            self.status_var.set("Plotted selected files")
            self._show_stats_window() # This is the new window that Adrien is working on
        except Exception as error:
            messagebox.showerror("Plot error", f"An error occurred while plotting:\n{error}")
            self.status_var.set("Plot error")

    def _select_sim_param_files(self):
        files = filedialog.askopenfilenames(
            parent=self.root,
            title="Select files to analyze",
            initialdir=self.initial_dir,
            filetypes=[("CSV files", "*.csv"), ("All files", "*.*")]
        )
        if not files:
            messagebox.showwarning("Upload Simulation Parameters",
                                   "Please upload all Simulation Parameter files to chart.")
            self.sim_param_files_selected = False
            return

        self.sim_param_paths = list(files)
        self.sim_param_files_selected = True  # ← removed _refresh_file_listbox() call

    def _create_wind_altitude_window(self):
        if hasattr(self, '_wind_altitude_window') and self._wind_altitude_window.winfo_exists():
            self._wind_altitude_window.destroy()

        if not self.sim_param_files_selected:
            messagebox.showwarning("No simulation parameter files selected")
            return None

        #Need to fix the name
        self._wind_altitude_window = tk.Toplevel(self.root)
        self._wind_altitude_window.title(f"Wind-Altitude Breakdown")
        self._wind_altitude_window.resizable(True, True)

        frame = ttk.Frame(self._wind_altitude_window, padding=12)
        frame.rowconfigure(0, weight=1)
        frame.columnconfigure(0, weight=1)

        # Calling of the outliers backend functions
        outlier_indices = get_outlier_indices(self.active_paths[0])  # or whichever file is relevant
        launch_data = store_launch_information(outlier_indices, sim_parameters_file=self.sim_param_paths[0], historical_file=self.active_paths[0])
        sorted_by_max_speed, sorted_by_avg_speed = sort_outliers_winds(launch_data)

        fig, (ax1, ax2) = plt.subplots(1, 2, sharey=True, figsize=(14, 7))
        plotting_top_outliers(sorted_by_max_speed, sorted_by_avg_speed, ax1=ax1, ax2=ax2)
        self._wind_fig = fig

        frame.pack(fill="both", expand=True)

        canvas = FigureCanvasTkAgg(self._wind_fig, master=frame)
        canvas.get_tk_widget().grid(row=0, column=0, sticky="nsew")
        canvas.draw()

        # ADD: NavigationToolbar2Tk in a toolbar_frame below the canvas, same as _setup_matplotlib

        # BIND: self._wind_altitude_window.protocol("WM_DELETE_WINDOW", lambda: [plt.close(self._wind_fig), self._wind_altitude_window.destroy()])

    def _show_stats_window(self):
        # Guard against duplicate windows
        if hasattr(self, '_stats_win') and self._stats_win.winfo_exists():
            self._stats_win.destroy()

        self._stats_win = tk.Toplevel(self.root)
        self._stats_win.title("General Flight Statistics")
        self._stats_win.resizable(False, False)

        frame = ttk.Frame(self._stats_win, padding=12)
        frame.pack(fill="both", expand=True)

        raw_colours, _ = get_colors(self.file_paths)
        colours = [to_hex(c) for c in raw_colours]
        current_row = 0

        self._sim_param_files = {}  # key: i -> sim_params path

        for i, file_path in enumerate(self.active_paths):
            stats = coordinate_stats(file_path)
            file_name = os.path.basename(file_path)
            colour = colours[i]

            self._sim_param_files[i] = None

            # Header row: colour box + filename
            header_frame = ttk.Frame(frame)
            header_frame.grid(row=current_row, column=0, columnspan=2, sticky="w", pady=(10, 4))
            tk.Label(header_frame, background=colour, width=2, height=1).pack(side="left", padx=(0, 6))
            ttk.Label(header_frame, text=file_name, font=("", 10, "bold italic")).pack(side="left")
            current_row += 1

            # Stat rows
            rows = [
                ("Total Simulations", f"{stats.total_simulations}"),
                ("Mean Apogee", f"{stats.mean_apogee:,} ft"),
                ("Std Dev Apogee", f"{stats.std_apogee:.1f} ft"),
                ("Mean Landing Distance", f"{stats.mean_landing_distance:.1f} miles"),
                ("Std Dev Landing Dist.", f"{stats.std_landing_distance:.1f} miles"),
                ("Max Landing Distance", f"{stats.max_landing_distance:.1f} miles"),
                ("Accuracy (within 10 NM)", f"{stats.accuracy_launches * 100:.1f}%"),
                ("Mean Min Stability", f"{stats.mean_min_stability:.3f}"),
                ("Mean Lateral Velocity", f"{stats.mean_lateral_velocity:.2f} m/s"),
                ("Mean Wind Speed", f"{stats.mean_wind_speed:.2f} mph"),
            ]

            for label, value in rows:
                ttk.Label(frame, text=label, anchor="w").grid(row=current_row, column=0, sticky="w", pady=2,
                                                              padx=(0, 16))
                ttk.Label(frame, text=value, anchor="e", font=("", 10, "bold")).grid(row=current_row, column=1,
                                                                                     sticky="e", pady=2)
                current_row += 1

            # ── Per-block upload + graph buttons ────────────────────
            btn_frame = ttk.Frame(frame)
            btn_frame.grid(row=current_row, column=0, columnspan=2, sticky="ew", pady=(8, 2))

            sim_label = ttk.Label(btn_frame, text="No file", foreground="grey", font=("", 8))

            graph_btn = ttk.Button(
                btn_frame,
                text="Plot Wind-Altitude Chart",
                state="disabled",
                command=lambda i=i, fp=file_path: self._run_outlier_graph(i, fp)
            )

            def _upload_sim(i=i, lbl=sim_label, gb=graph_btn):
                path = filedialog.askopenfilename(
                    parent=self._stats_win,
                    title="Select Sim Parameters CSV",
                    filetypes=[("CSV files", "*.csv"), ("All files", "*.*")]
                )
                if path:
                    self._sim_param_files[i] = path
                    lbl.config(text=os.path.basename(path), foreground="black")
                    gb.config(state="normal")

            ttk.Button(btn_frame, text="Upload Sim Params", command=_upload_sim).pack(side="left", padx=(0, 8))
            sim_label.pack(side="left", padx=(0, 16))
            graph_btn.pack(side="left")

            current_row += 1
            # ────────────────────────────────────────────────────────

            # Divider between files, skip on last file
            if file_path != self.active_paths[-1]:
                ttk.Separator(frame, orient="horizontal").grid(
                    row=current_row, column=0, columnspan=2, sticky="ew", pady=6)
                current_row += 1

    def _run_outlier_graph(self, i, hist_path):
        sim_path = self._sim_param_files[i]

        # Guard against duplicate windows per block
        win_attr = f'_wind_win_{i}'
        if hasattr(self, win_attr) and getattr(self, win_attr).winfo_exists():
            getattr(self, win_attr).destroy()

        win = tk.Toplevel(self.root)
        win.title(f"Wind-Altitude Breakdown — {os.path.basename(hist_path)}")
        win.resizable(True, True)
        setattr(self, win_attr, win)

        frame = ttk.Frame(win, padding=12)
        frame.rowconfigure(0, weight=1)
        frame.columnconfigure(0, weight=1)
        frame.pack(fill="both", expand=True)

        outlier_indices = get_outlier_indices(hist_path)
        launch_data = store_launch_information(outlier_indices, sim_path, hist_path)
        sorted_max, sorted_avg = sort_outliers_winds(launch_data)

        fig, (ax1, ax2) = plt.subplots(1, 2, sharey=True, figsize=(14, 7))
        plotting_top_outliers(sorted_max, sorted_avg, ax1=ax1, ax2=ax2)

        canvas = FigureCanvasTkAgg(fig, master=frame)
        canvas.get_tk_widget().grid(row=0, column=0, sticky="nsew")
        canvas.draw()

        toolbar_frame = ttk.Frame(frame)
        toolbar_frame.grid(row=1, column=0, sticky="ew")
        NavigationToolbar2Tk(canvas, toolbar_frame).update()

        win.protocol("WM_DELETE_WINDOW", lambda: [plt.close(fig), win.destroy()])

    def save_file(self):
        """
        Function is called when 'Save plot' button is clicked; calls save_plot() to save the figure.
        :return:
        """
        if not self.file_paths:
            messagebox.showwarning("No files", "Please select files before plotting.")
            return

        try:
            LC_flag = self.plot_LC_ellipse.get()
            sigma_flag = self.plot_sigma_ellipses.get()
            confidence_flag = self.plot_confidence_ellipse.get()

            # Validate confidence level if checkbox is set
            confidence_level = None
            if confidence_flag:
                value = (self.confidence_level.get() or "").strip()
                if not value:
                    confidence_level = 0.95  # default value
                else:
                    try:
                        confidence_level = float(value)
                        if not (0.0 < confidence_level < 1.0):
                            raise ValueError("Confidence level must be between 0 and 1 (exclusive).")
                    except Exception as e:
                        messagebox.showerror("Invalid confidence level",
                                             f"Please enter a valid confidence value between 0 and 1 (e.g. 0.95).\n\nError: {e}")
                        self.status_var.set("Invalid confidence value")
                        return

            plot_title = self.plot_title.get()
            default_name = make_safe_filename(plot_title).name

            save_path = filedialog.asksaveasfilename(
                title="Save plot as...",
                initialfile=default_name,
                defaultextension=".png",
                filetypes=[("PNG image", "*.png"), ("JPEG image", "*.jpg;*.jpeg"), ("All files", "*.*")]
            )

            # If user cancelled the dialog, return
            if not save_path:
                self.status_var.set("Save cancelled")
                return

            # Confirm overwrite if file exists
            if os.path.exists(save_path):
                if not messagebox.askyesno("Confirm overwrite", f"File already exists:\n{save_path}\n\nOverwrite?"):
                    self.status_var.set("Save cancelled")
                    return

            save_plot(
                file_paths=self.file_paths,
                plot_title=plot_title,
                output_path=save_path,
                plot_LC_ellipse=LC_flag,
                plot_sigma_ellipses=sigma_flag,
                plot_confidence_ellipse=confidence_flag,
                confidence=confidence_level)
            self.status_var.set(f"Plot saved: '{self.plot_title.get()}'")

        except Exception as error:
            messagebox.showerror("Plot error", f"An error occurred while plotting:\n{error}")
            self.status_var.set("Plot error")

    def clear_all(self):
        """
        Clears all inputs and plots, resets to default values.
        :return:
        """
        self.file_paths = []
        for widget in self.file_check_frame.winfo_children():
            widget.destroy()
        self.file_vars = []
        self.ax.clear()
        self.ax.relim()
        self.fig.tight_layout()
        self.canvas.draw()
        self.status_var.set("Cleared files and plot")


def main():
    root = tk.Tk()

    initial_dir = "/Users/luca/PycharmProjects/Rocketry Dispersion Zone Analysis/Plugin Exports/"
    FilePlotApp(root, initial_dir=initial_dir)
    root.geometry("1200x800")
    root.mainloop()


if __name__ == "__main__":
    main()
