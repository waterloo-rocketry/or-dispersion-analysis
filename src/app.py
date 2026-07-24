import os
import tkinter as tk
import matplotlib.pyplot as plt

from matplotlib.figure import Figure
from matplotlib.colors import to_hex
from tkinter import ttk, filedialog, messagebox
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk

from data_engine import all_names, coordinate_stats, analyze_outlier_winds
from plotting import (get_colors, make_safe_filename, plot_data, save_plot, set_default_style,
                      plot_outlier_analysis, TOP_OUTLIERS_COUNT)

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
        self.file_vars = []
        self.outlier_date_summary = {}
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

        # self.sim_param_files_selected = None

        # Allow main window to expand
        self.root.columnconfigure(0, weight=1)
        self.root.rowconfigure(0, weight=1)
        main_frame.columnconfigure(0, weight=1)
        main_frame.rowconfigure(0, weight=1)

        # Inset main frame
        middle_frame = ttk.Frame(main_frame)
        middle_frame.grid(row=0, column=0, sticky="nsew", pady=(0, 0))

        # Let left (buttons/files) and right (plot) share space dynamically
        middle_frame.columnconfigure(0, weight=0, minsize=120)
        middle_frame.columnconfigure(1, weight=1)
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
        # self.file_canvas.config(yscrollcommand=scrollbar.set)

        h_scrollbar = ttk.Scrollbar(list_frame, orient="horizontal", command=self.file_canvas.xview)
        h_scrollbar.grid(row=1, column=0, sticky="ew")

        self.file_canvas.config(yscrollcommand=scrollbar.set, xscrollcommand=h_scrollbar.set)

        # This inner frame is what actually holds the checkboxes
        self.file_check_frame = ttk.Frame(self.file_canvas)

        # Save the window ID to a variable
        self.canvas_window = self.file_canvas.create_window((0, 0), window=self.file_check_frame, anchor="nw")

        # This line makes the scroll region update whenever checkboxes are added
        self.file_check_frame.bind("<Configure>", lambda e: self.file_canvas.configure(
            scrollregion=self.file_canvas.bbox("all")))

        # Force the inner frame to be at least as wide as the canvas
        self.file_canvas.bind("<Configure>", lambda e: self.file_canvas.itemconfig(self.canvas_window, width=e.width))

        # Set plot title
        set_title = ttk.Frame(middle_frame)
        set_title.grid(row=2, column=0, sticky="nsew", padx=(0, 8), pady=(6, 6))
        set_title.columnconfigure(0, weight=1)
        set_title.rowconfigure((0, 1), weight=0)

        self.plot_title = tk.StringVar(value="Blah blah, plot title, etc. etc.")
        self.title_label = ttk.Label(set_title, text="Set Plot Title:")
        self.title_entry = ttk.Entry(set_title, textvariable=self.plot_title, width=40)

        self.title_label.grid(row=0, column=0, sticky="w", padx=8, pady=(6, 0))
        self.title_entry.grid(row=1, column=0, sticky="w", padx=8, pady=(2, 0))

        # Configurable inputs for plots
        checkbox_stack = ttk.Frame(middle_frame)
        checkbox_stack.grid(row=3, column=0, sticky="nsew", padx=(0, 8), pady=(0, 6))
        checkbox_stack.columnconfigure(0, weight=1)
        checkbox_stack.rowconfigure((0, 1, 2, 3), weight=0)

        self.plot_LC_ellipse = tk.BooleanVar(value=True)
        self.plot_sigma_ellipses = tk.BooleanVar(value=False)
        self.plot_confidence_ellipse = tk.BooleanVar(value=False)
        self.plot_top_outliers = tk.BooleanVar(value=True)

        self.LC_ellipse_box = ttk.Checkbutton(checkbox_stack, text="Plot LC Ellipse", variable=self.plot_LC_ellipse)
        self.sigma_ellipse_box = ttk.Checkbutton(checkbox_stack, text="Plot Sigma Ellipses",
                                                 variable=self.plot_sigma_ellipses)
        self.confidence_ellipse_box = ttk.Checkbutton(checkbox_stack, text="Plot Confidence Ellipse",
                                                      variable=self.plot_confidence_ellipse)
        self.top_outliers_box = ttk.Checkbutton(
            checkbox_stack, text=f"Highlight Top {TOP_OUTLIERS_COUNT} Outliers",
            variable=self.plot_top_outliers)

        self.LC_ellipse_box.grid(row=0, column=0, sticky="w", padx=4, pady=(12, 2))
        self.sigma_ellipse_box.grid(row=1, column=0, sticky="w", padx=4, pady=2)
        self.confidence_ellipse_box.grid(row=2, column=0, sticky="w", padx=4, pady=(2, 0))
        self.top_outliers_box.grid(row=3, column=0, sticky="w", padx=4, pady=(2, 0))

        self.confidence_level = tk.StringVar(value="0.95")
        self.confidence_label = ttk.Label(checkbox_stack, text="Confidence level (e.g. 0.95):")
        self.confidence_entry = ttk.Entry(checkbox_stack, textvariable=self.confidence_level, width=10)

        self.confidence_label.grid(row=4, column=0, sticky="w", padx=8, pady=(6, 0))
        self.confidence_entry.grid(row=4, column=0, sticky="e", padx=8, pady=(6, 0))
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
        """
        set_default_style()

        self.fig = Figure(figsize=(10, 8), dpi=50)
        self.ax = self.fig.add_subplot(111)

        self.canvas = FigureCanvasTkAgg(self.fig, master=self.plot_container)
        self.canvas_widget = self.canvas.get_tk_widget()
        self.canvas_widget.grid(row=0, column=0, sticky="nsew")

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
        for var in (self.plot_LC_ellipse, self.plot_sigma_ellipses, self.plot_confidence_ellipse, self.plot_top_outliers, self.plot_title):
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
        files = filedialog.askopenfilenames(
            parent=self.root, title="Select files to plot", initialdir=self.initial_dir,
            filetypes=[("CSV files", "*.csv"), ("All files", "*.*")]
        )
        if not files:
            return
        self.file_paths = list(files)
        self.launch_names = list(all_names(self.file_paths))
        self._refresh_file_listbox()


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
            var = tk.BooleanVar(value=True)

            file_name = os.path.basename(file_path)
            display_name = file_name if len(file_name) <= 38 else file_name[:35] + "..."

            cb = ttk.Checkbutton(self.file_check_frame, text=display_name, variable=var)
            cb.pack(anchor="w", fill="x", expand=True, padx=4, pady=2)

            self.file_vars.append(var)


    def _get_validated_confidence(self, confidence_flag):
        """
        Validates and returns the user-entered confidence level as a float in (0, 1).
        Returns the default (0.95) if the checkbox is unset or the field is blank.
        Shows an error dialog and returns None if the entered value is invalid,
        so callers can bail out early rather than crashing on a bad float() cast.
        :param confidence_flag: bool, whether the confidence ellipse checkbox is checked
        :return:                validated float, or None if invalid (error dialog already shown)
        """
        if not confidence_flag:
            return 0.95

        value = (self.confidence_level.get() or "").strip()
        if not value:
            return 0.95

        try:
            confidence_level = float(value)
            if not (0.0 < confidence_level < 1.0):
                raise ValueError("Confidence level must be between 0 and 1 (exclusive).")
            return confidence_level
        except Exception as e:
            messagebox.showerror(
                title="Invalid confidence level",
                message=f"Please enter a valid confidence value between 0 and 1 (e.g. 0.95).\n\nError: {e}")
            self.status_var.set("Invalid confidence value")
            return None


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

        LC_flag = self.plot_LC_ellipse.get()
        sigma_flag = self.plot_sigma_ellipses.get()
        confidence_flag = self.plot_confidence_ellipse.get()
        confidence_level = self._get_validated_confidence(confidence_flag)
        if confidence_level is None:
            return

        self.outlier_date_summary = plot_data(
            file_paths=self.active_paths,
            plot_title=self.plot_title.get(),
            fig=self.fig,
            ax=self.ax,
            plot_LC_ellipse=LC_flag,
            plot_sigma_ellipses=sigma_flag,
            plot_confidence_ellipse=confidence_flag,
            confidence=confidence_level,
            plot_top_outliers=self.plot_top_outliers.get()
        ) or {}
        self.fig.tight_layout()
        self.canvas.draw()
        self._show_stats_window()


    def _show_stats_window(self):
        if hasattr(self, '_stats_win') and self._stats_win.winfo_exists():
            self._stats_win.destroy()

        self._stats_win = tk.Toplevel(self.root)
        self._stats_win.title("General Flight Statistics")
        self._stats_win.resizable(False, False)

        frame = ttk.Frame(self._stats_win, padding=12)
        frame.pack(fill="both", expand=True)

        raw_colours, _ = get_colors(self.active_paths)
        colours = [to_hex(c) for c in raw_colours]
        self._sim_param_files = {}

        for i, file_path in enumerate(self.active_paths):
            stats = coordinate_stats(file_path)
            file_name = os.path.basename(file_path)
            display_header = file_name if len(file_name) <= 40 else file_name[:37] + "..."
            self._sim_param_files[i] = None

            # Header
            header_frame = ttk.Frame(frame)
            header_frame.pack(fill="x", pady=(10, 4))
            tk.Label(header_frame, background=colours[i], width=2, height=1).pack(side="left", padx=(0, 6))
            ttk.Label(header_frame, text=display_header, font=("", 11, "bold italic")).pack(side="left")

            # Stats Grid Layout
            grid_frame = ttk.Frame(frame)
            grid_frame.pack(fill="x", padx=16)

            rows = [
                ("Total Simulations", f"{stats.total_simulations}"),
                ("Mean Apogee", f"{stats.mean_apogee:,} ft"),
                ("Std Dev Apogee", f"{stats.std_apogee:.1f} ft"),
                ("Mean Landing Distance", f"{stats.mean_landing_distance:.1f} NM"),
                ("Std Dev Landing Dist.", f"{stats.std_landing_distance:.1f} NM"),
                ("Max Landing Distance", f"{stats.max_landing_distance:.1f} NM"),
                ("Avg Landing Coordinates", f"({stats.avg_lat}, {stats.avg_lon})"),
                ("Accuracy (within 10 NM)", f"{stats.accuracy_launches * 100:.1f}%"),
                ("Mean Min Stability", f"{stats.mean_min_stability:.3f}"),
                ("Mean Lateral Velocity", f"{stats.mean_lateral_velocity:.2f} m/s"),
                ("Mean Wind Speed", f"{stats.mean_wind_speed:.2f} mph"),
            ]

            for row_idx, (label, value) in enumerate(rows):
                ttk.Label(grid_frame, text=label).grid(row=row_idx, column=0, sticky="w", pady=2)
                ttk.Label(grid_frame, text=value, font=("", 11, "bold")).grid(row=row_idx, column=1, sticky="e", pady=2, padx=(20, 0))

            # Buttons
            btn_frame = ttk.Frame(frame)
            btn_frame.pack(pady=(8, 8))

            sim_label = ttk.Label(btn_frame, text="No file", foreground="grey", font=("", 8))

            graph_btn = ttk.Button(
                btn_frame, text="Plot Wind-Altitude Chart",
                command=lambda idx=i, fp=file_path: self._run_outlier_graph(idx, fp)
            )

            def _upload_sim(idx=i, lbl=sim_label, gb=graph_btn):
                path = filedialog.askopenfilename(parent=self._stats_win, title="Select Sim Parameters CSV")
                if path:
                    self._sim_param_files[idx] = path
                    sim_name = os.path.basename(path)
                    lbl.config(text=sim_name if len(sim_name) <= 22 else sim_name[:19] + "...", foreground="white")
                    gb.grid(row=1, column=0, columnspan=2, pady=(8, 0))

            ttk.Button(btn_frame, text="Upload Sim Params", command=_upload_sim).grid(row=0, column=0, padx=(0, 8))
            sim_label.grid(row=0, column=1)

            if file_path != self.active_paths[-1]:
                ttk.Separator(frame, orient="horizontal").pack(fill="x", pady=6)


    def _run_outlier_graph(self, i, hist_path):
        """Builds a divided window with overlay plotting and statistical analysis."""
        sim_path = self._sim_param_files[i]
        win_attr = f'_wind_win_{i}'

        if hasattr(self, win_attr) and getattr(self, win_attr).winfo_exists():
            getattr(self, win_attr).destroy()

        win = tk.Toplevel(self.root)
        win.title(f"Outlier Wind Analysis — {os.path.basename(hist_path)}")
        win.geometry("1050x650")
        setattr(self, win_attr, win)

        # Structure window: Left (Plot 75%), Right (Sidebar 25%)
        main_frame = ttk.Frame(win, padding=12)
        main_frame.pack(fill="both", expand=True)
        main_frame.columnconfigure(0, weight=3)
        main_frame.columnconfigure(1, weight=1)
        main_frame.rowconfigure(0, weight=1)

        plot_frame = ttk.Frame(main_frame)
        plot_frame.grid(row=0, column=0, sticky="nsew", padx=(0, 12))

        stats_frame = ttk.LabelFrame(main_frame, text="Outlier Analytics", padding=14)
        stats_frame.grid(row=0, column=1, sticky="nsew")

        outliers, summary = analyze_outlier_winds(hist_path, sim_path)

        # Wind-based Outlier Stats (its own sub-frame so its grid layout
        # doesn't conflict with the pack-based layout used for the section below)
        wind_stats_frame = ttk.Frame(stats_frame)
        wind_stats_frame.pack(fill="x")

        if summary.get("total_outliers", 0) > 0:
            stats_data = [
                ("Total Outliers Detected", f"{summary['total_outliers']} flights"),
                ("Worst Shear Altitude", f"{summary['overall_max_speed_alt']:,} m"),
                ("Peak Avg Layer Speed", f"{summary['overall_max_speed']:.1f} kn"),
            ]

            for row_idx, (label, val) in enumerate(stats_data):
                ttk.Label(wind_stats_frame, text=label, foreground="grey").grid(row=row_idx * 2, column=0, sticky="w", pady=(12, 0))
                ttk.Label(wind_stats_frame, text=val, font=("", 14, "bold")).grid(row=row_idx * 2 + 1, column=0, sticky="w", pady=(0, 4))

        else:
            ttk.Label(wind_stats_frame, text="No outliers >10 NM found.", font=("", 10, "italic")).pack(pady=20)

        # Top-Outlier Landing Dates (from the main dispersion plot's
        # outlier highlighting) - shown here instead of printed to the console
        outlier_dates = self.outlier_date_summary.get(hist_path, {})
        if outlier_dates:
            ttk.Separator(stats_frame, orient="horizontal").pack(fill="x", pady=(16, 8))
            dates_frame = ttk.Frame(stats_frame)
            dates_frame.pack(fill="x")

            ttk.Label(dates_frame, text=f"Top {TOP_OUTLIERS_COUNT} Outlier Landing Dates",
                     foreground="grey").grid(row=0, column=0, columnspan=2, sticky="w", pady=(0, 6))

            for row_idx, (date, count) in enumerate(sorted(outlier_dates.items()), start=1):
                ttk.Label(dates_frame, text=date).grid(row=row_idx, column=0, sticky="w", pady=1)
                ttk.Label(dates_frame, text=f"{count}", font=("", 10, "bold")).grid(
                    row=row_idx, column=1, sticky="e", padx=(20, 0), pady=1)

        # Render Plot
        fig = Figure(figsize=(8, 6))
        ax = fig.add_subplot(111)

        plot_outlier_analysis(ax, outliers, summary)

        canvas = FigureCanvasTkAgg(fig, master=plot_frame)
        canvas.get_tk_widget().pack(fill="both", expand=True)
        canvas.draw()

        toolbar_frame = ttk.Frame(plot_frame)
        toolbar_frame.pack(fill="x", side="bottom")
        NavigationToolbar2Tk(canvas, toolbar_frame).update()

        def on_close():
            fig.clear()
            win.destroy()

        win.protocol("WM_DELETE_WINDOW", on_close)


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

            confidence_level = self._get_validated_confidence(confidence_flag)
            if confidence_level is None:
                return

            plot_title = self.plot_title.get()
            default_name = make_safe_filename(plot_title).name

            save_path = filedialog.asksaveasfilename(
                title="Save plot as...",
                initialfile=default_name,
                defaultextension=".png",
                filetypes=[("PNG image", "*.png"), ("JPEG image", "*.jpg;*.jpeg"), ("All files", "*.*")]
            )

            # If user canceled the dialog, return
            if not save_path:
                self.status_var.set("Save cancelled")
                return

            # Confirm overwrite if file exists
            if os.path.exists(save_path):
                if not messagebox.askyesno(
                    title="Confirm overwrite",
                    message=f"File already exists:\n{save_path}\n\nOverwrite?"
                ):
                    self.status_var.set("Save cancelled")
                    return

            save_plot(
                file_paths=self.file_paths,
                plot_title=plot_title,
                output_path=save_path,
                plot_LC_ellipse=LC_flag,
                plot_sigma_ellipses=sigma_flag,
                plot_confidence_ellipse=confidence_flag,
                confidence=confidence_level,
                plot_top_outliers=self.plot_top_outliers.get()
            )
            self.status_var.set(f"Plot saved: '{self.plot_title.get()}'")

        except Exception as error:
            messagebox.showerror(
                title="Plot error",
                message=f"An error occurred while plotting:\n{error}"
            )
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
