import os
import threading

import pandas as pd

from PySide6.QtCore import Qt, QSettings, QObject, QThread, Signal
from PySide6.QtWidgets import (
    QMainWindow, QWidget, QVBoxLayout, QHBoxLayout, QGridLayout,
    QPushButton, QLabel, QLineEdit, QCheckBox, QGroupBox, QScrollArea,
    QFrame, QFileDialog, QMessageBox, QDialog, QApplication
)

from matplotlib.figure import Figure
from matplotlib.colors import to_hex
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg, NavigationToolbar2QT

from data_engine import (
    all_names,
    coordinate_stats,
    analyze_outlier_winds,
    _read_csv,
    _find_col,
    LC_WAIVER_RADIUS_NM,
    compute_water_landings
)
from plotting import (
    get_colors,
    make_safe_filename,
    plot_data,
    save_plot,
    set_default_style,
    plot_outlier_analysis,
    TOP_OUTLIERS_COUNT,
    LC_GEOGRAPHY_DIR
)


# Centralized dark theme, applied once at the QApplication level (see main.py) so every
# window - the main app, the stats dock, and the pop-up outlier windows - gets consistent
# styling automatically. Previously, individual labels had colors like "white" or "grey"
# hardcoded inline, which assumed a dark background; if the OS was in light mode those
# labels could become unreadable. Giving the app its own dark background here makes that
# assumption safe everywhere instead of accidentally depending on the OS theme.
APP_STYLESHEET = """
QMainWindow, QDialog, QWidget {
    background-color: #1e1f22;
    color: #e6e6e6;
}
QGroupBox {
    border: 1px solid #3c3f41;
    border-radius: 6px;
    margin-top: 10px;
    padding-top: 8px;
    font-weight: bold;
}
QGroupBox::title {
    subcontrol-origin: margin;
    left: 10px;
    padding: 0 4px;
    color: #cfd2d6;
}
QPushButton {
    background-color: #3c3f41;
    color: #e6e6e6;
    border: 1px solid #4b4f52;
    border-radius: 5px;
    padding: 6px 10px;
}
QPushButton:hover {
    background-color: #4b4f52;
}
QPushButton:pressed {
    background-color: #2c2e30;
}
QPushButton:checked {
    background-color: #5a8fd6;
    border: 1px solid #7aa6e0;
}
QPushButton:disabled {
    color: #7a7a7a;
}
QLineEdit {
    background-color: #2b2d30;
    color: #e6e6e6;
    border: 1px solid #4b4f52;
    border-radius: 4px;
    padding: 3px;
}
QCheckBox {
    color: #e6e6e6;
}
QScrollArea {
    background-color: #1e1f22;
    border: none;
}
QLabel#mutedLabel {
    color: #9a9a9a;
}
QLabel#simFileLabel {
    color: #9a9a9a;
    font-size: 8pt;
}
QLabel#simFileLabelActive {
    color: #f0f0f0;
    font-size: 8pt;
}
"""


class _CsvLoaderWorker(QObject):
    """
    Loads a batch of CSV files on a background thread so large files don't freeze the UI.
    Only handles the disk-I/O part - matplotlib rendering and Qt widget updates still happen
    back on the GUI thread (via the `finished` signal) since neither is thread-safe.
    """
    finished = Signal(dict, str)  # (path -> DataFrame, error message or "" on success)

    def __init__(self, paths, cache, cache_lock):
        super().__init__()
        self._paths = paths
        self._cache = cache
        self._cache_lock = cache_lock

    def run(self):
        try:
            result = {}
            for path in self._paths:
                with self._cache_lock:
                    if path not in self._cache:
                        self._cache[path] = _read_csv(path)
                    result[path] = self._cache[path]
            self.finished.emit(result, "")
        except Exception as e:
            self.finished.emit({}, str(e))


class OutlierWindow(QDialog):
    """
    Non-modal window for the per-file outlier wind analysis plot. Clears its
    matplotlib Figure on close (mirrors the original Tkinter
    WM_DELETE_WINDOW -> fig.clear() handler) so unused figures don't linger.
    """
    def __init__(self, parent, fig):
        super().__init__(parent)
        self._fig = fig
        # Non-modal dialogs are only hidden on close by default, not destroyed - without
        # this, repeatedly opening wind-analysis windows across several plot/replot cycles
        # would accumulate hidden-but-still-alive QDialog objects over the session.
        self.setAttribute(Qt.WidgetAttribute.WA_DeleteOnClose)

    def closeEvent(self, event):
        self._fig.clear()
        super().closeEvent(event)


class FilePlotApp(QMainWindow):
    def __init__(self, initial_dir=None):
        """
        Initialize class.
        :param initial_dir: user-specified initial directory, optional. If omitted, falls back
                            to the last directory the user selected files from (persisted across
                            runs via QSettings), then to the user's home directory on first launch.
        """
        super().__init__()
        self.setWindowTitle("Waterloo Rocketry Dispersion Zone Analysis")

        self.settings = QSettings("WaterlooRocketry", "DispersionZoneAnalysis")
        self.initial_dir = initial_dir or self.settings.value("last_directory", os.path.expanduser("~"))

        self.file_paths = []
        self.launch_names = []
        self.file_checks = []
        self.active_paths = []
        self.outlier_date_summary = {}
        self._sim_param_files = {}
        self._csv_cache = {}
        self._csv_cache_lock = threading.Lock()
        self._water_cache = {}
        self._wind_windows = {}
        self._loader_thread = None
        self._loader_worker = None
        self._build_ui()
        self._setup_matplotlib()
        self._setup_stats_window()
        self.resize(1200, 800)


    def _build_ui(self):
        """
        Creates GUI with the following configuration:
            > Left Panel:
                |--> 'Select files' button to select .csv files
                |--> 'Plot data' button to plot data from .csv files
                |--> 'Clear all' button to clear plot
                |--> 'Save plot' button to save generate plot as a .png image
                |--> 'Export Stats' button to export the current stats panel as a .csv file
                |--> 'Show Stats Panel' toggle to open/close the docked stats sidebar
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
        central = QWidget()
        self.setCentralWidget(central)
        grid = QGridLayout(central)
        grid.setContentsMargins(8, 8, 8, 8)

        # Left column stacked buttons
        button_stack = QWidget()
        button_layout = QVBoxLayout(button_stack)
        button_layout.setContentsMargins(0, 0, 8, 6)

        self.select_btn = QPushButton("Select Files")
        self.plot_btn = QPushButton("Plot")
        self.clear_btn = QPushButton("Clear")
        self.save_file_btn = QPushButton("Save Plot")
        self.export_stats_btn = QPushButton("Export Stats")
        self.stats_toggle_btn = QPushButton("Show Stats Window")
        self.stats_toggle_btn.setCheckable(True)

        for btn in (self.select_btn, self.plot_btn, self.clear_btn, self.save_file_btn,
                    self.export_stats_btn, self.stats_toggle_btn):
            button_layout.addWidget(btn)

        grid.addWidget(button_stack, 0, 0)

        # Left column file list (scrollable, checkbox per file)
        list_group = QGroupBox("Selected files")
        list_group.setMinimumWidth(220)
        list_group_layout = QVBoxLayout(list_group)

        self.file_scroll = QScrollArea()
        self.file_scroll.setWidgetResizable(True)
        self.file_list_widget = QWidget()
        self.file_list_layout = QVBoxLayout(self.file_list_widget)
        self.file_list_layout.addStretch()  # keeps checkboxes anchored to the top
        self.file_scroll.setWidget(self.file_list_widget)
        list_group_layout.addWidget(self.file_scroll)

        grid.addWidget(list_group, 1, 0)

        # Set plot title
        title_widget = QWidget()
        title_layout = QVBoxLayout(title_widget)
        title_layout.setContentsMargins(0, 6, 8, 6)
        self.title_entry = QLineEdit("Blah blah, plot title, etc. etc.")
        title_layout.addWidget(QLabel("Set Plot Title:"))
        title_layout.addWidget(self.title_entry)

        grid.addWidget(title_widget, 2, 0)

        # Configurable inputs for plots
        checkbox_stack = QWidget()
        checkbox_layout = QVBoxLayout(checkbox_stack)
        checkbox_layout.setContentsMargins(0, 0, 8, 6)

        self.LC_ellipse_box = QCheckBox("Plot LC Ellipse")
        self.LC_ellipse_box.setChecked(True)
        self.sigma_ellipse_box = QCheckBox("Plot Sigma Ellipses")
        self.confidence_ellipse_box = QCheckBox("Plot Confidence Ellipse")
        self.water_landings_box = QCheckBox("Plot Water Landings (cyan)")
        self.water_landings_box.setChecked(True)
        self.top_outliers_box = QCheckBox(f"Highlight Top {TOP_OUTLIERS_COUNT} Outliers (yellow)")
        self.top_outliers_box.setChecked(True)

        for box in (
                self.LC_ellipse_box,
                self.sigma_ellipse_box,
                self.confidence_ellipse_box,
                self.water_landings_box,
                self.top_outliers_box
        ):
            checkbox_layout.addWidget(box)

        self.confidence_container = QWidget()
        conf_layout = QHBoxLayout(self.confidence_container)
        conf_layout.setContentsMargins(8, 6, 0, 0)
        self.confidence_entry = QLineEdit("0.95")
        self.confidence_entry.setFixedWidth(40)
        conf_layout.addWidget(QLabel("Confidence (e.g. 0.95):"))
        conf_layout.addWidget(self.confidence_entry)
        conf_layout.addStretch()
        self.confidence_container.setVisible(False)
        checkbox_layout.addWidget(self.confidence_container)

        grid.addWidget(checkbox_stack, 3, 0)

        # Right column plot area
        plot_group = QGroupBox("Plot")
        self.plot_layout = QVBoxLayout(plot_group)
        grid.addWidget(plot_group, 0, 1, 4, 1)

        grid.setColumnStretch(0, 0)
        grid.setColumnStretch(1, 1)
        grid.setRowStretch(1, 1)

        # Bottom status bar
        self.statusBar().showMessage("Ready")

        self._register_input_connections()
        self.setMinimumSize(700, 400)


    def _setup_matplotlib(self):
        """
        Creates matplotlib Figure + Axes and embeds into the Qt plot panel.
        """
        set_default_style()

        self.fig = Figure(figsize=(10, 8), dpi=100)
        self.ax = self.fig.add_subplot(111)

        self.canvas = FigureCanvasQTAgg(self.fig)
        self.toolbar = NavigationToolbar2QT(self.canvas, self)

        self.plot_layout.addWidget(self.canvas)
        self.plot_layout.addWidget(self.toolbar)


    def _setup_stats_window(self):
        """
        Builds the toggleable 'Flight Statistics' pop-up window.
        Uses a single reusable QDialog so it floats over the main UI without
        messing with the plot's aspect ratio, while avoiding desktop clutter.
        """
        self.stats_window = QDialog(self)
        self.stats_window.setWindowTitle("Flight Statistics")
        self.stats_window.setMinimumSize(400, 600)

        # Setup layout to hold the scroll area
        layout = QVBoxLayout(self.stats_window)
        layout.setContentsMargins(0, 0, 0, 0)

        self.stats_scroll = QScrollArea()
        self.stats_scroll.setWidgetResizable(True)
        layout.addWidget(self.stats_scroll)

        # Uncheck the toggle button if the user closes the window via the 'X' button
        self.stats_window.finished.connect(lambda _: self.stats_toggle_btn.setChecked(False))

        # Tie the toggle button to the window's visibility
        def toggle_window(checked):
            if checked:
                self.stats_window.show()
            else:
                self.stats_window.hide()

        self.stats_toggle_btn.toggled.connect(toggle_window)


    def _on_plot_option_changed(self, *_args):
        """
        Updates status bar message when one of the user plot option widgets is updated.
        :return:
        """
        self.statusBar().showMessage("Plot options changed — press 'Plot' to apply")


    def _toggle_confidence_entry(self, checked):
        """
        Shows and hides confidence level user input box depending on if
        self.confidence_ellipse_box is checked.
        :return:
        """
        self.confidence_container.setVisible(checked)
        self._on_plot_option_changed()


    def _register_input_connections(self):
        """
        Wires widget signals to their change handlers.
        :return:
        """
        for box in (self.LC_ellipse_box, self.sigma_ellipse_box, self.water_landings_box, self.top_outliers_box):
            box.toggled.connect(self._on_plot_option_changed)
        self.confidence_ellipse_box.toggled.connect(self._toggle_confidence_entry)
        self.title_entry.textChanged.connect(self._on_plot_option_changed)

        self.select_btn.clicked.connect(self.select_files)
        self.plot_btn.clicked.connect(self.plot_selected)
        self.clear_btn.clicked.connect(self.clear_all)
        self.save_file_btn.clicked.connect(self.save_file)
        self.export_stats_btn.clicked.connect(self.export_stats)


    def _load_csv(self, file_path):
        """
        Returns a cached DataFrame for file_path, reading it from disk only the first time it's
        requested. The same file previously got re-read from disk separately by the main plot,
        the stats panel, and the outlier-wind analysis - each triggering its own full CSV parse.
        Routing every read through this cache means each file is only ever parsed once per
        selection.
        """
        with self._csv_cache_lock:
            if file_path not in self._csv_cache:
                self._csv_cache[file_path] = _read_csv(file_path)
            return self._csv_cache[file_path]


    def _get_water_landing(self, file_path):
        """
        Returns a cached (water_count, water_prob, water_indices) result for file_path,
        computing it only the first time it's requested for the current file selection.
        This spatial join was previously recomputed independently for the main plot's
        cyan water-landing overlay and again for the stats panel's "Water Landing
        Probability" figure - once each, every time either was refreshed.
        """
        if file_path not in self._water_cache:
            data = self._load_csv(file_path)
            file_name = os.path.basename(file_path)
            lat_series = _find_col(data, 'Landing Latitude', file_name)
            lon_series = _find_col(data, 'Landing Longitude', file_name)
            lakes_file = LC_GEOGRAPHY_DIR / "lakes.geojson"
            self._water_cache[file_path] = compute_water_landings(data, lat_series, lon_series, lakes_file, file_name)
        return self._water_cache[file_path]


    def _set_busy(self, busy):
        """Disables the actions that trigger CSV loads and shows a busy cursor/status message."""
        for btn in (self.select_btn, self.plot_btn, self.clear_btn, self.save_file_btn, self.export_stats_btn):
            btn.setEnabled(not busy)
        if busy:
            QApplication.setOverrideCursor(Qt.CursorShape.WaitCursor)
            self.statusBar().showMessage("Loading files…")
        else:
            QApplication.restoreOverrideCursor()


    def _load_csvs_async(self, paths, on_done):
        """
        Loads `paths` on a background QThread (with a busy cursor and disabled controls, since
        parsing large CSVs was previously blocking the whole UI), then calls on_done(data_by_path)
        back on the GUI thread once finished.
        """
        self._set_busy(True)

        # Store the callback as an instance variable so the new handler method can access it
        self._current_on_done = on_done

        self._loader_thread = QThread(self)
        self._loader_worker = _CsvLoaderWorker(paths, self._csv_cache, self._csv_cache_lock)
        self._loader_worker.moveToThread(self._loader_thread)
        self._loader_thread.started.connect(self._loader_worker.run)

        # Connect to a formal class method so Qt natively knows to route it to the main thread
        self._loader_worker.finished.connect(self._on_worker_finished)

        self._loader_thread.finished.connect(self._loader_thread.deleteLater)
        self._loader_thread.start()


    def _on_worker_finished(self, data_by_path, error):
        """
        Handles the completion of the background worker explicitly on the main GUI thread.
        """
        self._set_busy(False)
        self._loader_thread.quit()
        self._loader_thread.wait()

        if error:
            QMessageBox.critical(
                self, "File load error",
                f"Could not load one or more CSV files:\n{error}"
            )
            self.statusBar().showMessage("File load error")
            return

        # Execute the stored callback (e.g., plotting or stats generation) safely on the main thread
        if hasattr(self, '_current_on_done') and self._current_on_done:
            self._current_on_done(data_by_path)


    def select_files(self):
        """
        Function enabling user-selected files using QFileDialog.
        :return:
        """
        files, _ = QFileDialog.getOpenFileNames(
            self, "Select files to plot", self.initial_dir,
            "CSV files (*.csv);;All files (*.*)"
        )
        if not files:
            return
        self.file_paths = list(files)
        self.launch_names = list(all_names(self.file_paths))
        self._csv_cache = {}  # selection changed - drop any stale cached DataFrames
        self._water_cache = {}

        # `active_paths` (and the stats panel showing them) reflect whatever was last
        # plotted, which is now a stale selection - previously this was left in place,
        # so "Export Stats" after selecting new files but before re-plotting would
        # silently export stats for the old selection (and try to re-read old files that
        # may no longer exist at that path). Clearing it here means "Export Stats" falls
        # through to its existing "please plot at least one file" guard instead.
        self.active_paths = []
        self.stats_window.hide()
        self.stats_toggle_btn.setChecked(False)

        # Remember this directory so next launch starts here instead of the user
        # having to navigate back to it manually.
        self.initial_dir = os.path.dirname(self.file_paths[0])
        self.settings.setValue("last_directory", self.initial_dir)

        self._refresh_file_listbox()


    def _refresh_file_listbox(self):
        """
        Adds selected .csv filenames to the file display window so the user can see what files have been selected.
        :return:
        """
        # Wipe all existing checkboxes (keep the trailing stretch item)
        while self.file_list_layout.count() > 1:
            item = self.file_list_layout.takeAt(0)
            widget = item.widget()
            if widget is not None:
                widget.deleteLater()
        self.file_checks = []

        # Create one checkbox per file
        for file_path in self.file_paths:
            file_name = os.path.basename(file_path)
            display_name = file_name if len(file_name) <= 38 else file_name[:35] + "..."

            cb = QCheckBox(display_name)
            cb.setChecked(True)
            self.file_list_layout.insertWidget(self.file_list_layout.count() - 1, cb)
            self.file_checks.append(cb)


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

        value = self.confidence_entry.text().strip()
        if not value:
            return 0.95

        try:
            confidence_level = float(value)
            if not (0.0 < confidence_level < 1.0):
                raise ValueError("Confidence level must be between 0 and 1 (exclusive).")
            return confidence_level
        except Exception as e:
            QMessageBox.critical(
                self, "Invalid confidence level",
                f"Please enter a valid confidence value between 0 and 1 (e.g. 0.95).\n\nError: {e}"
            )
            self.statusBar().showMessage("Invalid confidence value")
            return None


    def plot_selected(self):
        """
        Function is called when 'Plot' button is clicked; loads the checked files (in the
        background) and then calls plot_data() to post-process data and display on the graph.
        :return:
        """
        self.active_paths = [p for p, cb in zip(self.file_paths, self.file_checks) if cb.isChecked()]
        if not self.active_paths:
            QMessageBox.warning(self, "No files checked", "Please check at least one file to plot.")
            return

        LC_flag = self.LC_ellipse_box.isChecked()
        sigma_flag = self.sigma_ellipse_box.isChecked()
        confidence_flag = self.confidence_ellipse_box.isChecked()
        confidence_level = self._get_validated_confidence(confidence_flag)
        if confidence_level is None:
            return

        plot_title = self.title_entry.text()
        top_outliers_flag = self.top_outliers_box.isChecked()
        water_landings_flag = self.water_landings_box.isChecked()


        def _do_plot(data_by_path):
            try:
                water_by_path = {p: self._get_water_landing(p) for p in self.active_paths}
                self.outlier_date_summary = plot_data(
                    file_paths=self.active_paths,
                    plot_title=plot_title,
                    fig=self.fig,
                    ax=self.ax,
                    data_by_path=data_by_path,
                    water_by_path=water_by_path,
                    plot_LC_ellipse=LC_flag,
                    plot_sigma_ellipses=sigma_flag,
                    plot_confidence_ellipse=confidence_flag,
                    confidence=confidence_level,
                    plot_top_outliers=top_outliers_flag,
                    plot_water_landings=water_landings_flag
                ) or {}
                self.canvas.draw()
                self._populate_stats_panel()
                self.statusBar().showMessage("Plot updated")
            except Exception as error:
                QMessageBox.critical(
                    self, "Plot error",
                    f"An error occurred while plotting:\n{error}"
                )
                self.statusBar().showMessage("Plot error")

        self._load_csvs_async(self.active_paths, _do_plot)


    def _populate_stats_panel(self):
        """Rebuilds the docked 'Flight Statistics' panel's contents for the current active_paths."""
        old_widget = self.stats_scroll.takeWidget()
        if old_widget is not None:
            old_widget.deleteLater()

        inner = QWidget()
        inner_layout = QVBoxLayout(inner)

        raw_colours, _ = get_colors(self.active_paths)
        colours = [to_hex(c) for c in raw_colours]
        self._sim_param_files = {}

        for i, file_path in enumerate(self.active_paths):
            file_name = os.path.basename(file_path)
            data = self._load_csv(file_path)
            stats = coordinate_stats(data, file_label=file_name)
            display_header = file_name if len(file_name) <= 40 else file_name[:37] + "..."
            self._sim_param_files[i] = None

            # Header
            header_layout = QHBoxLayout()
            swatch = QLabel()
            swatch.setFixedSize(18, 18)
            swatch.setStyleSheet(f"background-color: {colours[i]}; border: 1px solid black;")
            header_label = QLabel(display_header)
            font = header_label.font()
            font.setBold(True)
            font.setItalic(True)
            font.setPointSize(11)
            header_label.setFont(font)
            header_layout.addWidget(swatch)
            header_layout.addWidget(header_label)
            header_layout.addStretch()
            inner_layout.addLayout(header_layout)

            water_count, water_prob, _ = self._get_water_landing(file_path)

            # Stats grid layout
            stats_grid = QGridLayout()
            rows = [
                ("Total Simulations", f"{stats.total_simulations}"),
                ("Mean Apogee", f"{stats.mean_apogee:,} ft"),
                ("Std Dev Apogee", f"{stats.std_apogee:.1f} ft"),
                ("Mean Landing Distance", f"{stats.mean_landing_distance:.1f} NM @ {stats.theta % 360}°"),
                ("Std Dev Landing Dist.", f"{stats.std_landing_distance:.1f} NM"),
                ("Max Landing Distance", f"{stats.max_landing_distance:.1f} NM"),
                ("Avg Landing Coordinates", f"({stats.avg_lat}, {stats.avg_lon})"),
                (f"Accuracy (within {LC_WAIVER_RADIUS_NM} NM)", f"{stats.accuracy_launches * 100:.1f}%"),
                ("Water Landing Probability", f"{water_prob:.2f}% ({water_count} sims)"),
                ("Mean Min Stability", f"{stats.mean_min_stability:.3f}"),
                ("Mean Lateral Velocity", f"{stats.mean_lateral_velocity:.2f} m/s"),
                ("Mean Wind Speed", f"{stats.mean_wind_speed:.2f} kn"),
            ]
            for row_idx, (row_label, value) in enumerate(rows):
                stats_grid.addWidget(QLabel(row_label), row_idx, 0)
                value_label = QLabel(value)
                vfont = value_label.font()
                vfont.setBold(True)
                vfont.setPointSize(11)
                value_label.setFont(vfont)
                stats_grid.addWidget(value_label, row_idx, 1, Qt.AlignmentFlag.AlignRight)
            inner_layout.addLayout(stats_grid)

            # Buttons
            btn_layout = QHBoxLayout()
            sim_label = QLabel("No file")
            sim_label.setObjectName("simFileLabel")

            graph_btn = QPushButton("Plot Wind-Altitude Chart")
            graph_btn.setVisible(False)
            graph_btn.clicked.connect(
                lambda _checked=False, idx=i, fp=file_path: self._run_outlier_graph(idx, fp)
            )


            def _upload_sim(_checked=False, idx=i, lbl=sim_label, gb=graph_btn):
                path, _ = QFileDialog.getOpenFileName(self, "Select Sim Parameters CSV")
                if path:
                    self._sim_param_files[idx] = path
                    sim_name = os.path.basename(path)
                    lbl.setText(sim_name if len(sim_name) <= 22 else sim_name[:19] + "...")
                    lbl.setObjectName("simFileLabelActive")
                    lbl.style().unpolish(lbl)
                    lbl.style().polish(lbl)
                    gb.setVisible(True)

            upload_btn = QPushButton("Upload Sim Params")
            upload_btn.clicked.connect(_upload_sim)

            btn_layout.addWidget(upload_btn)
            btn_layout.addWidget(sim_label)
            btn_layout.addStretch()
            inner_layout.addLayout(btn_layout)
            inner_layout.addWidget(graph_btn)

            if file_path != self.active_paths[-1]:
                separator = QFrame()
                separator.setFrameShape(QFrame.Shape.HLine)
                separator.setFrameShadow(QFrame.Shadow.Sunken)
                inner_layout.addWidget(separator)

        inner_layout.addStretch()
        self.stats_scroll.setWidget(inner)
        self.stats_window.show()
        self.stats_toggle_btn.setChecked(True)


    def _run_outlier_graph(self, i, hist_path):
        """Builds a divided window with overlay plotting and statistical analysis."""
        sim_path = self._sim_param_files[i]

        if i in self._wind_windows and self._wind_windows[i].isVisible():
            self._wind_windows[i].close()

        def _do_graph(data_by_path):
            sim_results = data_by_path[hist_path]
            sim_params = data_by_path[sim_path]

            fig = Figure(figsize=(8, 6))
            win = OutlierWindow(self, fig)
            win.setWindowTitle(f"Outlier Wind Analysis — {os.path.basename(hist_path)}")
            win.resize(1050, 650)
            self._wind_windows[i] = win

            # WA_DeleteOnClose means `win` is actually destroyed (asynchronously) after
            # closing, not just hidden - drop the now-dangling dict entry when that
            # happens. Guarded by identity so this can't accidentally remove a *newer*
            # window that reused the same index `i` in the meantime.
            win.destroyed.connect(
                lambda _obj=None, idx=i, w=win: (
                    self._wind_windows.pop(idx, None) if self._wind_windows.get(idx) is w else None
                )
            )

            # Structure window: Left (Plot 75%), Right (Sidebar 25%)
            main_layout = QHBoxLayout(win)

            plot_widget = QWidget()
            plot_layout = QVBoxLayout(plot_widget)

            stats_group = QGroupBox("Outlier Analytics")
            stats_layout = QVBoxLayout(stats_group)

            main_layout.addWidget(plot_widget, 3)
            main_layout.addWidget(stats_group, 1)

            outliers, summary = analyze_outlier_winds(sim_results, sim_params)

            # Wind-based Outlier Stats
            if summary.get("total_outliers", 0) > 0:
                stats_data = [
                    ("Total Outliers Detected", f"{summary['total_outliers']} flights"),
                    ("Worst Shear Altitude", f"{summary['overall_max_speed_alt']:,} m"),
                    ("Peak Avg Layer Speed", f"{summary['overall_max_speed']:.1f} kn"),
                ]
                for stat_label, val in stats_data:
                    label_widget = QLabel(stat_label)
                    label_widget.setObjectName("mutedLabel")
                    stats_layout.addWidget(label_widget)
                    value_widget = QLabel(val)
                    vfont = value_widget.font()
                    vfont.setBold(True)
                    vfont.setPointSize(14)
                    value_widget.setFont(vfont)
                    stats_layout.addWidget(value_widget)
            else:
                info_label = QLabel("No outliers >10 NM found.")
                info_label.setStyleSheet("font-style: italic;")
                stats_layout.addWidget(info_label, alignment=Qt.AlignmentFlag.AlignHCenter)

            # Top-Outlier Landing Dates (from the main dispersion plot's outlier highlighting)
            outlier_dates = self.outlier_date_summary.get(hist_path, {})
            if outlier_dates:
                separator = QFrame()
                separator.setFrameShape(QFrame.Shape.HLine)
                stats_layout.addWidget(separator)

                title_label = QLabel(f"Top {TOP_OUTLIERS_COUNT} Outlier Landing Dates")
                title_label.setObjectName("mutedLabel")
                stats_layout.addWidget(title_label)

                dates_grid = QGridLayout()
                for row_idx, (date, count) in enumerate(sorted(outlier_dates.items())):
                    dates_grid.addWidget(QLabel(date), row_idx, 0)
                    count_label = QLabel(str(count))
                    cfont = count_label.font()
                    cfont.setBold(True)
                    count_label.setFont(cfont)
                    dates_grid.addWidget(count_label, row_idx, 1, Qt.AlignmentFlag.AlignRight)
                stats_layout.addLayout(dates_grid)

            stats_layout.addStretch()

            # Render Plot
            ax = fig.add_subplot(111)
            plot_outlier_analysis(ax, outliers, summary)

            canvas = FigureCanvasQTAgg(fig)
            toolbar = NavigationToolbar2QT(canvas, win)
            plot_layout.addWidget(canvas)
            plot_layout.addWidget(toolbar)

            win.show()

        self._load_csvs_async([hist_path, sim_path], _do_graph)


    def save_file(self):
        """
        Function is called when 'Save Plot' button is clicked; calls save_plot() to save the figure.
        :return:
        """
        if not self.file_paths:
            QMessageBox.warning(self, "No files", "Please select files before plotting.")
            return

        LC_flag = self.LC_ellipse_box.isChecked()
        sigma_flag = self.sigma_ellipse_box.isChecked()
        confidence_flag = self.confidence_ellipse_box.isChecked()

        confidence_level = self._get_validated_confidence(confidence_flag)
        if confidence_level is None:
            return

        plot_title = self.title_entry.text()
        default_name = make_safe_filename(plot_title).name

        save_path, _ = QFileDialog.getSaveFileName(
            self, "Save plot as...", default_name,
            "PNG image (*.png);;JPEG image (*.jpg *.jpeg);;All files (*.*)"
        )

        # If user canceled the dialog, return
        if not save_path:
            self.statusBar().showMessage("Save cancelled")
            return

        # Confirm overwrite if file exists
        if os.path.exists(save_path):
            reply = QMessageBox.question(
                self, "Confirm overwrite",
                f"File already exists:\n{save_path}\n\nOverwrite?",
                QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No
            )
            if reply != QMessageBox.StandardButton.Yes:
                self.statusBar().showMessage("Save cancelled")
                return


        def _do_save(data_by_path):
            try:
                water_by_path = {p: self._get_water_landing(p) for p in self.file_paths}
                save_plot(
                    file_paths=self.file_paths,
                    plot_title=plot_title,
                    output_path=save_path,
                    data_by_path=data_by_path,
                    water_by_path=water_by_path,
                    plot_LC_ellipse=LC_flag,
                    plot_sigma_ellipses=sigma_flag,
                    plot_confidence_ellipse=confidence_flag,
                    confidence=confidence_level,
                    plot_top_outliers=self.top_outliers_box.isChecked(),
                    plot_water_landings=self.water_landings_box.isChecked()
                )
                self.statusBar().showMessage(f"Plot saved: '{plot_title}'")
            except Exception as error:
                QMessageBox.critical(
                    self, "Plot error",
                    f"An error occurred while plotting:\n{error}"
                )
                self.statusBar().showMessage("Plot error")

        self._load_csvs_async(self.file_paths, _do_save)


    def export_stats(self):
        """
        Exports the per-file stats for the current active_paths (i.e. whatever's in the stats
        panel from the last Plot) to a .csv file.
        :return:
        """
        if not self.active_paths:
            QMessageBox.warning(self, "No stats", "Please plot at least one file before exporting stats.")
            return

        save_path, _ = QFileDialog.getSaveFileName(
            self, "Export stats as...", "flight_stats.csv", "CSV files (*.csv);;All files (*.*)"
        )
        if not save_path:
            self.statusBar().showMessage("Export cancelled")
            return

        def _do_export(data_by_path):
            try:
                rows = []
                for file_path in self.active_paths:
                    file_name = os.path.basename(file_path)
                    data = data_by_path[file_path]
                    stats = coordinate_stats(data, file_label=file_name)
                    rows.append({
                        "File": file_name,
                        "Total Simulations": stats.total_simulations,
                        "Mean Apogee (ft)": stats.mean_apogee,
                        "Std Dev Apogee (ft)": stats.std_apogee,
                        "Mean Landing Distance (NM)": stats.mean_landing_distance,
                        "Std Dev Landing Distance (NM)": stats.std_landing_distance,
                        "Max Landing Distance (NM)": stats.max_landing_distance,
                        "Avg Landing Latitude": stats.avg_lat,
                        "Avg Landing Longitude": stats.avg_lon,
                        f"Accuracy (within {LC_WAIVER_RADIUS_NM} NM) %": stats.accuracy_launches * 100,
                        "Mean Min Stability": stats.mean_min_stability,
                        "Mean Lateral Velocity (m/s)": stats.mean_lateral_velocity,
                        "Mean Wind Speed (kn)": stats.mean_wind_speed,
                    })

                pd.DataFrame(rows).to_csv(save_path, index=False)
                self.statusBar().showMessage(f"Stats exported: '{os.path.basename(save_path)}'")
            except Exception as error:
                QMessageBox.critical(
                    self, "Export error",
                    f"An error occurred while exporting stats:\n{error}"
                )
                self.statusBar().showMessage("Export error")

        # Previously this read each CSV synchronously on the GUI thread via self._load_csv,
        # unlike every other CSV-touching action in the app - meaning a cache miss here
        # (e.g. right after a fresh selection) could freeze the UI on a large file. Routing
        # through _load_csvs_async keeps this consistent with plot_selected/save_file.
        self._load_csvs_async(self.active_paths, _do_export)


    def clear_all(self):
        """
        Clears all inputs and plots, resets to default values.
        :return:
        """
        self.file_paths = []
        while self.file_list_layout.count() > 1:
            item = self.file_list_layout.takeAt(0)
            widget = item.widget()
            if widget is not None:
                widget.deleteLater()
        self.file_checks = []
        self.active_paths = []
        self._csv_cache = {}
        self._water_cache = {}
        self._sim_param_files = {}

        # Close any open outlier-wind windows rather than leaving them floating around
        # showing stats for files that "Clear" just dropped from the app entirely.
        for win in list(self._wind_windows.values()):
            if win.isVisible():
                win.close()
        self._wind_windows = {}

        self.ax.clear()
        self.canvas.draw()

        old_widget = self.stats_scroll.takeWidget()
        if old_widget is not None:
            old_widget.deleteLater()
        self.stats_window.hide()
        self.stats_toggle_btn.setChecked(False)

        self.statusBar().showMessage("Cleared files and plot")
