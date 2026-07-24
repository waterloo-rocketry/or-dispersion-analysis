import os

from PySide6.QtCore import Qt
from PySide6.QtWidgets import (
    QMainWindow, QWidget, QVBoxLayout, QHBoxLayout, QGridLayout,
    QPushButton, QLabel, QLineEdit, QCheckBox, QGroupBox, QScrollArea,
    QFrame, QFileDialog, QMessageBox, QDialog
)

from matplotlib.figure import Figure
from matplotlib.colors import to_hex
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg, NavigationToolbar2QT

from data_engine import all_names, coordinate_stats, analyze_outlier_winds
from plotting import (get_colors, make_safe_filename, plot_data, save_plot, set_default_style,
                      plot_outlier_analysis, TOP_OUTLIERS_COUNT)


class OutlierWindow(QDialog):
    """
    Non-modal window for the per-file outlier wind analysis plot. Clears its
    matplotlib Figure on close (mirrors the original Tkinter
    WM_DELETE_WINDOW -> fig.clear() handler) so unused figures don't linger.
    """
    def __init__(self, parent, fig):
        super().__init__(parent)
        self._fig = fig

    def closeEvent(self, event):
        self._fig.clear()
        super().closeEvent(event)


class FilePlotApp(QMainWindow):
    def __init__(self, initial_dir=None):
        """
        Initialize class.
        :param initial_dir: user-specified initial directory, optional - returns root directory
        """
        super().__init__()
        self.setWindowTitle("Waterloo Rocketry Dispersion Zone Analysis")
        self.initial_dir = initial_dir or "."
        self.file_paths = []
        self.launch_names = []
        self.file_checks = []
        self.active_paths = []
        self.outlier_date_summary = {}
        self._sim_param_files = {}
        self._stats_win = None
        self._wind_windows = {}
        self._build_ui()
        self._setup_matplotlib()
        self.resize(1200, 800)

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
        for btn in (self.select_btn, self.plot_btn, self.clear_btn, self.save_file_btn):
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
        self.top_outliers_box = QCheckBox(f"Highlight Top {TOP_OUTLIERS_COUNT} Outliers")
        self.top_outliers_box.setChecked(True)

        for box in (self.LC_ellipse_box, self.sigma_ellipse_box,
                    self.confidence_ellipse_box, self.top_outliers_box):
            checkbox_layout.addWidget(box)

        self.confidence_container = QWidget()
        conf_layout = QHBoxLayout(self.confidence_container)
        conf_layout.setContentsMargins(8, 6, 0, 0)
        self.confidence_entry = QLineEdit("0.95")
        self.confidence_entry.setFixedWidth(60)
        conf_layout.addWidget(QLabel("Confidence level (e.g. 0.95):"))
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
        for box in (self.LC_ellipse_box, self.sigma_ellipse_box, self.top_outliers_box):
            box.toggled.connect(self._on_plot_option_changed)
        self.confidence_ellipse_box.toggled.connect(self._toggle_confidence_entry)
        self.title_entry.textChanged.connect(self._on_plot_option_changed)

        self.select_btn.clicked.connect(self.select_files)
        self.plot_btn.clicked.connect(self.plot_selected)
        self.clear_btn.clicked.connect(self.clear_all)
        self.save_file_btn.clicked.connect(self.save_file)

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
        Function is called when 'Plot' button is clicked; calls plot_data() to post-process data and display on
        the graph.
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

        self.outlier_date_summary = plot_data(
            file_paths=self.active_paths,
            plot_title=self.title_entry.text(),
            fig=self.fig,
            ax=self.ax,
            plot_LC_ellipse=LC_flag,
            plot_sigma_ellipses=sigma_flag,
            plot_confidence_ellipse=confidence_flag,
            confidence=confidence_level,
            plot_top_outliers=self.top_outliers_box.isChecked()
        ) or {}
        # self.fig.tight_layout()
        self.canvas.draw()
        self._show_stats_window()

    def _show_stats_window(self):
        if self._stats_win is not None:
            self._stats_win.close()

        self._stats_win = QDialog(self)
        self._stats_win.setWindowTitle("General Flight Statistics")
        self._stats_win.resize(420, 640)

        outer_layout = QVBoxLayout(self._stats_win)
        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        inner = QWidget()
        inner_layout = QVBoxLayout(inner)
        scroll.setWidget(inner)
        outer_layout.addWidget(scroll)

        raw_colours, _ = get_colors(self.active_paths)
        colours = [to_hex(c) for c in raw_colours]
        self._sim_param_files = {}

        for i, file_path in enumerate(self.active_paths):
            stats = coordinate_stats(file_path)
            file_name = os.path.basename(file_path)
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

            # Stats grid layout
            stats_grid = QGridLayout()
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
                stats_grid.addWidget(QLabel(label), row_idx, 0)
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
            sim_label.setStyleSheet("color: grey; font-size: 8pt;")

            graph_btn = QPushButton("Plot Wind-Altitude Chart")
            graph_btn.setVisible(False)
            graph_btn.clicked.connect(
                lambda _checked=False, idx=i, fp=file_path: self._run_outlier_graph(idx, fp)
            )

            def _upload_sim(_checked=False, idx=i, lbl=sim_label, gb=graph_btn):
                path, _ = QFileDialog.getOpenFileName(self._stats_win, "Select Sim Parameters CSV")
                if path:
                    self._sim_param_files[idx] = path
                    sim_name = os.path.basename(path)
                    lbl.setText(sim_name if len(sim_name) <= 22 else sim_name[:19] + "...")
                    lbl.setStyleSheet("color: white; font-size: 8pt;")
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
        self._stats_win.show()

    def _run_outlier_graph(self, i, hist_path):
        """Builds a divided window with overlay plotting and statistical analysis."""
        sim_path = self._sim_param_files[i]

        if i in self._wind_windows and self._wind_windows[i].isVisible():
            self._wind_windows[i].close()

        fig = Figure(figsize=(8, 6))
        win = OutlierWindow(self, fig)
        win.setWindowTitle(f"Outlier Wind Analysis — {os.path.basename(hist_path)}")
        win.resize(1050, 650)
        self._wind_windows[i] = win

        # Structure window: Left (Plot 75%), Right (Sidebar 25%)
        main_layout = QHBoxLayout(win)

        plot_widget = QWidget()
        plot_layout = QVBoxLayout(plot_widget)

        stats_group = QGroupBox("Outlier Analytics")
        stats_layout = QVBoxLayout(stats_group)

        main_layout.addWidget(plot_widget, 3)
        main_layout.addWidget(stats_group, 1)

        outliers, summary = analyze_outlier_winds(hist_path, sim_path)

        # Wind-based Outlier Stats
        if summary.get("total_outliers", 0) > 0:
            stats_data = [
                ("Total Outliers Detected", f"{summary['total_outliers']} flights"),
                ("Worst Shear Altitude", f"{summary['overall_max_speed_alt']:,} m"),
                ("Peak Avg Layer Speed", f"{summary['overall_max_speed']:.1f} kn"),
            ]
            for label, val in stats_data:
                label_widget = QLabel(label)
                label_widget.setStyleSheet("color: grey;")
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
            title_label.setStyleSheet("color: grey;")
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

    def save_file(self):
        """
        Function is called when 'Save Plot' button is clicked; calls save_plot() to save the figure.
        :return:
        """
        if not self.file_paths:
            QMessageBox.warning(self, "No files", "Please select files before plotting.")
            return

        try:
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

            save_plot(
                file_paths=self.file_paths,
                plot_title=plot_title,
                output_path=save_path,
                plot_LC_ellipse=LC_flag,
                plot_sigma_ellipses=sigma_flag,
                plot_confidence_ellipse=confidence_flag,
                confidence=confidence_level,
                plot_top_outliers=self.top_outliers_box.isChecked()
            )
            self.statusBar().showMessage(f"Plot saved: '{plot_title}'")

        except Exception as error:
            QMessageBox.critical(
                self, "Plot error",
                f"An error occurred while plotting:\n{error}"
            )
            self.statusBar().showMessage("Plot error")

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
        self.ax.clear()
        # self.fig.tight_layout()
        self.canvas.draw()
        self.statusBar().showMessage("Cleared files and plot")
