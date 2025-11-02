import os
import tkinter as tk
import matplotlib.pyplot as plt

from dispersion_backend import plot_data
from dispersion_backend import save_plot
from dispersion_backend import set_default_style
from dispersion_backend import make_safe_filename
from tkinter import ttk, filedialog, messagebox
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk


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

        # Left column stacked buttoms
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

        self.file_listbox = tk.Listbox(list_frame, height=12, width=40, selectmode=tk.EXTENDED)
        self.file_listbox.grid(row=0, column=0, sticky="nsew")

        scrollbar = ttk.Scrollbar(list_frame, orient="vertical", command=self.file_listbox.yview)
        scrollbar.grid(row=0, column=1, sticky="ns")
        self.file_listbox.config(yscrollcommand=scrollbar.set)

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
        self._refresh_file_listbox()
        self.status_var.set(f"Selected {len(self.file_paths)} files")

    def _refresh_file_listbox(self):
        """
        Adds selected .csv filenames to the file display window so the user can see what files have been selected.
        :return:
        """
        self.file_listbox.delete(0, tk.END)
        for file_path in self.file_paths:
            self.file_listbox.insert(tk.END, os.path.basename(file_path))

    def plot_selected(self):
        """
        Function is called when 'Plot data' button is clicked; calls plot_data() to post-process data and display on
        the graph.
        :return:
        """
        if not self.file_paths:
            messagebox.showwarning("No files", "Please select files before plotting.")
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

            plot_data(
                file_paths=self.file_paths,
                plot_title=self.plot_title.get(),
                ax=self.ax,
                plot_LC_ellipse=LC_flag,
                plot_sigma_ellipses=sigma_flag,
                plot_confidence_ellipse=confidence_flag,
                confidence=confidence_level
            )
            self.fig.tight_layout()
            self.canvas.draw()
            self.status_var.set("Plotted selected files")
        except Exception as error:
            messagebox.showerror("Plot error", f"An error occurred while plotting:\n{error}")
            self.status_var.set("Plot error")


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
        self.file_listbox.delete(0, tk.END)
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
