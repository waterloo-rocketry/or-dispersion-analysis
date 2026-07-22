import tkinter as tk
from app import FilePlotApp


def main():
    root = tk.Tk()

    # Root paths for testing environments
    initial_dir = "/Users/luca/PycharmProjects/Rocketry Dispersion Zone Analysis/Plugin Exports/"

    app = FilePlotApp(root, initial_dir=initial_dir)
    root.geometry("1200x800")
    root.mainloop()


if __name__ == "__main__":
    main()
