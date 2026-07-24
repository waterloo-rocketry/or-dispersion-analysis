import sys

from PySide6.QtWidgets import QApplication

from app import FilePlotApp


def main():
    app = QApplication(sys.argv)

    # Root paths for testing environments
    initial_dir = "/Users/luca/PycharmProjects/Rocketry Dispersion Zone Analysis/Plugin Exports/"

    window = FilePlotApp(initial_dir=initial_dir)
    window.resize(1200, 800)
    window.show()

    sys.exit(app.exec())


if __name__ == "__main__":
    main()