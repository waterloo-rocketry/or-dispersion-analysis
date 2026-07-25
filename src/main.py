import sys

from PySide6.QtWidgets import QApplication

from app import FilePlotApp, APP_STYLESHEET


def main():
    app = QApplication(sys.argv)
    app.setStyleSheet(APP_STYLESHEET)

    # No hardcoded path: FilePlotApp remembers the last directory you picked files
    # from (via QSettings) and reopens there next time, defaulting to your home
    # directory on first run. Pass initial_dir="..." here if you want to override that.
    window = FilePlotApp()
    window.resize(1200, 800)
    window.show()

    sys.exit(app.exec())


if __name__ == "__main__":
    main()
