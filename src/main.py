import sys
import os

# Ensure src directory is on the import path
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from GUI import DatabaseGUI


def main():
    """
    Launches the Tkinter GUI application.
    """
    print("Chemical Compound Database — Drug Discovery")
    print("=" * 48)
    print("Starting GUI...")
    print()

    app = DatabaseGUI()
    app.mainloop()


if __name__ == "__main__":
    main()