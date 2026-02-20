import tkinter as tk
from tkinter import ttk, filedialog, messagebox
import os
import tempfile

import database as db
from src.parameters import (compute_properties, db_row, FILTER_COLUMNS, FILTER_DESCRIPTIONS)
from src.parser import parse_sdf

# RDKit drawing for 2D structure images 
from rdkit import Chem
from rdkit.Chem import Draw, AllChem
DRAW_AVAILABLE = True

# Pillow for displaying images in Tkinter ---
from PIL import Image, ImageTk
PIL_AVAILABLE = True


class DatabaseGUI(tk.Tk):
    """
    GUI API for the Chemical Compound Database.
    Attributes:
        conn:   SQLite connection
        all_data:   list of all compound rows.
        current_data:   Currently displayed rows (may be filtered).
        sort_column:    Active sort column name.
        sort_ascending: Current sort direction.
        active_filter:  Name of active filter, or None.
    """
    # Table columns: (db_key, display_header, width)
    COLUMNS = [
        ("compound_id", "ID",   50),
        ("name",    "CompoundName", 120),
        ("smiles",  "SMILES",   180),
        ("formula", "Formula",  100),
        ("molecular_weight",    "MW (Da)",  50),
        ("logp",    "LogP",     50),
        ("hbd", "HBDonors",  50),
        ("hba", "HBAcceptors",  50),
        ("polar_surface_area",  "PSA",  50),
        ("rotatable_bonds", "RotatableBonds", 50),
        ("aromatic_rings",  "AromaticRings",    50),
        ("heavy_atom_count",    "HeavyAtoms",   50),
        ("qed", "QED",  50),
    ]

    def __init__(self):
        """Initialise the GUI window, database, and all widgets."""
        super().__init__()
        self.title("Chemical Compound Database — Drug Discovery")
        self.geometry("1420x760")
        self.minsize(1050, 600)

        # Database setup 
        script_dir = os.path.dirname(os.path.abspath(__file__))
        db_path = os.path.join(script_dir, "compounds.db")
        schema_path = os.path.join(script_dir, "database.sql")
        self.conn = db.connect(db_path)
        db.init_schema(self.conn, schema_path)

        # State 
        self.all_data = []
        self.current_data = []
        self.sort_column = "compound_id"
        self.sort_ascending = True
        self.active_filter = None

        #  Build GUI 
        self._create_menu()
        self._create_toolbar()
        self._create_filter_bar()
        self._create_main_area()
        self._create_status_bar()

        # Load existing data 
        self._refresh_data()
        self.protocol("WM_DELETE_WINDOW", self._on_close)

# GUI Construction
    def _create_menu(self):
            """Top menu bar."""
            menubar = tk.Menu(self)
            self.config(menu=menubar)

            file_menu = tk.Menu(menubar, tearoff=0)
            file_menu.add_command(label="Load SDF File", command=self._load_sdf)
            file_menu.add_command(label="Clear Database", command=self._clear_db)
            file_menu.add_separator()
            file_menu.add_command(label="Exit", command=self._on_close)
            menubar.add_cascade(label="File", menu=file_menu)

            help_menu = tk.Menu(menubar, tearoff=0)
            help_menu.add_command(label="Filter Criteria", command=self._show_filter_help)
            help_menu.add_command(label="About", command=self._show_about)
            menubar.add_cascade(label="Help", menu=help_menu)

    def _create_toolbar(self):
            """Toolbar with Load/Clear buttons and compound count."""
            bar = ttk.Frame(self, padding=5)
            bar.pack(side=tk.TOP, fill=tk.X)

            ttk.Button(bar, text="Load SDF File", command=self._load_sdf).pack(side=tk.LEFT, padx=5)
            ttk.Button(bar, text="Clear Database", command=self._clear_db).pack(side=tk.LEFT, padx=5)

            self.lbl_count = ttk.Label(bar, text="Compounds: 0",font=("TkDefaultFont", 10, "bold"))
            self.lbl_count.pack(side=tk.RIGHT, padx=10)

    def _create_filter_bar(self):
        """
        Filter buttons for all 5 drug discovery criteria
        Filter pass/fail is pre-computed during import and stored in the parameters table;filtering is SQL WHERE clause
        """
        frame = ttk.LabelFrame(self, text="Drug Discovery Filters", padding=5)
        frame.pack(side=tk.TOP, fill=tk.X, padx=5, pady=(0, 5))

        ttk.Button(frame, text="Show All", command=self._show_all).pack(side=tk.LEFT, padx=3)
        ttk.Separator(frame, orient=tk.VERTICAL).pack(side=tk.LEFT, fill=tk.Y, padx=5)

        self.filter_buttons = {}
        for name in FILTER_COLUMNS:
            btn = ttk.Button(frame, text=name,
                             command=lambda n=name: self._apply_filter(n))
            btn.pack(side=tk.LEFT, padx=3)
            self.filter_buttons[name] = btn

        self.lbl_filter = ttk.Label(frame, text="No filter active",foreground="grey")
        self.lbl_filter.pack(side=tk.RIGHT, padx=10)