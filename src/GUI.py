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


class CompoundDatabaseGUI(tk.Tk):
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
