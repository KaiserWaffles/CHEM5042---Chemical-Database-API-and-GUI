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

    def _create_main_area(self):
        """Compound table (left) and detail panel (right)."""
        pane = ttk.PanedWindow(self, orient=tk.HORIZONTAL)
        pane.pack(side=tk.TOP, fill=tk.BOTH, expand=True, padx=5, pady=5)

        # Table
        tbl_frame = ttk.Frame(pane)
        pane.add(tbl_frame, weight=3)

        col_ids = [c[0] for c in self.COLUMNS]
        self.tree = ttk.Treeview(tbl_frame, columns=col_ids, show="headings", selectmode="browse")
        for col_name, header, width in self.COLUMNS:
            self.tree.heading(col_name, text=header,command=lambda c=col_name: self._sort_by(c))
            self.tree.column(col_name, width=width, minwidth=35)

        vsb = ttk.Scrollbar(tbl_frame, orient=tk.VERTICAL, command=self.tree.yview)
        hsb = ttk.Scrollbar(tbl_frame, orient=tk.HORIZONTAL, command=self.tree.xview)
        self.tree.configure(yscrollcommand=vsb.set, xscrollcommand=hsb.set)
        self.tree.grid(row=0, column=0, sticky="nsew")
        vsb.grid(row=0, column=1, sticky="ns")
        hsb.grid(row=1, column=0, sticky="ew")
        tbl_frame.grid_rowconfigure(0, weight=1)
        tbl_frame.grid_columnconfigure(0, weight=1)

        self.tree.bind("<<TreeviewSelect>>", self._on_select)

        # Detail panel
        detail = ttk.LabelFrame(pane, text="Compound Details", padding=10)
        pane.add(detail, weight=1)

        self.detail_text = tk.Text(detail, wrap=tk.WORD, height=16, width=34, state=tk.DISABLED, font=("Courier", 10))
        self.detail_text.pack(fill=tk.X, pady=(0, 8))

        self.filter_frame = ttk.LabelFrame(detail, text="Filter Results", padding=5)
        self.filter_frame.pack(fill=tk.X, pady=(0, 8))

        self.mol_frame = ttk.LabelFrame(detail, text="2D Structure", padding=5)
        self.mol_frame.pack(fill=tk.BOTH, expand=True)
        self.mol_label = ttk.Label(self.mol_frame, text="Select a compound\nto view structure", anchor=tk.CENTER)
        self.mol_label.pack(fill=tk.BOTH, expand=True)

    def _create_status_bar(self):
        """Bottom status bar."""
        self.status = ttk.Label(self, text="Ready. Load an SDF file to begin.",
                                relief=tk.SUNKEN, anchor=tk.W, padding=3)
        self.status.pack(side=tk.BOTTOM, fill=tk.X)

    # loading data
    def _load_sdf(self):
        """
        Load an SDF file, parse, calculate parameters and store in DB parse
            - User selects .sdf file via dialog
            - parser.parse_sdf() reads molecules with RDKit
            - parameters.compute_properties() calculates all parameters
            - parameters.db_row() converts to dict for database
            - database.insert_compound() and insert_properties() store rows
            - GUI table refreshes
        """
        filepath = filedialog.askopenfilename(title="Select SDF File", filetypes=[("SDF files", "*.sdf"), ("All files", "*.*")]
        )
        if not filepath:
            return
        try:
            self.status.config(text=f"Loading {os.path.basename(filepath)}...")
            self.update()

            # Clear existing data to avoid duplicates
            db.clear_database(self.conn)

            # Parse and import within a usage
            count = 0
            errors = 0
            for mol_data in parse_sdf(filepath):
                try:
                    # Calculate all parameters (parameters.py)
                    props = compute_properties(mol_data.mol)

                    # Insert compound identity (database.py)
                    cid = db.insert_compound(
                        self.conn,
                        name=mol_data.name,
                        smiles=mol_data.smiles,
                        formula=mol_data.formula
                    )

                    # Insert parameters + filter flags (database.py)
                    db.insert_properties(self.conn, cid, db_row(props))
                    count += 1

                except Exception as e:
                    print(f"  WARNING: Skipping {mol_data.name}: {e}")
                    errors += 1

            self.conn.commit()
            # Refresh GUI
            self._refresh_data()
            self.active_filter = None
            self.lbl_filter.config(text="No filter active", foreground="grey")

            msg = f"Loaded {count} compounds from {os.path.basename(filepath)}"
            if errors:
                msg += f" ({errors} skipped)"
            self.status.config(text=msg)
            messagebox.showinfo("Import Complete", msg)

        except Exception as e:
            try:
                self.conn.rollback()
            except Exception:
                pass
            messagebox.showerror("Import Error", f"Failed:\n\n{e}")
            self.status.config(text="Error loading file.")

    def _clear_db(self):
        """Clear database after confirmation."""
        if db.get_compound_count(self.conn) == 0:
            messagebox.showinfo("Empty", "Database is already empty.")
            return
        if messagebox.askyesno("Confirm", "Remove all compounds?"):
            db.clear_database(self.conn)
            self._refresh_data()
            self.active_filter = None
            self.lbl_filter.config(text="No filter active", foreground="grey")
            self.status.config(text="Database cleared.")
        
    def _refresh_data(self):
        """Reload all data from database and update table."""
        rows = db.get_all_compounds(self.conn)
        # Convert sqlite3.Row to plain tuples for easier sorting
        self.all_data = [tuple(r) for r in rows]
        self.current_data = list(self.all_data)
        self._populate_table(self.current_data)
        self.lbl_count.config(text=f"Compounds: {len(self.all_data)}")

    def _populate_table(self, data):
        """Clear and fill the Treeview with data rows."""
        for item in self.tree.get_children():
            self.tree.delete(item)

        # Only show the columns we defined (first 13 fields)
        num_cols = len(self.COLUMNS)
        for row in data:
            vals = []
            for i in range(num_cols):
                v = row[i] if i < len(row) else ""
                if v is None:
                    vals.append("N/A")
                elif isinstance(v, float):
                    vals.append(f"{v:.2f}")
                else:
                    vals.append(str(v))
            self.tree.insert("", tk.END, values=vals)

        total = len(self.all_data)
        showing = len(data)
        if self.active_filter:
            self.status.config(
                text=f"Showing {showing}/{total} (Filter: {self.active_filter})")
        else:
            self.status.config(text=f"Showing all {showing} compounds")
        
    # Sorting feature
    def _sort_by(self, col_name):
        """
        Sort table by clicked column header; click again to reverse.
        """
        if self.sort_column == col_name:
            self.sort_ascending = not self.sort_ascending
        else:
            self.sort_column = col_name
            self.sort_ascending = True

        col_names = [c[0] for c in self.COLUMNS]
        idx = col_names.index(col_name)

        def key(row):
            v = row[idx] if idx < len(row) else None
            if v is None:
                return (1, "")
            if isinstance(v, str):
                return (0, v.lower())
            return (0, v)

        self.current_data.sort(key=key, reverse=not self.sort_ascending)
        self._populate_table(self.current_data)

        arrow = "▲" if self.sort_ascending else "▼"
        for cn, hdr, _ in self.COLUMNS:
            self.tree.heading(cn, text=(hdr + arrow) if cn == col_name else hdr)
    
    # Filtering feature
    def _apply_filter(self, filter_name):
        """
        Filter compounds using pre-computed pass/fail columns.
        Filters are calculated at import time and stored in the
        parameters table: column check, works instantly even with 1,000,000 compounds.
        """
        if not self.all_data:
            messagebox.showinfo("No Data", "Load an SDF file first.")
            return

        col = FILTER_COLUMNS[filter_name]

        # Index of the filter column in our query results
        # Query order: compound_id(0), name(1), smiles(2), formula(3),
        #   mw(4), logp(5), hbd(6), hba(7), psa(8), rotb(9), rings(10),
        #   hatom(11), qed(12), ro5_viol(13), ro5_pass(14), ro3_pass(15),
        #   leadlike(16), veber(17), bioavail(18)
        col_map = {
            "ro5_pass": 14, "ro3_pass": 15, "leadlike_pass": 16,
            "veber_pass": 17, "bioavailability_pass": 18,
        }
        idx = col_map[col]

        self.current_data = [r for r in self.all_data if r[idx] == 1]
        self.active_filter = filter_name

        passed = len(self.current_data)
        total = len(self.all_data)
        pct = (100 * passed / total) if total > 0 else 0

        self.lbl_filter.config(
            text=f"{filter_name}: {passed}/{total} pass ({pct:.1f}%)",
            foreground="dark green" if passed > 0 else "red"
        )
        self._populate_table(self.current_data)

    # Detail view feature
    def _on_select(self, event):
        """Show compound details + filter summary + 2D structure."""
        sel = self.tree.selection()
        if not sel:
            return
        vals = self.tree.item(sel[0])["values"]
        if not vals:
            return

        # Find full row in cached data
        cid = vals[0]
        row = None
        for r in self.all_data:
            if str(r[0]) == str(cid):
                row = r
                break
        if row is None:
            return

        # Properties text
        self.detail_text.config(state=tk.NORMAL)
        self.detail_text.delete("1.0", tk.END)

        def fmt(v, dp=2):
            if v is None: return "N/A"
            if isinstance(v, float): return f"{v:.{dp}f}"
            return str(v)
        
        lines = [
            f"ID:   {row[0]}",
            f"Name: {row[1]}",
            f"SMILES:   {row[2]}",
            f"Formula:  {row[3]}",
            "-" * 35,
            f"MW (Da):  {fmt(row[4])}",
            f"LogP: {fmt(row[5])}",
            f"HBD:  {row[6]}",
            f"HBA:  {row[7]}",
            f"PSA (A^2):    {fmt(row[8])}",
            f"Rotatable Bonds:  {row[9]}",
            f"Aromatic Rings:   {row[10]}",
            f"Heavy Atoms:  {row[11]}",
            f"QED:  {fmt(row[12], 3)}",
            "-" * 35,
            f"Ro5 Violations:   {row[13]}",
        ]
        self.detail_text.insert(tk.END, "\n".join(lines))
        self.detail_text.config(state=tk.DISABLED)

        # Filter pass/fail labels
        for w in self.filter_frame.winfo_children():
            w.destroy()

        # row indices: ro5_pass=14, ro3=15, leadlike=16, veber=17, bioavail=18
        filter_results = [
            ("Lipinski Ro5",    row[14]),
            ("Rule of 3 (Fragment)",    row[15]),
            ("Lead-likeness",   row[16]),
            ("Veber Rules", row[17]),
            ("Bioavailability", row[18]),
        ]
        for fname, val in filter_results:
            passed = (val == 1)
            colour = "green" if passed else "red"
            symbol = "PASS" if passed else "FAIL"
            ttk.Label(self.filter_frame, text=f"  [{symbol}]  {fname}",
                     foreground=colour).pack(anchor=tk.W, pady=1)
        # 2D Structure
        self._show_structure(row[2])

    # View 2D structure feature
    def _show_structure(self, smiles):
        """Render 2D molecular structure from SMILES."""
        if not DRAW_AVAILABLE or not PIL_AVAILABLE:
            self.mol_label.config(
                text="Install RDKit + Pillow\nfor 2D structure view", image="")
            return
        if not smiles:
            self.mol_label.config(text="No SMILES", image="")
            return
        try:
            mol = Chem.MolFromSmiles(str(smiles))
            if mol is None:
                self.mol_label.config(text="Invalid SMILES", image="")
                return
            AllChem.Compute2DCoords(mol)
            img = Draw.MolToImage(mol, size=(280, 280))
            photo = ImageTk.PhotoImage(img)
            self.mol_label.config(image=photo, text="")
            self.mol_label.image = photo
        except Exception as e:
            self.mol_label.config(text=f"Error: {str(e)[:50]}", image="")
