import sqlite3
from pathlib import Path
from typing import Dict

def connect(db_path: str | Path) -> sqlite3.Connection:
    """
    Create SQLite connection with pragmatic settings;
    WAL improves concurrent reads (GUI) while writing during import.
    """
    conn = sqlite3.connect(str(db_path))
    conn.row_factory = sqlite3.Row
    conn.execute("PRAGMA foreign_keys = ON;")
    conn.execute("PRAGMA journal_mode = WAL;")
    conn.execute("PRAGMA synchronous = NORMAL;")
    return conn

def init_schema(conn: sqlite3.Connection, schema_path: str | Path) -> None:
    """Create tables and indexes from schema file."""
    sql = Path(schema_path).read_text(encoding="utf-8")
    conn.executescript(sql)
    conn.commit()

def insert_compound(conn: sqlite3.Connection, name: str, smiles: str, formula: Optional[str]) -> int:
    """Insert compound identity row; returns compound_id."""
    cur = conn.execute(
        "INSERT INTO compound(name, smiles, formula, source_file) VALUES (?, ?, ?);",
        (name, smiles, formula),
    )
    return int(cur.lastrowid)


def insert_properties(conn: sqlite3.Connection, compound_id: int, props: Dict[str, object]) -> None:
    conn.execute(
        """
        INSERT INTO parameters(
            compound_id,
            molecular_weight, logp, hbd, hba,
            polar_surface_area, rotatable_bonds, aromatic_rings, heavy_atom_count,
            qed,
            ro5_violation, ro5_pass, ro3_pass, leadlike_pass, veber_pass, bioavailability_pass
        ) VALUES (
            :compound_id,
            :molecular_weight, :logp, :hbd, :hba,
            :polar_surface_area, :rotatable_bonds, :aromatic_rings, :heavy_atom_count,
            :qed,
            :ro5_violation, :ro5_pass, :ro3_pass, :leadlike_pass, :veber_pass, :bioavailability_pass
        );
        """,
        {
            "compound_id": compound_id,
            "molecular_weight": props["mw"],
            "logp": props["logp"],
            "hbd": props["hbd"],
            "hba": props["hba"],
            "polar_surface_area": props["tpsa"],
            "rotatable_bonds": props["rotb"],
            "aromatic_rings": props["rings"],
            "heavy_atom_count": props["hatom"],
            "qed": props["qed"],
            "ro5_violation": props["ro5_violations"],
            "ro5_pass": props["ro5_pass"],
            "ro3_pass": props["ro3_pass"],
            "leadlike_pass": props["leadlike_pass"],
            "veber_pass": props["veber_pass"],
            "bioavailability_pass": props["bioavail_pass"],
        },
    )


