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
    """
    Insert compound identity row into the compound database table
    """
    cur = conn.execute(
        "INSERT INTO compound(name, smiles, formula, source_file) VALUES (?, ?, ?);",
        (name, smiles, formula),
    )
    return int(cur.lastrowid)


def insert_properties(conn: sqlite3.Connection, compound_id: int, props: Dict[str, object]) -> None:
    """
    Insert calculated parameters and filter results for a compound.

        compound_id: FK linking to compound table.
        props:       Dictionary from parameters.db_row() containing
                     all descriptor values and filter flags.
    """
    conn.execute(
        """
        INSERT INTO parameters(
            compound_id,
            molecular_weight, logp, hbd, hba,
            polar_surface_area, rotatable_bonds, aromatic_rings,
            heavy_atom_count, qed,
            ro5_violation, ro5_pass, ro3_pass, leadlike_pass,
            veber_pass, bioavailability_pass
        ) VALUES (
            :compound_id,
            :mw, :logp, :hbd, :hba,
            :tpsa, :rotb, :rings,
            :hatom, :qed,
            :ro5_violations, :ro5_pass, :ro3_pass, :leadlike_pass,
            :veber_pass, :bioavail_pass
        );
        """,
        {"compound_id": compound_id, **props},
    )

# QUERIES

def get_compounds(conn: sqlite3.Connection) -> List[sqlite3.Row]:
    """
    Select all compounds with associated paramets
    Inner join so only compounds with calculated parameters exists
    """
    cur = conn.execute("""
        SELECT c.compound_id, c.name, c.smiles, c.formula,
               p.molecular_weight, p.logp, p.hbd, p.hba,
               p.polar_surface_area, p.rotatable_bonds,
               p.aromatic_rings, p.heavy_atom_count, p.qed,
               p.ro5_violation, p.ro5_pass, p.ro3_pass,
               p.leadlike_pass, p.veber_pass, p.bioavailability_pass
        FROM compound c
        INNER JOIN parameters p ON c.compound_id = p.compound_id
        ORDER BY c.compound_id;
    """)
    return cur.fetchall()

def get_compound_counts(conn: sqlite3.Connection) -> int:
    """
    Get total nu ber of compounds in the database
    library of 100 or 1,000,000+ compounds possibleg
    """
    cur = conn.execute("SELECT COUNT(*) FROM compound;")
    return cur.fetchone()