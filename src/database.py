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
