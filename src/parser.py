import io
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors, Draw, QED
from rdkit.Chem import rdDepictor
from dataclasses import dataclass
from typing import Iterator, Optional

# Structure for RDkit
@dataclass(frozen=True)
class ParsedMol:
    """Parsed compound record from SDF; used for description."""
    name: str
    mol: Chem.Mol
    smiles: str
    formula: Optional[str]


# API
def parse_sdf(filepath: str) -> list[dict]:
    """
    Parse .SDF file and return a list of compound parameters dicts; Each valid compound is tidied by RDKit, properties calculated, and an optional molecular weight (MW) ceiling filter applied.  Molecules unable to be parsed or exceed MW_MAX_THRESHOLD skipped

    Returns: list[dict]
        One dict per compound, with KEYS matching the SQLite schema
    """

