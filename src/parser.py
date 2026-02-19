import io
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors, Draw, QED
from rdkit.Chem import rdDepictor
from dataclasses import dataclass
from typing import Iterator, Optional
from pathlib import Path

# Structure for RDkit
@dataclass(frozen=True)
class ParsedMol:
    """Parsed compound record from SDF; used for description."""
    name: str
    mol: Chem.Mol
    smiles: str
    formula: Optional[str]


# API
def extract_sdf_name(mol: Chem.Mol, fallback: str) -> str:
    """Extract name from common SDF fields; fallback to compound name."""
    for key in ("_Name", "Name", "ID", "CompoundID", "MOL_NAME"):
        if mol.HasProp(key):
            value = mol.GetProp(key).strip()
            if value:
                return value
    return fallback



