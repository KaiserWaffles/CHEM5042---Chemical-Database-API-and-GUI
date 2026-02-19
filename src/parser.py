from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
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


def parse_sdf(sdf_path: str | Path) -> Iterator[ParsedMol]:
    """
    Read molecules from an SDF file; yields ParsedMol objects;
    Invalid records (mol is None) are skipped.
    """
    sdf_path = Path(sdf_path)
    supplier = Chem.SDMolSupplier(str(sdf_path), removeHs=False)

    for i, mol in enumerate(supplier, start=1):
        if mol is None:
            continue

        name = extract_sdf_name(mol, fallback=f"comp{i}")

        # SMILES is derived from structure in SDF
        smiles = Chem.MolToSmiles(mol, isomericSmiles=True)

        formula = None
        try:
            formula = rdMolDescriptors.CalcMolFormula(mol)
        except Exception:
            formula = None

        yield ParsedMol(name=name, mol=mol, smiles=smiles, formula=formula)
