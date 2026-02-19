from dataclasses import dataclass
from typing import Iterator, Optional
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

# Structure for RDkit
@dataclass(frozen=True)
class ParsedMol:
    """
    Parsed compound record from SDF; 
        name:    Compound name (from SDF _Name or auto-generated).
        mol:     RDKit Mol object for property calculation.
        smiles:  Canonical SMILES string.
        formula: Molecular formula (e.g. 'C9H8O4'), or None on failure.
    """
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
    removeHs=True so that explicit hydrogens don't interfere with descriptor calculations and SMILES generation.
    Args:
        sdf_path: Path to the .sdf file (e.g. 'Molecule2.sdf').
    Yields:
        ParsedMol: One per successfully parsed molecule.
    """
    sdf_path = Path(sdf_path)
    if not sdf_path.exists():
        raise FileNotFoundError(f"SDF file not found: {sdf_path}")
    
    # sanitize catches invalid structures
    supplier = Chem.SDMolSupplier(str(sdf_path), removeHs = False, sanitize = True)

    skipped = 0
    parsed = 0

    for i, mol in enumerate(supplier, start=1):
        if mol is None:
            print(f"  WARNING: Molecule at index {i} could not be parsed — skipping.")
            skipped += 1
            continue

        try:
            name = extract_sdf_name(mol, fallback=f"comp{i}")

            # Generate canonical SMILES (isomeric preserves stereochemistry)
            smiles = Chem.MolToSmiles(mol, isomericSmiles=True)

            # Calculate molecular formula
            formula = None
            try:
                formula = rdMolDescriptors.CalcMolFormula(mol)
            except Exception:
                formula = None

            parsed += 1
            yield ParsedMol(name=name, mol=mol, smiles=smiles, formula=formula)

        except Exception as e:
            print(f"  WARNING: Error processing molecule {i}: {e} — skipping.")
            skipped += 1
            continue
        
    print(f"Parsed {parsed} compounds from: {sdf_path.name}")
    if skipped > 0:
        print(f"Skipped {skipped} molecules due to parsing errors.")