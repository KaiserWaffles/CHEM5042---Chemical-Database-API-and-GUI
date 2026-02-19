from dataclasses import dataclass
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors, QED

# structure for RDKit
@dataclass(frozen=True)
class Parameters:
    mw: float
    logp: float
    hbd: int
    hba: int
    # rotatble bonds
    rotb: int 
    # aromatic rings
    rings: int
    # polar surface area
    tpsa: float
    qed: float
    ro5_violations: int
    ro5_pass: int
    ro3_pass: int
    leadlike_pass: int
    veber_pass: int
    bioavail_pass: int
