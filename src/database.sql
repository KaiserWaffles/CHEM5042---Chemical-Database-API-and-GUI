--Two tables to enable future proofing of existing compounds with multiple sets of parameters/properties
--LOOK UP WHAT ATTRIBUTES ARE NEEDED
CREATE TABLE IF NOT EXISTS compound (
  compound_id INTEGER,
  name VARCHAR (255) NOT NULL,
  smiles VARCHAR (255) NOT NULL UNIQUE,
  formula VARCHAR (255),
  PRIMARY KEY (compound_id)
);

--READ PAPERS ON WHAT GOES IN PARAMETERS
CREATE TABLE IF NOT EXISTS parameters (
    compound_id         INTEGER PRIMARY KEY,
    -- Core molecular descriptors (from RDKit)
    molecular_weight  REAL  NOT NULL,
    logp  REAL  NOT NULL,
    hbd INTEGER NOT NULL,
    hba INTEGER NOT NULL,
    polar_surface_area  REAL  NOT NULL,
    rotatable_bonds INTEGER NOT NULL, 
    aromatic_rings  INTEGER NOT NULL,
    heavy_atom_count  INTEGER NOT NULL,
    qed REAL  NOT NULL,

    -- Pre-computed filter results (stored for fast querying)
    ro5_violation INTEGER NOT NULL, -- Count of Ro5 violations (0-4)
    ro5_pass  INTEGER NOT NULL, -- 1 if violations <= 1, else 0
    ro3_pass  INTEGER NOT NULL, -- 1 if passes Rule of 3 (fragment)
    leadlike_pass INTEGER NOT NULL, -- 1 if passes lead-likeness
    veber_pass  INTEGER NOT NULL, -- 1 if passes Veber rules
    bioavailability_pass  INTEGER NOT NULL, -- 1 if passes bioavailability
    FOREIGN KEY (compound_id) REFERENCES compound(compound_id) ON DELETE CASCADE
);

-- Indexes for query performance
-- Essential when scaling to 1,000+ compounds
CREATE INDEX IF NOT EXISTS idx_params_mw
    ON parameters(molecular_weight);

CREATE INDEX IF NOT EXISTS idx_params_logp
    ON parameters(logp);

CREATE INDEX IF NOT EXISTS idx_params_ro5
    ON parameters(ro5_pass);

CREATE INDEX IF NOT EXISTS idx_params_veber
    ON parameters(veber_pass);
