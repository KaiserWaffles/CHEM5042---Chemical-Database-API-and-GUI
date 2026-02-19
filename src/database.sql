--Two tables to enable future proofing of existing compounds with multiple sets of parameters/properties
--LOOK UP WHAT ATTRIBUTES ARE NEEDED
CREATE TABLE compound (
  compound_id INTEGER AUTOINCREMENT,
  name VARCHAR (255), NOT NULL
  smiles VARCHAR (255) NOT NULL UNIQUE,
  formula VARCHAR (255),
  PRIMARY KEY (compound_id)
);

--READ PAPERS ON WHAT GOES IN PARAMETERS
CREATE TABLE parameters (
  compound_id INTEGER, 
  molecular_weight REAL NOT NULL,
  logp REAL NOT NULL,
  hbd INTEGER NOT NULL,
  hba INTEGER NOT NULL,
  polar_surface_area REAL NOT NULL,
  rotatable_bonds INTEGER NOT NULL,
  aromatic_rings INTEGER NOT NULL,
  heavy_atom_count INTEGER NOT NULL,

  ro5_violation INTEGER NOT NULL, --NOT BOOLEAN; counts if how many conditions are violated
  ro5_pass INTEGER NOT NULL,
  ro3_pass INTEGER NOT NULL,
  leadlike_pass INTEGER NOT NULL,
  veber_pass INTEGER NOT NULL,
  bioavailability_pass INTEGER NOT NULL,
  PRIMARY KEY (compound_id)
  FOREIGN KEY (compound_id) REFERENCES compound(compound_id) ON DELETE CASCADE
);