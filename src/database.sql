--Two tables to enable future proofing of existing compounds with multiple sets of parameters/properties
--LOOK UP WHAT ATTRIBUTES ARE NEEDED
CREATE TABLE compound (
  compound_id INTEGER AUTOINCREMENT,
  name VARCHAR (255),
  smiles VARCHAR (255) NOT NULL UNIQUE,
  formula VARCHAR (255),
  PRIMARY KEY (compound_id)
);

--READ PAPERS ON WHAT GOES IN PARAMETERS
CREATE TABLE parameters (
  compound_id INTEGER,
  bioavailability VARCHAR,
  PRIMARY KEY (compound_id)
  FOREIGN KEY (compound_id) REFERENCES compound(compound_id) ON DELETE CASCADE
);