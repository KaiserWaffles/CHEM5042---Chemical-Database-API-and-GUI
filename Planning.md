Steps:
"""
1 - Create ERD for relational database schema
2 - Create repo skeleton for modules
    > main.py =         process program and design
    > sdf_parser.py =   reads & parse .sdf file (RDKit)
    > database.py =     SQLite schema
    > parameters.py =   derived parameters calculation and filters
    > GUI.py =          Tkinter
3 - Plan flowchart 
4 - Write Report alongside steps 1 and 2
    > Methods
        > Resources used
        > approach for creating the database
    > Results
        > How database can be searched
        > How database is used (filtering)
        > Screenshots of GUI and results for the set of 100 compounds
"""


TASK 1: Implement a database of chemical compounds

ENSURE:
    1. Create a design for the program in outline
    2. Document the code
    3. Define filtering criteria
    4. Consider what challenges you would face in populating database with 1,000 or 1,000,000 compounds and how they could be overcome


A database for the storage and retrieval of a series of chemical compounds WITH derived parameters relevant to their use in drug discovery
    1. Read and parse the chemical structures from a single .sdf file [Python]
    2. Implement a relational database  to populate with molecular description and associated parameters calculated for each molecule [SQL]
    3. Populate database with the library provided (100 compounds) [SQL or Python?]
    4. Implement functionality to calculate the derived parameters you wish to store and populate specific fields within your database [ RDKit ]
    5. Produce a GUI interface to allow viewing of the compounds in your database together with appropriate parameters calculated for each molecule [Python]
    6. Implement functionality to rank the entries in your database to allow filtering of the compounds by parameters such as LogP or Molecular Weight [Python]
    7. Filter compounds using the widely accepted criteria such as Lipinkski’s rules, lead likeness and/or bioavailability [Python]


Implement programming in Python, building upon BIOL4292.

Create a suitable SQLite database to contain chemical data and associated fields

Produce GUI front end to the database, allowing visualisation of relevant data such as compound name, formula, molecular weight etc. and sort on the different fields

A library of compounds (comp1, comp2 …etc) is provided in a multiple compound format  on Moodle. Files are in .sdf format as described in lectures and are not typical drug like molecules but rather bioactive compounds.  This may cause issues with ‘Molecular Weight’ or other factors for certain compounds that may need to be excluded from the compounds stored within the database

RDKIT: Open-Source Cheminformatics Software (rdkit.org) will help to read and calculate properties for the compounds in the .sdf file

TKinter is recommended for developing the GUI (pythonguides.com/python-tkinter-table-tutorial)


Personal Questions to look up 
    • What makes a “typical drug like molecule”
    • What are “bioactive compounds”
    • What is “lead likeness” as a form of filtering compounds



TASK 2: Generate flow diagram to represent how you would generate multiple conformers of a given molecule to sample conformational space


Consider how from coordinates of a chemical in an .sdf format, one could generate multiple conformers of the molecule. Only consider rotation of parts of the molecule relative to one another around “rotatable” bonds. Pay attention to how one could restrict the number of conformations to a minimum while still sampling conformational space effectively

ENSURE:
    1. Flow diagram produced is easy to follow and should cover the main steps that you would consider when writing a program
    2. Pay particular attention to how one could restrict the number of conformations
    3. Consider how one might store the coordinates once produced

Jmol (jmol.sourceforge.net) will visualise the compound and label all atoms. Include an image of this in the report as it will help describe the process of generating multiple conformations. The flow chart will need to work for your molecule, (if time, consider how to make it more generic)


Personal Questions to answer
    • Define and differentiate a chemical and a molecule
    • Define what a conformer is 


                Task 3: Report Writeup

Brief introduction to “In-Silicon Valley Drug Discovery” to give sense of why the database of chemicals is being produced (TASK 1).
Methods section to outline how you approached producing your searchable database of compounds, what resources used etc.
Results section to explain how your database can be searched, its look and feel (include images of how the GUI looks and the results you have for the set of 100 compounds)

Task 2 is standalone section; just a methods section to outline the strategy for generating multiple conformations of a chemical

Finish report with
Conclusions and references found useful. 

Typically 5-10 pages (slightly longer if lots of figures).


FINAL TASK: Discuss the database and talk through what I have done through this assignment and planning process with Dr Lapthorn.
