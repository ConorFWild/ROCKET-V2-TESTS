from pathlib import Path

TEST_STRUCTURE_PATH = Path('./test_data/XX01ZVNS2B-x3819/dimple.pdb')
TEST_STRUCTURE_WITH_LIGAND_PATH = Path('./test_data/XX01ZVNS2B-x3819/refine.pdb')
TEST_REFLECTIONS_PATH = Path('./test_data/XX01ZVNS2B-x3819/dimple.mtz')
TEST_GRID_PATH = Path('./test_data/XX01ZVNS2B-x3819/dimple_2mFo-DFc_map.ccp4')
TEST_GRID_PATH_2 = Path('./test_data/XX01ZVNS2B-x3819/dimple_2mFo-DFc_map_2.ccp4')
TEST_CIF_PATH = Path('./test_data/XX01ZVNS2B-x3819/ASAP-0031519-001.cif')
MONOMER_JSON = Path('./test_data/monomers.json')
RESIDUES = ["ALA",
                 "ARG",
                 "ASN",
                 "ASP",
                 "CYS",
                 "GLN",
                 "GLU",
                 "HIS",
                 "ILE",
                 "LEU",
                 "LYS",
                 "MET",
                 "PHE",
                 "PRO",
                 "SER",
                 "THR",
                 "TRP",
                 "TYR",
                 "VAL",
                 "GLY",
                 ]

COMMANDS = {
    ''
}

MOLPROBITY_COMMAND = 'molprobity.molprobity'

class MolprobityConstants:
    RAMACHANDRAN_OUTLIERS = r'Ramachandran outliers =\s+([0-9.]+)'
    RAMACHANDRAN_FAVORED = r'\s+favored =\s+([0-9.]+)'
    ROTAMER_OUTLIERS = r'\s+Rotamer outliers\s+=\s+([0-9.]+)'
    BETA_DEVIATIONS = r'\s+C-beta deviations\s+=\s+([0-9.]+)'
    CLASHSCORE = r'\s+Clashscore\s+=\s+([0-9.]+)'
    RMS_BONDS = r'\s+RMS\(bonds\)\s+=\s+([0-9.]+)'
    RMS_ANGLES = r'\s+RMS\(angles\)\s+=\s+([0-9.]+)'
    MOLPROBITY = r'\s+MolProbity score\s+=\s+([0-9.]+)'
    RESOLUTION = r'\s+Resolution\s+=\s+([0-9.]+)'
    RWORK = r'\s+R-work\s+=\s+([0-9.]+)'
    RFREE = r'\s+R-free\s+=\s+([0-9.]+)'
    WAVELENGTH = r'\s+Wavelength\s+=\s+([0-9.]+)'
    RESOLUTION_FULL=  r'\s+Resolution range\s+=\s+([0-9.]+)\s+-\s+([0-9.]+)\s+\(([0-9.]+)\s+-\s+([0-9.]+)\)'
    NUM_REFLECTIONS= r'\s+Number of reflections\s+=\s+([0-9.]+)\s+\(([0-9.]+)\)'
    NUM_REFLECTIONS_NON_OUTLYING = r'\s+after outlier rejection\s+=\s+([0-9.]+)\s+\(([0-9.]+)\)'
    COMPLETENESS = r'\s+Completeness\s+=\s+([0-9.]+)%\s+\(([0-9.]+)%\)'
    MOLPROBITY_CSV_PATH = '/dls/data2temp01/cfw/pandda_paper/data/molprobity.csv'

PHENIX_MAP_CORR = "phenix.map_correlations"
