"""
Functions for getting the molprobity values of structures and comparing to distributions
"""
import re
from pathlib import Path
from typing import Union
import shutil
import subprocess

import numpy.typing as npt

from .types import TStructure, TAtomID, TGrid, TReflections, TMolprobity

from .constants import MolprobityConstants as MC
from . import constants

PANDDAS_PATH = '/dls/data2temp01/labxchem/data/2017/lb18145-17/processing/edanalyzer/output/pandda_new_score/panddas_new_score/'
TEST_MOLPROBITY_FILE = '/dls/data2temp01/labxchem/data/2017/lb18145-17/processing/edanalyzer/output/pandda_new_score/panddas_new_score/KAT2B/processed_datasets/KAT2B-x0150/molprobity.out'


class MolprobityStats(TMolprobity):
    def __init__(self, stats):
        for k, v in stats.items():
            self.__setattr__(k, v)


def parse_molprobity_file(path) -> TMolprobity:    
    # Get text 
    with open(path, 'r') as f:
        text = f.read()
    
    return MolprobityStats(
        {
            'ramachandran_outliers': float(re.findall(MC.RAMACHANDRAN_OUTLIERS, text)[0]),
            'ramachandran_favored': float(re.findall(MC.RAMACHANDRAN_FAVORED, text)[0]),
            'rotamer_outliers': float(re.findall(MC.ROTAMER_OUTLIERS, text)[0]),
            'beta_deviations': float(re.findall(MC.BETA_DEVIATIONS, text)[0]),
            'clashscore': float(re.findall(MC.CLASHSCORE, text)[0]),
            'rms_bonds': float(re.findall(MC.RMS_BONDS, text)[0]),
            'rms_angles': float(re.findall(MC.RMS_ANGLES, text)[0]),
            'molprobity': float(re.findall(MC.MOLPROBITY, text)[0]),
            'resolution': float(re.findall(MC.RESOLUTION, text)[0]),
            'rwork': float(re.findall(MC.RWORK, text)[1]),
            'rfree': float(re.findall(MC.RFREE, text)[1]),
            'wavelength': float(re.findall(MC.WAVELENGTH, text)[0]),
            'resolution_full_low': list(map(float, re.findall(MC.RESOLUTION_FULL, text)[0]))[0],
            'resolution_full_high': list(map(float, re.findall(MC.RESOLUTION_FULL, text)[0]))[1],
            'num_reflections': list(map(int, re.findall(MC.NUM_REFLECTIONS, text)[0]))[0],
            'num_reflections_unique': list(map(int, re.findall(MC.NUM_REFLECTIONS, text)[0]))[1],
            'num_reflections_non_outlying': list(map(int, re.findall(MC.NUM_REFLECTIONS_NON_OUTLYING, text)[0]))[0],
            'num_reflections_non_outlying_unique': list(map(int, re.findall(MC.NUM_REFLECTIONS_NON_OUTLYING, text)[0]))[1],
            'completeness': list(map(float, re.findall(MC.COMPLETENESS, text)[0]))[0],
            'completeness_unique': list(map(float, re.findall(MC.COMPLETENESS, text)[0]))[1],
    }
    )

def run_molprobity(st_path: Path, mtz_path: Path) -> TMolprobity:
    if shutil.which(constants.MOLPROBITY_COMMAND):
        p=subprocess.Popen(
        f'{constants.MOLPROBITY_COMMAND} {st_path} {mtz_path}',
        shell=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE
    )
        stdout, stderr = p.communicate()
        print(str(stdout))
        print(str(stderr))
    else:
        raise Exception(f'{constants.MOLPROBITY_COMMAND} not in path!')

def molprobity(
        structure: Union[TStructure, Path], 
        reflections: Union[TReflections, Path]
        ) -> TMolprobity:
    """
    Returns the Molprobity stats of the structure (probably using CCP4).
    Note -stats- part - include where sits on PDB distributions!
    """
    match (structure, reflections):
        case (TStructure(), TReflections()):
            mp = run_molprobity(structure.path, reflections.path)
        case (TStructure(), Path()):
            mp = run_molprobity(structure.path, reflections)
        case (Path(), TStructure()):
            mp = run_molprobity(structure, reflections.path)
        case (Path(), Path()):
            mp = run_molprobity(structure, reflections)
        case _:
            raise Exception(f"Structure must be a path or have the path attribute!")

    return mp
    ...