"""
Functions for calculating the RSCC between structures
"""
from pathlib import Path
import tempfile
import subprocess
import shutil
import re

import numpy.typing as npt
import gemmi

from .types import TStructure, TAtomID, TGrid, TReflections
from . import constants

PHENX_MODEL = "model"
PHENX_MAP_COEFFS_1 = "map_coeffs_1"
PHENX_MAP_COEFFS_2 = "map_coeffs_2"
PHENX_XMAP_1 = "map_in_1"
PHENX_XMAP_2 = "map_in_2"
PHENIX_MAP_CORR_RES_FILE = '_cc_per_residue.log'

PHENIX_MAP_CORRELATIONS = "phenix.map_correlations {dtype1}={path1} {dtype2}={path2} file_name_prefix={out_dir}/" 

def rscc_st_to_st(
        ref_structure: TStructure, 
        mov_structure: TStructure, 
        atom_selection: list[TAtomID],
        ) -> float:
    """
    Returns the RSCC between the predicted maps of two structures
    """
    ...


def rscc_st_to_grid(
    st: TStructure,
    xmap: TGrid,
    atom_selection: list[TAtomID],
) -> float:
    """
    Returns the RSCC between the predicted maps of a structure and a map
    """

def _get_path(st: TStructure | Path| None) -> Path | None:
    match st:
        case TStructure():
            return st.path
        case TReflections():
            return st.path
        case TGrid():
            return st.path
        case Path():
            return st
        case None:
            return None
        case _:
            raise Exception(f'Should have structure, path or None, got: {st}')

def parse_rscc_output(stdout):
    pattern = r'Map-map CC is ([\.0-9]+)'
    match = re.search(pattern, stdout)
    return float(match.group(1))

def parse_rscc_file(path):
    pattern = r"([\S]+)[\s]+([\S]+)[\s]+([\S]+)[\s]+([\S]+)[\s]+([\S]+)[\s]+([\S]+)[\s]+"
    with open(path , 'r') as f:  
        matches = re.findall(pattern, f.read(),)

    data = []
    for match in matches:
        data.append( {
            'chain':match[0],
            'res_name':match[1],
            'insertion': match[2],
            'rscc': float(match[3]),
            'b': float(match[4]),
            'occ': float(match[5]),
        })
    return data

def run_phenix_rscc_script(script, out_dir):
    if shutil.which(constants.PHENIX_MAP_CORR):
        p=subprocess.Popen(
            script,
            shell=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE
        )
        stdout, stderr = p.communicate()
    else:

        raise Exception(
            f'{constants.PHENIX_MAP_CORR} not in path!\n'
                        )
    
    try:
        rscc = parse_rscc_output(str(stdout))
    except Exception as e:
        rscc = None
    
    try:
        rsccs = parse_rscc_file(Path(out_dir) / '_cc_per_residue.log')
    except:
        rsccs = None


    match (rscc, rsccs):
        case (None, None):
            raise Exception(
                f'Did not get any RSCCs!\n'
                f'Phenix command: {script}\n'
                f'Phenix map.correlation stdout: {stdout}\n'
                f'Phenix map.correlation stderr: {stderr}\n'
            )

    return rscc, rsccs


def _make_structure_mtz_file(st_path, tmp_dir):
    st = gemmi.read_structure(str(st_path))
    out_path = Path(tmp_dir) / 'sf.mtz'
    script = f"phenix.fmodel {st_path} high_resolution={st.resolution} file_name={out_path}"
    if shutil.which('phenix.fmodel'):
        p=subprocess.Popen(
            script,
            shell=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE
        )
        stdout, stderr = p.communicate()
    else:
        raise Exception()

    if out_path.exists():
        return out_path
    else:
        raise Exception()

def get_rscc(
        d1: TStructure | TReflections | TGrid | Path,
        d2: TStructure | TReflections | TGrid | Path 
):

    dtypes = []
    for j, extension in enumerate([_get_path(d1).suffix, _get_path(d2).suffix]):
        match extension:
            case '.pdb' | '.cif' | '.mmcif':
                dtypes.append(PHENX_MODEL)
            case '.ccp4':
                if PHENX_XMAP_1 in dtypes:
                    dtypes.append(PHENX_XMAP_2)
                else:
                    dtypes.append(PHENX_XMAP_1)
            case '.mtz':
                if PHENX_MAP_COEFFS_1 in dtypes:
                    dtypes.append(PHENX_MAP_COEFFS_2)
                else:
                    dtypes.append(PHENX_MAP_COEFFS_1)
            case _:
                if j == 0:
                    raise Exception(f'{d1} is not a pdb, cif, mmcif, ccp4 or mtz file!')
                else:
                    raise Exception(f'{d2} is not a pdb, cif, mmcif, ccp4 or mtz file!')
                                    

    with tempfile.TemporaryDirectory(delete=False) as fp:
        match dtypes:
            case [x, y] if (x == y == PHENX_MODEL):
                mtz_file = _make_structure_mtz_file(_get_path(d2), fp)
                script = PHENIX_MAP_CORRELATIONS.format(
                                                    dtype1=PHENX_MODEL,
                                                    path1=_get_path(d1),
                                                    dtype2=PHENX_MAP_COEFFS_1,
                                                    path2=mtz_file,
                                                    out_dir=fp
                                                )
                rscc, rsccs = run_phenix_rscc_script(script, fp)
            case _:
                script = PHENIX_MAP_CORRELATIONS.format(
                                    dtype1=dtypes[0],
                                    path1=_get_path(d1),
                                    dtype2=dtypes[1],
                                    path2=_get_path(d2),
                                    out_dir=fp
                                )
                rscc, rsccs = run_phenix_rscc_script(script, fp)


    return rscc, rsccs