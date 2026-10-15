"""
Functions for calculating the RSCC between structures
"""
from pathlib import Path
import tempfile
import subprocess
import shutil
import re

import numpy.typing as npt

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
    pattern = r'Map-map CC is ([\S]+)'
    match = re.search(pattern, stdout)
    print(match)
    return float(match.groups(0))

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
        print(e)
        rscc = None
    
    try:
        rsccs = parse_rscc_file(Path(out_dir) / '_cc_per_residue.log')
        print(rscc)
    except:
        rsccs = None
        # raise Exception(
        #     f'Molprobity command: "{script}\n'
        #     f'Molprobity stdout: {stdout}\n'
        #     f'Molprobity stderr: {stderr}'
        #                 )

    return rscc, rsccs

def get_rscc(
        st: TStructure | Path| None = None,
        mtz: TReflections| Path| None= None,
        xmap: TGrid| Path| None= None,
        st_2: TStructure| Path| None= None,
        mtz_2:TReflections| Path| None= None,
        xmap_2: TGrid| Path | None= None
):
    st = _get_path(st)
    mtz = _get_path(mtz)
    xmap = _get_path(xmap)
    st_2 = _get_path(st_2)
    mtz_2 = _get_path(mtz_2)
    xmap_2 = _get_path(xmap_2)

    with tempfile.TemporaryDirectory(delete=False) as fp:

        print([st, mtz, xmap, st_2, mtz_2, xmap_2])
        match (st, mtz, xmap, st_2, mtz_2, xmap_2):
            case (Path(), Path(), None, None, None, None):
                script = PHENIX_MAP_CORRELATIONS.format(
                    dtype1=PHENX_MODEL,
                    path1=st,
                    dtype2=PHENX_MAP_COEFFS_1,
                    path2=mtz,
                    out_dir=fp
                )
                print(script)
                rscc, rsccs = run_phenix_rscc_script(script, fp)
                print(rscc)
                print(rsccs)
            case (Path(), None, Path(), None, None, None):
                script = PHENIX_MAP_CORRELATIONS.format(
                                    dtype1=PHENX_MODEL,
                                    path1=st,
                                    dtype2=PHENX_XMAP_1,
                                    path2=xmap,
                                    out_dir=fp
                                )
                print(script)
                run_phenix_rscc_script(script, fp)
            case (Path(), None, None, Path(), None, None,):
                ...
            case (None, Path(), None, None, Path(), None,):
                        ...
            case (None, Path(), None, None, None, Path(),):
                ...
            case (None, None, Path(), None, None, Path(),):
                script = PHENIX_MAP_CORRELATIONS.format(
                    dtype1=PHENX_XMAP_1,
                    path1=xmap,
                    dtype2=PHENX_XMAP_2,
                    path2=xmap_2,
                    out_dir=fp
                )
                print(script)
                run_phenix_rscc_script(script, fp)
            case (_, _, _, _, _, _):
                raise Exception(f'Wrong input types!')