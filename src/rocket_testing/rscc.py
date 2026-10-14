"""
Functions for calculating the RSCC between structures
"""
from pathlib import Path
import tempfile

import numpy.typing as npt

from .types import TStructure, TAtomID, TGrid, TReflections

PHENX_MODEL = "model"
PHENX_MAP_COEFFS_1 = "map_coeffs_1"
PHENX_MAP_COEFFS_2 = "map_coeffs_1"
PHENX_XMAP_1 = "map_in_1"
PHENX_XMAP_1 = "map_in_2"

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

            case (Path(), None, Path(), None, None, None):
                ""
                ...
            case (Path(), None, None, Path(), None, None,):
                ...
            case (None, Path(), None, None, Path(), None,):
                        ...
            case (None, Path(), None, None, None, Path(),):
                ...
            case (None, None, Path(), None, None, Path(),):
                ...
            case (_, _, _, _, _, _):
                raise Exception(f'Wrong input types!')