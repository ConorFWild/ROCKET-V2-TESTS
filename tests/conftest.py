import pytest

import gemmi

from rocket_testing import constants 
from rocket_testing.types import TStructure, TAtomID, TGrid, TReflections, Structure, Reflections, TMolprobity, Grid
from rocket_testing.selections import get_atom_ids
from rocket_testing.molprobity import MolprobityStats

@pytest.fixture
def structure() -> TStructure:
    return Structure(constants.TEST_STRUCTURE_PATH)

@pytest.fixture
def offset_structure(structure) -> TStructure:
    offset_structure = structure.clone()
    for chain in offset_structure[0]:
            for res in chain:
                for atom in res:
                    atom.pos = atom.pos + gemmi.Position(0.1,0.1,0.1)
    return offset_structure

@pytest.fixture
def reflections() -> TReflections:
    return Reflections(constants.TEST_REFLECTIONS_PATH)

@pytest.fixture
def grid() -> TGrid:
    return Grid(constants.TEST_GRID_PATH)

@pytest.fixture
def all_atom_selection(structure) -> list[TAtomID]:
    return list(get_atom_ids(structure))


@pytest.fixture
def chain_selection(structure)  -> list[TAtomID]:
    return [_atom_id for _atom_id in get_atom_ids(structure) if _atom_id.chain == 'A']

@pytest.fixture
def residue_selection(structure) -> list[TAtomID]:
    return [_atom_id for _atom_id in get_atom_ids(structure) if (_atom_id.chain == 'A') & (_atom_id.insertion =='63')]



@pytest.fixture
def residues_selection(structure) -> list[TAtomID]:
    residue_insertions = ["59", "60", "61"]
    return [_atom_id for _atom_id in get_atom_ids(structure) if (_atom_id.chain == 'A') & (_atom_id.insertion in residue_insertions)]

@pytest.fixture
def residues_selection_swapped(structure) -> list[TAtomID]:
    residue_insertions = ["62",]
    return [_atom_id for _atom_id in get_atom_ids(structure) if (_atom_id.chain == 'A') & (_atom_id.insertion in residue_insertions)]


@pytest.fixture
def atom_selection(structure) -> list[TAtomID]:
    return [
        _atom_id 
        for _atom_id 
        in get_atom_ids(structure) 
        if (_atom_id.chain == 'A') & (_atom_id.insertion == "63") & (_atom_id.name == 'CA')]

@pytest.fixture
def incompatible_selection(structure) -> list[TAtomID]:
    return [
        _atom_id 
        for _atom_id 
        in get_atom_ids(structure) 
        if (_atom_id.chain == 'A') & (_atom_id.insertion == "60") & (_atom_id.name == 'CA')]



@pytest.fixture
def perturbed_structure() -> TStructure:
    return ...

@pytest.fixture
def structure_swapped(structure):
    _structure_swapped = structure.clone()
    pos_1 = structure[0]['A']['62'][0]["OE1"][0].pos
    pos_2 = structure[0]['A']['62'][0]["OE2"][0].pos
    _structure_swapped[0]['A']['62'][0]["OE1"][0].pos = pos_2
    _structure_swapped[0]['A']['62'][0]["OE2"][0].pos = pos_1
    return _structure_swapped


@pytest.fixture
def cut_structure(structure) -> TStructure:
    cut_structure = structure.clone()
    chains_to_delete = []
    for model in cut_structure:
        for chain in model:
            chains_to_delete.append(chain.name)

    for chain_name in chains_to_delete:
        del cut_structure[0][chain_name]

    for model in structure:
        for chain in model:
            new_chain = gemmi.Chain(chain.name)
            for res in chain:
                if (chain.name == 'A') and (res.seqid.num == 60):
                    continue
                new_chain.add_residue(res)
            cut_structure[0].add_chain(new_chain) 
    return cut_structure

@pytest.fixture
def known_rfree() -> TStructure:
    return ...

@pytest.fixture
def known_molprobity() -> TMolprobity:
    return MolprobityStats(
        {'ramachandran_outliers': 0.0, 
         'num_reflections_non_outlying_unique': 1464, 
         'resolution_full_high': 1.91, 
         'num_reflections_unique': 1464, 
         'rms_bonds': 0.014, 
         'completeness_unique': 92.08, 
         'num_reflections_non_outlying': 16354, 
         'molprobity': 1.76, 
         'resolution_full_low': 41.703, 
         'resolution': 1.91, 
         'clashscore': 5.62, 
         'wavelength': 0.9221, 
         'rfree': 0.2839, 
         'rms_angles': 1.85, 
         'beta_deviations': 0.0, 
         'rwork': 0.2483, 
         'rotamer_outliers': 1.88, 
         'ramachandran_favored': 96.34, 
         'completeness': 99.16, 
         'num_reflections': 16357,
         }
         )