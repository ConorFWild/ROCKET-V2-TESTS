import pytest

import gemmi

from rocket_testing import constants 
from rocket_testing.types import TStructure, TAtomID, TGrid, TReflections, Structure
from rocket_testing.selections import get_atom_ids

@pytest.fixture
def structure() -> TStructure:
    return Structure(constants.TEST_STRUCTURE_PATH)

@pytest.fixture
def offset_structure(structure) -> TStructure:
    offset_structure = structure.clone()
    for chain in offset_structure.st[0]:
            for res in chain:
                for atom in res:
                    atom.pos = atom.pos + gemmi.Position(0.1,0.1,0.1)
    return offset_structure

@pytest.fixture
def reflections() -> TReflections:
    return ...


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
    pos_1 = structure.st[0]['A']['62'][0]["OE1"][0].pos
    pos_2 = structure.st[0]['A']['62'][0]["OE2"][0].pos
    _structure_swapped.st[0]['A']['62'][0]["OE1"][0].pos = pos_2
    _structure_swapped.st[0]['A']['62'][0]["OE2"][0].pos = pos_1
    return _structure_swapped


@pytest.fixture
def cut_structure(structure) -> TStructure:
    cut_structure = structure.clone()
    chains_to_delete = []
    for model in cut_structure.st:
        for chain in model:
            chains_to_delete.append(chain.name)

    for chain_name in chains_to_delete:
        del cut_structure.st[0][chain_name]

    for model in structure.st:
        for chain in model:
            new_chain = gemmi.Chain(chain.name)
            for res in chain:
                if (chain.name == 'A') and (res.seqid.num == 60):
                    continue
                new_chain.add_residue(res)
            cut_structure.st[0].add_chain(new_chain) 
    return cut_structure

@pytest.fixture
def known_rfree() -> TStructure:
    return ...

@pytest.fixture
def known_molprobity() -> TStructure:
    return ...