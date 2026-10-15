"""
Unit test that symrmsds are working correctly
"""
import pytest

import numpy as np

from rocket_testing.types import TStructure, TAtomID, TRestraints
from rocket_testing.symrmsd import symmetric_rmsd

def test_symrmsd_self_all_atom_with_ligand(
        structure_with_ligand: TStructure, 
        restraints: TRestraints, 
        all_atom_selection_with_ligand: list[TAtomID],
        ):
    # Test the usual all atom RMSD ignoring atom altlocs
    structure_copy: TStructure = structure_with_ligand.clone()
    rmsd, distances = symmetric_rmsd(structure_with_ligand, structure_copy, all_atom_selection_with_ligand, {restraints.name: restraints})
    assert np.allclose(rmsd, 0.0)

def test_symrmsd_self_all_atom(structure: TStructure, all_atom_selection: list[TAtomID]):
    # Test the usual all atom RMSD ignoring atom altlocs
    structure_copy: TStructure = structure.clone()
    rmsd, distances = symmetric_rmsd(structure, structure_copy, all_atom_selection)

    assert np.allclose(rmsd, 0.0)

def test_symrmsd_self_chain(structure: TStructure, chain_selection: list[TAtomID]):
    # Test all atom within chain RMSD
    structure_copy: TStructure = structure.clone()
    rmsd, distances = symmetric_rmsd(structure, structure_copy, chain_selection)

    assert np.allclose(rmsd, 0.0)

def test_symrmsd_self_res(structure: TStructure, residue_selection: list[TAtomID]):
    # Test symrmsd: res to res
    structure_copy: TStructure = structure.clone()
    rmsd, distances = symmetric_rmsd(structure, structure_copy, residue_selection)

    assert np.allclose(rmsd, 0.0)
    assert len(distances) == 1
    assert len(distances[('A', "63", "LYS")]) == 9

def test_symrmsd_self_atom(structure: TStructure, atom_selection: list[TAtomID]):
     
    # Test symrmsd: atom to atom
    structure_copy: TStructure = structure.clone()
    rmsd, distances = symmetric_rmsd(structure, structure_copy, atom_selection)

    assert np.allclose(rmsd, 0.0)
    assert len(distances) == 1
    assert len(distances[('A', "63", "LYS")]) == 1

def test_symrmsd_self_sym_atom_swap(structure, structure_swapped: TStructure, residues_selection_swapped: list[TAtomID]):
    # Test the core functionalty - two chemically symmetric atoms swapped in place with no increase in RMSD! 
    structure_copy: TStructure = structure.clone()
    rmsd, distances = symmetric_rmsd(structure, structure_swapped, residues_selection_swapped)

    assert np.allclose(rmsd, 0.0)
    assert len(distances) == 1
    assert len(distances[('A', "62", "GLU")]) == 9
    
def test_symrmsd_self_partial(structure: TStructure, cut_structure: TStructure, residues_selection: list[TAtomID]):

    # Test the expected behaviour when part of the selection is missing i.e. skipping that subset
    rmsd, distances = symmetric_rmsd(structure, cut_structure, residues_selection)

    assert len(distances) == 2
    assert np.allclose(rmsd, 0.0)

    
@pytest.mark.xfail()
def test_symrmsd_self_incompatible(structure: TStructure, cut_structure: TStructure, incompatible_selection: list[TAtomID]):
    # TStructure to itself: TStructure to incompatible structure - should fail

    rmsd, distances = symmetric_rmsd(structure, cut_structure, incompatible_selection)

def test_symrmsd_self_offset(structure: TStructure, offset_structure: TStructure, residue_selection: list[TAtomID]):

    # TStructure to self: offset - should be distance of offset!
    structure_copy: TStructure = structure.clone()
    
    rmsd, distances = symmetric_rmsd(structure_copy, offset_structure, residue_selection, align=False)

    assert np.allclose(rmsd, np.sqrt(np.sum(np.square([0.1, 0.1, 0.1]))))

def test_symrmsd_self_offset_align(structure: TStructure, offset_structure: TStructure, residue_selection: list[TAtomID]):
    # TStructure to self: offset - should be zero as offset is alignable!
    structure_copy: TStructure = structure.clone()
    
    rmsd, distances = symmetric_rmsd(structure_copy, offset_structure, residue_selection, align=True)

    assert np.allclose(rmsd, 0.0)

@pytest.mark.skip('Not needed - test case covered by offset!')
def test_symrmsd_self_perturbed(structure: TStructure, perturbed_structure: TStructure, residue_selection: list[TAtomID]):
    # TStructure to self pertrubed - should be small non-zero

    rmsd, distances = symmetric_rmsd(structure, perturbed_structure, residue_selection)
 
    assert 0 < rmsd
    