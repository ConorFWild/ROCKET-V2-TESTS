"""
Functions for making atomic selections from structural models and 
calculating graph symmetric RMSDs
"""

from functools import cache
import json
from typing import Optional, Iterable, Mapping

import numpy.typing as npt
# from spyrmsd import rmsd
import gemmi
from Bio.Align import Alignment
from Bio import Align
import networkx
import numpy as np
import networkx.algorithms.isomorphism as iso


from .selections import get_atoms, get_atom_ids
from .types import TStructure, TResID, TAtomID, AtomID
from . import constants

@cache
def get_monomer_graphs():
    with open(constants.MONOMER_JSON, 'r') as f:
        monomer_json = json.load(f)

    monomer_graphs = {}
    for monomer, monomer_data in monomer_json.items():
        monomer_graphs[monomer] = networkx.cytoscape_graph(monomer_data)
        monomer_graphs[monomer].remove_nodes_from([x for x in monomer_graphs[monomer].nodes if monomer_graphs[monomer].nodes[x]['Z'] == 1])
        monomer_graphs[monomer].remove_nodes_from(['OXT'])  # Graphs are generated from CCP3 MONLIB, which has OXTs

    return monomer_graphs

def get_sequence(ref_structure, ref_chain):
    # Get the sequences AS ATTESTED IN THE ACTUAL MODEL
    seq = []
    index_to_res = {}
    j = 0
    for residue in ref_structure.st[0][ref_chain]:
        if residue.name not in constants.RESIDUES:
            continue
        seq.append(residue.name)
        index_to_res[j] = residue
        j += 1

    seq_one_letter = [gemmi.find_tabulated_residue(_res).one_letter_code for _res in seq]
    return ''.join(seq_one_letter), index_to_res

def iterate_alignments(alignments, ref_structure, mov_structure, index_to_ref_res, index_to_mov_res):
    ref_indicies = alignments[0].indices[0]
    mov_indicies = alignments[0].indices[1]

    # Iterate the aligments - when there is a gap index is -1 so we just continue!
    for ref_index, mov_index in zip(ref_indicies, mov_indicies):
        match ref_index, mov_index:
            case (-1, x):
                ...
            case (x, -1):
                ...
            case (x, y):
                ref_res = index_to_ref_res[x]
                mov_res = index_to_mov_res[y]
                yield ref_res, mov_res


def get_res_atoms_in_selection(res, chain, atom_selection):
    atoms = {}
    atom_names = set([atom.name for atom in res])
    for atom_name in atom_names:

        atom = res[atom_name][0]
        atom_id = AtomID(chain, str(res.seqid.num), atom.name, atom.altloc)
        if atom_id in atom_selection:
            atoms[atom.name] = atom

    return atoms

def get_coords_from_atoms(atoms):
    
    return np.array(
        [
            [atom.pos.x, atom.pos.y, atom.pos.z]
            for atom 
            in atoms
        ])

def get_res_by_id(st, resid):
    return st.st[0][resid[0]][resid[1]]

def get_atom_distances(isomorphism, ref_atoms, mov_atoms):

    distances = []
    for ref_atom_name, ref_atom in ref_atoms.items():
        isomorphic_atom_name = isomorphism[ref_atom_name]
        distances.append(ref_atom.pos.dist(mov_atoms[isomorphic_atom_name].pos))

    return distances

def symmetric_rmsd(
        ref_structure: TStructure, 
        mov_structure: TStructure, 
        atom_selection: Optional[list[TAtomID]]=None,
        chain_map = None,
        align=True
        ) -> tuple[float, npt.ArrayLike]:
    """
    Returns the RMSD over the selection and the set of distances used
    """

    # Fail early if some basic checks fail
    assert isinstance(ref_structure, TStructure)
    assert isinstance(mov_structure, TStructure)
    assert isinstance(chain_map, Mapping) | isinstance(chain_map, type(None))
    assert isinstance(align, bool)        
    assert isinstance(atom_selection, Iterable)

    # Default selection of all atoms
    if not atom_selection:
        atom_selection = get_atom_ids(ref_structure)

    # Rename end oxygens in structures - this should probably create clones for internal use...
    for atom in get_atoms(ref_structure):
        if atom.name == 'OXT':
            atom.name = 'O'
    for atom in get_atoms(mov_structure):
        if atom.name == 'OXT':
            atom.name = 'O'

    # Get the protein monomer graphs - these canonically define ordering of residue atoms
    monomer_graphs = get_monomer_graphs()

    # Get isomorphisms for each monomer - only use heavy atoms!
    isomorphisms = {}
    for monomer, monomer_graph in monomer_graphs.items():
        matcher = iso.GraphMatcher(monomer_graph, monomer_graph, node_match=iso.categorical_node_match('Z', 0))
        isomorphisms[monomer] = []
        for isomorphism in matcher.isomorphisms_iter():
            isomorphisms[monomer].append(isomorphism)

    # Check chains match if no chain map provided
    if chain_map is None:
        chain_map = {}
        for ref_chain in ref_structure.st[0]:
            try:
                mov_chain = mov_structure.st[0][ref_chain.name]
                chain_map[ref_chain.name] = ref_chain.name
            except:
                raise Exception(
                    f'Chain {ref_chain.name} in ref structure from {ref_structure.path} not present in mov structure {mov_structure.path}\n'
                    # f'Move chains are {[x.name for x in mov_structure.st[0]]}\n'
                    f'Ref chains are {[_chain.name for _chain in ref_structure.st[0]]} while mov chains are {[_mov_chain.name for _mov_chain in mov_structure.st[0]]}'
                    )

    # Get the distances between all atoms in the selection
    min_distances = {}
    for ref_chain_name, mov_chain_name in chain_map.items():
        # Align structures by residue
        if align:
            aligned_mov_structure = mov_structure.clone()
            ref_pol = ref_structure.st[0][ref_chain_name].get_polymer()
            mov_pol = mov_structure.st[0][mov_chain_name].get_polymer()
            ptype = ref_pol.check_polymer_type()
            sup = gemmi.calculate_superposition(ref_pol, mov_pol, ptype, gemmi.SupSelect.CaP)
            aligned_mov_structure.st[0][mov_chain_name].get_polymer().transform_pos_and_adp(sup.transform)
        else:
            aligned_mov_structure = mov_structure

        # Get sequence alignment
        aligner = Align.PairwiseAligner()
        ref_seq, index_to_ref_res = get_sequence(ref_structure, ref_chain_name)
        mov_seq, index_to_mov_res = get_sequence(aligned_mov_structure, mov_chain_name)
        if len(ref_seq) | len(mov_seq) == 0:  # Skip unalignable chains: i.e. water chains
            continue
        alignments = aligner.align(
            ref_seq, 
            mov_seq,
            )

        # Factorize by residue - the symmetric motifs of residues are independent!
        # This means the sum of the minima of each residue is the -global- minima!
        # Very handy! The isomorphism graph of the whole protein is huuuuuuuge!
        for ref_residue, mov_residue in iterate_alignments(alignments, ref_structure, aligned_mov_structure, index_to_ref_res, index_to_mov_res):
   
            # Get selection atoms
            ref_atoms = get_res_atoms_in_selection(ref_residue, ref_chain_name, atom_selection)
            if len(ref_atoms) == 0:
                continue
            mov_atoms = {
                _ref_atom_name: mov_residue[_ref_atom_name][0] 
                for _ref_atom_name in ref_atoms}

            # Get all the possible graph isomorphisms
            delta_sets = []
            for isomorphism in isomorphisms[ref_residue.name]:
                atom_distances = get_atom_distances(isomorphism, ref_atoms, mov_atoms)
                delta_sets.append(atom_distances)
            min_distances[(ref_chain_name, str(ref_residue.seqid.num), ref_residue.name)] = min(delta_sets, key=lambda _delta_set: np.sum(_delta_set))

    # Get the final RMSD
    rmsd = np.sqrt(np.mean(np.square(np.concat([_delta for _delta in min_distances.values()], axis=None))))
    return rmsd, min_distances
    ...