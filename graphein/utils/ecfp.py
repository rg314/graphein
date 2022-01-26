"""Implementation of Extended-Connectivity Fingerprint algorithm for Protein Graphs (David Rogers and Mathew Hahn 2010).

This is based on the MorganFingerprints.cpp code from rdkit:
    https://github.com/rdkit/rdkit/blob/master/Code/GraphMol/Fingerprints/MorganFingerprints.cpp
and the ecfp.py code from pinky:
    https://github.com/ubccr/pinky/blob/1961bbf156e91022a04e829a41d47ad7bd5ca540/pinky/fingerprints/ecfp.py
"""
# Graphein
# Author: Arian Jamasb <arian@jamasb.io>
# License: MIT
# Code Repository: https://github.com/a-r-j/graphein

from copy import deepcopy
from typing import Dict, List, Optional

import networkx as nx
from bitarray import bitarray


class ECFP:
    def __init__(
            self,
            radius: int = 3,
            features: Optional[List[str]] = ["residue_number", "b_factor"]):
        self.radius = radius
        self.features = features

    @staticmethod
    def gen_hash(lst):
        return hash(tuple(lst))

    def get_invariants(self, g):
        """Generate initial atom identifiers using atomic invariants"""
        node_ids: Dict = {}
        for n, d in g.nodes(data=True):
            components: List = []
            components.append(round(d[feat]) for feat in self.features)
            node_ids[n] = self.gen_hash(components)

        return node_ids

    def ecfp(self, g: nx.Graph):
        node_ids = self.get_invariants(g)

        fp: Dict = {}
        for i in g.nodes():
            fp[i] = fp.get(i, 0) + 1

        node_index_map: Dict[str, int] = {n: i for i, n in enumerate(g.nodes)}
        index_node_map: Dict[str, int] = {i: n for i, n in enumerate(g.nodes)}
        edge_index_map: Dict[str, int] = {e: i for i, e in enumerate(g.edges)}
        index_edge_map: Dict[str, int] = {i: e for i, e in enumerate(g.edges)}

        neighbourhoods: List = []
        atom_neighbourhoods = [g.number_of_edges() * bitarray('0') for a in range(g.number_of_nodes())]
        dead_atoms = len(g.nodes()) * bitarray('0')

        for layer in range(1, self.radius+1):
            round_ids: Dict = {}
            round_atom_neighbourhoods = deepcopy(atom_neighbourhoods)
            neighborhoods_this_round: List = []

            for i, a in enumerate(g.nodes()):
                if dead_atoms[i]: continue

                nbsr: List = []
                for j, adj_node in enumerate([n for n in g.neighbors(a)]):
                    b_index = node_index_map[adj_node]
                    try:
                        bond_index = edge_index_map[(a, adj_node)]
                        round_atom_neighbourhoods[i][bond_index] = True
                        bond_type = list(g.get_edge_data(a, adj_node)["kind"])[0]
                    except:
                        bond_index = edge_index_map[(adj_node, a)]
                        round_atom_neighbourhoods[i][bond_index] = True
                        bond_type = list(g.get_edge_data(a, adj_node)["kind"])[0]

                    round_atom_neighbourhoods[i] |= atom_neighbourhoods[b_index]
                    nbsr.append((bond_type, node_ids[adj_node]))

                nbsr = sorted(nbsr)
                nbsr = [item for sublist in nbsr for item in sublist]
                nbsr.insert(0, node_ids[a])
                nbsr.insert(0, layer)

                round_ids[a] = self.gen_hash(nbsr)
                neighborhoods_this_roub  nd.append(
                    (round_atom_neighbourhoods[i], round_ids[a], i)
                )
                print(neighborhoods_this_round)
            for lst in neighborhoods_this_round:
                if lst[0] not in neighbourhoods:
                    fp[lst[1]] = fp.get(lst[1], 0) + 1
                    neighbourhoods.append(lst[0])
                else:
                    dead_atoms[lst[2]] = True

            node_ids = round_ids
            atom_neighbourhoods = deepcopy(round_atom_neighbourhoods)
        return fp


if __name__ == "__main__":
    from graphein.protein.edges.distance import add_hydrogen_bond_interactions
    from graphein.protein.graphs import construct_graph

    ecfp = ECFP(radius=4)
    g = construct_graph(pdb_code="3eiy", edge_construction_funcs=[add_hydrogen_bond_interactions])
    f = ecfp.ecfp(g)
    #f = ecfp.get_invariants(g)
