import pytest
from hypothesis import given
from hypothesis.strategies import text

from pathlib import Path

from graphein.rna.edges import (
    add_all_dotbracket_edges,
    add_base_pairing_interactions,
    add_phosphodiester_bonds,
)
from graphein.rna.graphs import (
    RNA_BASES,
    SUPPORTED_DOTBRACKET_NOTATION,
    construct_rna_graph,
)

TEST_SEQUENCE = "UUGGAGUACACAACCUGUACACUCUUUC"
TEST_DOTBRACKET = "..(((((..(((...)))..)))))..."


def test_construct_rna_graph():
    g = construct_rna_graph(
        dotbracket=TEST_DOTBRACKET,
        sequence=TEST_SEQUENCE,
        edge_construction_funcs=[
            add_base_pairing_interactions,
            add_phosphodiester_bonds,
        ],
    )

    h = construct_rna_graph(
        dotbracket=TEST_DOTBRACKET,
        sequence=TEST_SEQUENCE,
        edge_construction_funcs=[add_all_dotbracket_edges],
    )

    assert g.edges() == h.edges()
    assert g.nodes() == h.nodes()

    # Check number of nodes and edges
    assert len(g.nodes()) == len(TEST_SEQUENCE)
    assert len(g.nodes()) == len(TEST_DOTBRACKET)

    # Check node features are in alphabet
    for n, d in g.nodes(data=True):
        assert d["dotbracket_symbol"] in SUPPORTED_DOTBRACKET_NOTATION
        assert d["nucleotide"] in RNA_BASES

    # Count edges for each type
    phosphodiesters = 0
    bp = 0
    for u, v, d in g.edges(data=True):
        if d["attr"] == "phosphodiester_bond":
            phosphodiesters += 1
        elif d["attr"] == "base_pairing":
            bp += 1

    assert phosphodiesters == len(TEST_SEQUENCE) - 1
    assert bp == TEST_DOTBRACKET.count("(")


def test_pdb_rna_graph():
    rna = construct_rna_graph(
        pdb_code="2jyf"
    )

def test_construct_graph():
    """Example-based test that graph construction works correctly.

    Uses 4hhb PDB file as an example test case.
    """
    file_path = Path(__file__).parent / "test_data/2jyf.pdb"
    G = construct_rna_graph(pdb_path=str(file_path))


if __name__ == "__main__":
    test_pdb_rna_graph()
    test_construct_graph()
