import networkx as nx


def add_sequence_to_nodes(n, d):
    """
    Maps UniProt ACC to UniProt ID. Retrieves sequence from UniProt and adds it to the node
    :param n: Graph node.
    :param d: Graph attribute dictionary.
    """
    from bioservices import HGNC, UniProt

    h = HGNC(verbose=False)
    u = UniProt(verbose=False)

    d["uniprot_ids"] = h.fetch("symbol", d["protein_id"])["response"]["docs"][
        0
    ]["uniprot_ids"]

    # Todo these API calls should probably be batched
    # Todo mapping with bioservices to support other protein IDs?

    for id in d["uniprot_ids"]:
        d[f"sequence_{id}"] = u.get_fasta_sequence(id)


if __name__ == "__main__":
    from functools import partial

    import matplotlib.pyplot as plt

    from graphein.ppi.config import PPIGraphConfig
    from graphein.ppi.edges import add_biogrid_edges, add_string_edges
    from graphein.ppi.graphs import compute_ppi_graph

    protein_list = [
        "CDC42",
        "CDK1",
        "KIF23",
        "PLK1",
        "RAC2",
        "RACGAP1",
        "RHOA",
        "RHOB",
    ]

    config = PPIGraphConfig()
    kwargs = config.kwargs

    g = compute_ppi_graph(
        protein_list=protein_list,
        edge_construction_funcs=[
            partial(add_string_edges, kwargs=kwargs),
            partial(add_biogrid_edges, kwargs=kwargs),
        ],
    )

    for n, d in g.nodes(data=True):
        add_sequence_to_nodes(n, d)

    print(nx.get_node_attributes(g, "sequence"))
    edge_colors = [
        "r"
        if g[u][v]["kind"] == {"string"}
        else "b"
        if g[u][v]["kind"] == {"biogrid"}
        else "y"
        for u, v in g.edges()
    ]

    print(nx.info(g))
    nx.draw(g, with_labels=True, edge_color=edge_colors)
    plt.show()
