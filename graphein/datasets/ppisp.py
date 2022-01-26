"""Class for working with the PPISP Dataset"""
# Graphein
# Author: Arian Jamasb <arian@jamasb.io>
# License: MIT
# Code Repository: https://github.com/a-r-j/graphein
import logging
import multiprocessing
import pickle
import traceback
from functools import lru_cache, partial
from pathlib import Path
from typing import List, Optional, Tuple

import networkx as nx
import numpy as np
import pandas as pd
from Bio.PDB import PDBList

from graphein.protein.config import ProteinGraphConfig
from graphein.protein.graphs import construct_graph

log = logging.getLogger(__name__)


class PPISP:
    def __init__(
        self,
        df: Optional[pd.DataFrame] = None,
        config: ProteinGraphConfig = ProteinGraphConfig(),
        num_cores: int = 32,
    ):
        self.ROOT_DIR: Path = (
            Path(__file__).parent.parent.parent.resolve()
            / "datasets"
            / "ppisp"
        )
        self._NUM_CORES: int = num_cores

        if df is not None:
            self.df: pd.DataFrame = self.parse_dataset()
        else:
            self.df: pd.DataFrame = self._load_dataset()

        self.chain_list: List[str] = self.get_chains()
        self.pdb_list: List[str] = self.get_pdbs()
        self.pssms: List[np.array] = self._load_pssm()
        self.node_labels: List[np.array] = self._load_node_labels()
        self.config: ProteinGraphConfig = config

        self.bad_pdbs: List[str] = []

        if self.config is not None:
            self.graphs: List[nx.Graph] = self.construct_graphs()

    def _load_dataset(self) -> pd.DataFrame:
        file_path = self.ROOT_DIR / "deepppisp_clean.csv"
        return pd.read_csv(file_path)

    def get_chains(self) -> List[str]:
        return list(self.df["chains"])

    def get_pdbs(self) -> List[str]:
        return list(self.df["pdb_code"])

    @lru_cache
    def _load_pssm(self) -> List[np.array]:
        """
        Loads PSSM (position-specific scoring matrix) features from pickled files

        :return: List of PSSMs for each protein (len protein x 20)
        :rtype: List[np.array]
        """
        with open(self.ROOT_DIR / "dset186_pssm_data.pkl", "rb") as f:
            dset_186_pssms = pickle.load(f)

        with open(self.ROOT_DIR / "dset164_pssm_data.pkl", "rb") as f:
            dset_164_pssms = pickle.load(f)

        with open(self.ROOT_DIR / "dset72_pssm_data.pkl", "rb") as f:
            dset_72_pssms = pickle.load(f)

        pssms = dset_186_pssms + dset_164_pssms + dset_72_pssms
        return [np.array(p) for p in pssms]

    @lru_cache
    def _load_node_labels(self) -> List[np.array]:
        """
        Loads node labels from pickled files

        :return: List of node labels for each protein
        :rtype: List[np.array]
        """
        with open(self.ROOT_DIR / "dset186_label.pkl", "rb") as f:
            dset186_labels = pickle.load(f)

        with open(self.ROOT_DIR / "dset164_label.pkl", "rb") as f:
            dset164_labels = pickle.load(f)

        with open(self.ROOT_DIR / "dset72_label.pkl", "rb") as f:
            dset72_labels = pickle.load(f)

        labels = dset186_labels + dset164_labels + dset72_labels
        return [np.array(l) for l in labels]

    @lru_cache
    def _load_sequences(self) -> List[np.array]:
        """
        Loads protein sequences from pickled files

        :return: List of node labels for each protein
        :rtype: List[np.array]
        """
        with open(self.ROOT_DIR / "dset186_sequence_data.pkl", "rb") as f:
            dset186_seqs = pickle.load(f)

        with open(self.ROOT_DIR / "dset164_sequence_data.pkl", "rb") as f:
            dset164_seqs = pickle.load(f)

        with open(self.ROOT_DIR / "dset72_sequence_data.pkl", "rb") as f:
            dset72_seqs = pickle.load(f)

        seqs = dset186_seqs + dset164_seqs + dset72_seqs
        return [np.array(s) for s in seqs]

    def split_data(self, strategy: str):
        if strategy == "deep_ppisp":
            train_df = self.df.loc[self.df["train"] == 1.0]
            test_df = self.df.loc[self.df["train"] == 0.0]
        else:
            raise ValueError(
                f"Unsupported splitting strategy: {strategy}. Use one of: 'deep_ppisp'"
            )
        self.train = PPISP(train_df)
        self.test = PPISP(test_df)

    def download_pdbs(self, path: str):
        """
        Downloads dataset PDBs to a specified directories

        :param path: Path to desired output location
        :type path: str
        """
        pdbl = PDBList()
        pdbl.download_pdb_files(pdb_codes=self.pdb_list, pdir=path)

    def __len__(self) -> int:
        """Returns length of the dataset

        :returns: Dataset length
        :rtype: int
        """
        return len(self.df)

    def construct_graphs(self) -> List[nx.Graph]:
        """
        Constructs graphs for protein chains in the dataset.

        :param config: Config specifying protein graph construction parameters
        :type config: graphein.protein.ProteinGraphConfig
        :return: List of protein structure graphs
        :rtype: List[nx.Graph]
        """
        pool = multiprocessing.Pool(self._NUM_CORES)
        graph_list = list(
            pool.map(
                self._graph_constructor,
                [
                    (pdb, self.chain_list[i])
                    for i, pdb in enumerate(self.pdb_list)
                ],
            )
        )
        pool.close()
        pool.join()
        return graph_list

    def _graph_constructor(self, args: Tuple[str, str]):
        """
        Partialed graph constructor for multiprocessing

        :param args: Tuple of pdb code and chain to build graph of
        :type args: Tuple[str, str]
        :return: Protein structure graph
        :rtype: nx.Graph
        """
        log.info(
            f"Constructing graph for: {args[0]}. Chain selection: {args[1]}"
        )
        func = partial(construct_graph, config=self.config)
        try:
            result = func(pdb_code=args[0], chain_selection=args[1])
            return result
        except Exception as ex:
            log.info(
                f"Graph construction error (PDB={args[0]})! {traceback.format_exc()}"
            )
            log.info(ex)
            self.bad_pdbs.append(args[0])
            return None


if __name__ == "__main__":
    c = PPISP()
