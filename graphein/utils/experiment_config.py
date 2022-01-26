from enum import Enum
from pathlib import Path
from typing import Optional

from pydantic import BaseModel

from graphein.grn.config import GRNGraphConfig
from graphein.ppi.config import PPIGraphConfig
from graphein.protein.config import ProteinGraphConfig, ProteinMeshConfig


class REPRESENTATION_MODE(Enum):
    PROTEIN_GRAPH = "protein_graph"
    PROTEIN_MESH = "protein_mesh"
    RNA_GRAPH = "rna_graph"
    PPI_GRAPH = "ppi_graph"
    GRN_GRAPH = "grn_graph"


class ExperimentConfig(BaseModel):
    representation_mode: REPRESENTATION_MODE
    protein_graph_config: Optional[ProteinGraphConfig] = None
    protein_mesh_config: Optional[ProteinMeshConfig] = None
    ppi_graph_config: Optional[PPIGraphConfig] = None
    grn_graph_config: Optional[GRNGraphConfig] = None

    out_dir: Path

    class Config:
        arbitrary_types_allowed: bool = True
