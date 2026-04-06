"""Shared test fixtures for pyNebulosa tests."""

import anndata
import numpy as np
import pytest
import scipy.sparse as sp


@pytest.fixture
def adata_simple():
    """Small synthetic AnnData with known expression patterns and UMAP embedding.

    200 cells, 5 genes. Gene 'marker1' is expressed in cells 0-49 (cluster in
    top-right of UMAP). Gene 'marker2' is expressed in cells 100-149 (cluster
    in bottom-left).
    """
    rng = np.random.default_rng(42)
    n_cells = 200
    n_genes = 5

    # Create UMAP-like embedding with 4 clusters
    umap = np.zeros((n_cells, 2))
    # Cluster 1: top-right
    umap[:50] = rng.normal(loc=[3, 3], scale=0.5, size=(50, 2))
    # Cluster 2: top-left
    umap[50:100] = rng.normal(loc=[-3, 3], scale=0.5, size=(50, 2))
    # Cluster 3: bottom-left
    umap[100:150] = rng.normal(loc=[-3, -3], scale=0.5, size=(50, 2))
    # Cluster 4: bottom-right
    umap[150:200] = rng.normal(loc=[3, -3], scale=0.5, size=(50, 2))

    # Create expression matrix
    X = rng.poisson(lam=0.5, size=(n_cells, n_genes)).astype(float)
    # marker1: high expression in cluster 1
    X[:50, 0] = rng.poisson(lam=5, size=50).astype(float)
    # marker2: high expression in cluster 3
    X[100:150, 1] = rng.poisson(lam=5, size=50).astype(float)

    gene_names = ["marker1", "marker2", "gene3", "gene4", "gene5"]
    cell_names = [f"cell_{i}" for i in range(n_cells)]

    adata = anndata.AnnData(
        X=sp.csr_matrix(X),
        obs={"cluster": (["A"] * 50 + ["B"] * 50 + ["C"] * 50 + ["D"] * 50)},
    )
    adata.var_names = gene_names
    adata.obs_names = cell_names
    adata.obsm["X_umap"] = umap
    adata.obsm["X_pca"] = rng.normal(size=(n_cells, 10))

    return adata


@pytest.fixture
def adata_dense(adata_simple):
    """Same as adata_simple but with dense X matrix."""
    adata = adata_simple.copy()
    adata.X = np.asarray(adata.X.toarray())
    return adata
