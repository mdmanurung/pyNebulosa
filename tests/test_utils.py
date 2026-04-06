"""Tests for pynebulosa._utils — validation and data extraction helpers."""

import anndata
import numpy as np
import pytest

from pynebulosa._utils import _get_embeddings, _get_feature_data, _validate_dims


class TestValidateDims:
    def test_valid_two_dims(self):
        _validate_dims((0, 1))  # Should not raise

    def test_three_dims_raises(self):
        with pytest.raises(ValueError, match="Exactly 2"):
            _validate_dims((0, 1, 2))

    def test_one_dim_raises(self):
        with pytest.raises(ValueError, match="Exactly 2"):
            _validate_dims((0,))

    def test_empty_raises(self):
        with pytest.raises(ValueError, match="Exactly 2"):
            _validate_dims(())


class TestGetEmbeddings:
    def test_auto_detects_umap(self, adata_simple):
        emb, labels = _get_embeddings(adata_simple, None, (0, 1))
        assert emb.shape == (adata_simple.n_obs, 2)
        assert "UMAP" in labels[0]

    def test_fallback_to_last_obsm_key(self, adata_simple):
        """Covers line 73: fallback when X_umap is not present."""
        adata = adata_simple.copy()
        adata.obsm.pop("X_umap")
        # Only X_pca remains
        emb, labels = _get_embeddings(adata, None, (0, 1))
        assert emb.shape == (adata.n_obs, 2)
        assert "PCA" in labels[0]

    def test_empty_obsm_raises(self):
        """Covers line 65: no reductions at all."""
        adata = anndata.AnnData(X=np.ones((10, 5)))
        with pytest.raises(KeyError, match="No reductions"):
            _get_embeddings(adata, None, (0, 1))

    def test_explicit_reduction(self, adata_simple):
        emb, labels = _get_embeddings(adata_simple, "X_pca", (0, 1))
        assert emb.shape == (adata_simple.n_obs, 2)
        assert "PCA" in labels[0]

    def test_missing_reduction_raises(self, adata_simple):
        with pytest.raises(KeyError, match="not found"):
            _get_embeddings(adata_simple, "X_nonexistent", (0, 1))

    def test_out_of_range_dim_raises(self, adata_simple):
        with pytest.raises(ValueError, match="not present"):
            _get_embeddings(adata_simple, "X_umap", (0, 99))

    def test_dim_labels_correct(self, adata_simple):
        _, labels = _get_embeddings(adata_simple, "X_umap", (0, 1))
        assert labels == ("UMAP 1", "UMAP 2")

    def test_non_default_dims(self, adata_simple):
        emb, labels = _get_embeddings(adata_simple, "X_pca", (2, 5))
        assert emb.shape == (adata_simple.n_obs, 2)
        assert labels == ("PCA 3", "PCA 6")


class TestGetFeatureData:
    def test_gene_expression_sparse(self, adata_simple):
        result = _get_feature_data(adata_simple, ["marker1"], None)
        assert result.shape == (adata_simple.n_obs, 1)
        assert np.all(np.isfinite(result))

    def test_gene_expression_dense(self, adata_dense):
        result = _get_feature_data(adata_dense, ["marker1"], None)
        assert result.shape == (adata_dense.n_obs, 1)

    def test_multiple_genes(self, adata_simple):
        result = _get_feature_data(adata_simple, ["marker1", "marker2"], None)
        assert result.shape == (adata_simple.n_obs, 2)

    def test_metadata_feature(self, adata_simple):
        adata_simple.obs["score"] = np.random.default_rng(0).uniform(size=adata_simple.n_obs)
        result = _get_feature_data(adata_simple, ["score"], None)
        assert result.shape == (adata_simple.n_obs, 1)

    def test_missing_feature_raises(self, adata_simple):
        with pytest.raises(ValueError, match="not found"):
            _get_feature_data(adata_simple, ["NONEXISTENT"], None)

    def test_invalid_layer_raises(self, adata_simple):
        with pytest.raises(ValueError, match="Layer"):
            _get_feature_data(adata_simple, ["marker1"], "bad_layer")

    def test_valid_layer(self, adata_simple):
        """Covers line 152: extraction from a named layer."""
        adata_simple.layers["raw_counts"] = adata_simple.X.copy()
        result = _get_feature_data(adata_simple, ["marker1"], "raw_counts")
        assert result.shape == (adata_simple.n_obs, 1)

    def test_mixed_obs_and_var_features(self, adata_simple):
        """Covers line 164: feature in obs but not in var."""
        adata_simple.obs["score"] = np.random.default_rng(0).uniform(size=adata_simple.n_obs)
        # marker1 is in var, score is in obs
        result = _get_feature_data(adata_simple, ["marker1", "score"], None)
        assert result.shape == (adata_simple.n_obs, 2)
