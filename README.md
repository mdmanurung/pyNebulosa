# Nebulosa

[![PyPI](https://img.shields.io/pypi/v/nebulosa.svg)](https://pypi.org/project/nebulosa/)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

Single-cell data visualization using kernel gene-weighted density estimation.

## Motivation

Due to the sparsity observed in single-cell data (e.g. RNA-seq, ATAC-seq), the
visualization of cell features (e.g. gene, peak) is frequently affected and
unclear, especially when it is overlaid with clustering to annotate cell types.
Nebulosa aims to recover the signal from dropped-out features by incorporating
the similarity between cells allowing a "convolution" of the cell features.

## Installation

```bash
pip install nebulosa
```

## Quick start

```python
import scanpy as sc
import nebulosa as nb

adata = sc.datasets.pbmc3k_processed()

# Single gene density plot
nb.plot_density(adata, "NKG7")

# Multiple genes
nb.plot_density(adata, ["MS4A1", "NKG7", "CST3"])

# Joint density to identify co-expressing cells
nb.plot_density(adata, ["NKG7", "GNLY"], joint=True)
```

## How it differs from similar tools

| Tool | What it does |
|------|-------------|
| `scanpy.tl.embedding_density` | Cell density per category (not gene-weighted) |
| pyUCell | Rank-based gene signature scoring with KNN smoothing |
| **Nebulosa** | Gene-weighted KDE visualization on embeddings |

## API

| Function | Purpose |
|----------|---------|
| `nb.plot_density(adata, features)` | Main visualization function |
| `nb.calculate_density(weights, coords)` | Compute density values programmatically |
| `nb.wkde2d(x, y, w)` | Low-level weighted 2D KDE |

See the [tutorial notebook](vignettes/nebulosa_tutorial.ipynb) for a full walkthrough.

## Citation

If you use Nebulosa, please cite:

> Alquicira-Hernandez, J., Powell, J.E. Nebulosa recovers single-cell gene
> expression signals by kernel density estimation.
> *Bioinformatics*, 37(16), 2485-2487, 2021.
> https://doi.org/10.1093/bioinformatics/btab003
