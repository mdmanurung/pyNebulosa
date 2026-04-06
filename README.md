# pyNebulosa

[![PyPI](https://img.shields.io/pypi/v/pynebulosa.svg)](https://pypi.org/project/pynebulosa/)
[![Documentation](https://readthedocs.org/projects/pynebulosa/badge/?version=latest)](https://pynebulosa.readthedocs.io)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

**Single-cell data visualization using kernel gene-weighted density estimation.**

pyNebulosa is a Python package (ported from the
[original R/Bioconductor package](https://bioconductor.org/packages/Nebulosa))
that recovers gene expression signals lost to dropout in single-cell data by
using weighted kernel density estimation on low-dimensional embeddings. It
integrates natively with [scanpy](https://scanpy.readthedocs.io) and
[AnnData](https://anndata.readthedocs.io).

## Motivation

Standard feature plots in single-cell analysis are plagued by sparsity — many
cells show zero expression even for marker genes, making it hard to see
biological signal. pyNebulosa addresses this by computing a gene-weighted kernel
density estimate (KDE) on the embedding space, effectively smoothing expression
over neighbouring cells.

| Standard feature plot | pyNebulosa density plot |
|:-----:|:-----:|
| Sparse, hard to interpret | Smooth, recovers signal |

![Standard vs Density](docs/images/standard_vs_density.png)

## R vs Python

pyNebulosa was originally an R/Bioconductor package (Nebulosa) for Seurat and
SingleCellExperiment objects. This Python port provides an equivalent API for
the scanpy/AnnData ecosystem, producing identical density estimates.

![R vs Python comparison](docs/images/r_vs_python.png)

<table>
<tr><th>R (Seurat)</th><th>Python (scanpy)</th></tr>
<tr>
<td>

```r
library(Nebulosa)
library(Seurat)
pbmc <- pbmc3k.SeuratData::pbmc3k.final
plot_density(pbmc, "NKG7")
plot_density(pbmc, c("NKG7", "GNLY"),
             joint = TRUE)
```

</td>
<td>

```python
import scanpy as sc
import pynebulosa as nb
adata = sc.datasets.pbmc3k_processed()
nb.plot_density(adata, "NKG7")
nb.plot_density(adata, ["NKG7", "GNLY"],
                joint=True)
```

</td>
</tr>
</table>

## Installation

```bash
pip install pynebulosa
```

With scanpy (recommended):

```bash
pip install "pynebulosa[scanpy]"
```

## Quick start

```python
import scanpy as sc
import pynebulosa as nb

adata = sc.datasets.pbmc3k_processed()

# Single gene density plot
nb.plot_density(adata, "NKG7")
```

### Multi-marker panels

Visualize multiple markers at once to identify cell populations:

```python
nb.plot_density(adata, ["MS4A1", "NKG7", "CST3", "PPBP"], ncols=4)
```

![Multi-marker density](docs/images/multi_marker.png)

### Joint density

Find co-expressing cell populations by computing the product of individual
gene densities:

```python
nb.plot_density(adata, ["NKG7", "GNLY"], joint=True)
```

![Joint density](docs/images/joint_density.png)

## API

| Function | Purpose |
|----------|---------|
| `nb.plot_density(adata, features)` | Main visualization function |
| `nb.calculate_density(weights, coords)` | Compute density values programmatically |
| `nb.wkde2d(x, y, w)` | Low-level weighted 2D KDE |

### Key parameters

| Parameter | Description | Default |
|-----------|-------------|---------|
| `joint` | Compute joint density across features | `False` |
| `reduction` | Embedding key in `adata.obsm` | Auto-detect (`X_umap`) |
| `method` | `"wkde"` (custom) or `"ks"` (scipy) | `"wkde"` |
| `adjust` | Bandwidth adjustment factor | `1.0` |
| `cmap` | Matplotlib colormap | `"viridis"` |
| `layer` | Expression layer | `adata.X` |

## How it differs from similar tools

| Tool | Approach | Use case |
|------|----------|----------|
| `scanpy.tl.embedding_density` | Unweighted KDE per category | Where are cells of group X? |
| pyUCell | Rank-based scoring + KNN smoothing | Gene signature scoring |
| **pyNebulosa** | Gene-weighted KDE on embeddings | Recovering dropout signal |

## Documentation

- [Full documentation](https://pynebulosa.readthedocs.io) (ReadTheDocs)
- [Tutorial notebook](vignettes/nebulosa_tutorial.ipynb) (PBMC 3K walkthrough)

## Contributing

We welcome contributions of all kinds! If you'd like to contribute to
pyNebulosa, please:

1. Fork the repository and create a feature branch.
2. Make your changes and add tests where appropriate.
3. Ensure all tests pass with `pytest tests/ -v`.
4. Submit a pull request describing your changes.

Please report bugs and feature requests on the
[issue tracker](https://github.com/mdmanurung/pyNebulosa/issues).

## License

pyNebulosa is licensed under the [MIT License](LICENSE).

## Citation

If you use pyNebulosa, please cite:

> Alquicira-Hernandez, J., Powell, J.E. Nebulosa recovers single-cell gene
> expression signals by kernel density estimation.
> *Bioinformatics*, 37(16), 2485-2487, 2021.
> [doi:10.1093/bioinformatics/btab003](https://doi.org/10.1093/bioinformatics/btab003)
