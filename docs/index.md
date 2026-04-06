# pyNebulosa

**Single-cell data visualization using kernel gene-weighted density estimation.**

Due to the sparsity observed in single-cell data (e.g. RNA-seq, ATAC-seq), the
visualization of cell features (e.g. gene, peak) is frequently affected and
unclear, especially when overlaid with clustering to annotate cell types.
pyNebulosa recovers the signal from dropped-out features by incorporating the
similarity between cells, allowing a "convolution" of the cell features via
weighted kernel density estimation.

## Quick start

```python
import scanpy as sc
import pynebulosa as nb

adata = sc.datasets.pbmc3k_processed()

# Single gene density plot
nb.plot_density(adata, "NKG7")

# Multiple genes
nb.plot_density(adata, ["MS4A1", "NKG7", "CST3"])

# Joint density to identify co-expressing cells
nb.plot_density(adata, ["NKG7", "GNLY"], joint=True)
```

```{toctree}
:maxdepth: 2
:caption: Contents

installation
tutorial
api
comparison
```

## Citation

If you use pyNebulosa, please cite:

> Alquicira-Hernandez, J., Powell, J.E. Nebulosa recovers single-cell gene
> expression signals by kernel density estimation.
> *Bioinformatics*, 37(16), 2485-2487, 2021.
> [doi:10.1093/bioinformatics/btab003](https://doi.org/10.1093/bioinformatics/btab003)
