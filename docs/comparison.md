# Comparison with other tools

## Overview

| Tool | Approach | Use case |
|------|----------|----------|
| `scanpy.tl.embedding_density` | Unweighted KDE per cell group | Cell density by category |
| pyUCell | Rank-based scoring + KNN smoothing | Gene signature scoring |
| **Nebulosa** | Gene-weighted KDE on embeddings | Recovering dropout signal for visualization |

## scanpy.tl.embedding_density

Scanpy's `embedding_density` computes a standard (unweighted) kernel density
estimate for cells in each group/category. It answers: *"where are cells of
group X concentrated?"* It does **not** incorporate gene expression as weights.

Nebulosa, by contrast, uses gene expression values as weights in the KDE. This
means the density at each point reflects both the local cell concentration
**and** the expression level of the gene, effectively recovering signal lost
to dropout.

## pyUCell

[pyUCell](https://pyucell.readthedocs.io/) computes gene signature scores
using a rank-based approach (UCell), optionally smoothed over a KNN graph.
It is designed for **scoring** cells on multi-gene signatures and returns
a per-cell score.

Nebulosa focuses on **visualization**: it produces density-weighted embedding
plots that reveal expression patterns obscured by sparsity. While pyUCell
smooths scores over a graph, Nebulosa smooths expression over the 2D embedding
space using kernel density estimation.

## When to use Nebulosa

- You want to **visualize** where a gene or a small set of genes is expressed
  on an embedding (UMAP, t-SNE), but standard feature plots are too sparse
  to be informative.
- You want to find **co-expressing cell populations** using joint density of
  multiple markers.
- You want a direct **weighted KDE** approach rather than graph-based smoothing.
