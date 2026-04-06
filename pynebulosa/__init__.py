"""pyNebulosa: Single-cell data visualization using kernel gene-weighted density estimation.

Recovers single-cell gene expression signals by kernel density estimation,
addressing the dropout/sparsity problem in single-cell RNA-seq data.
"""

from pynebulosa._kde import calculate_density, wkde2d
from pynebulosa._plotting import plot_density

__all__ = ["plot_density", "calculate_density", "wkde2d"]
__version__ = "0.1.0"
