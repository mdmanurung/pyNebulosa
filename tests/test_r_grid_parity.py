"""Check the R grid scale and R findInterval mapping separately from bandwidth."""
from pathlib import Path
import numpy as np
from pynebulosa._kde import wkde2d,_get_dens

def test_fixed_r_bandwidth_reference():
    root=Path(__file__).parent/'fixtures'/'r_grid'
    coords=np.loadtxt(root/'coordinates.csv',delimiter=',')
    w=np.loadtxt(root/'weights.csv',delimiter=',')
    h=np.loadtxt(root/'bandwidth.csv',delimiter=',')
    grid=wkde2d(coords[:,0],coords[:,1],w/w.sum()*len(w),h=h)
    np.testing.assert_allclose(grid['z'],np.loadtxt(root/'grid.csv',delimiter=','),atol=1e-12,rtol=1e-10)
    cells=_get_dens(coords,grid['x'],grid['y'],grid['z'])
    np.testing.assert_allclose(cells,np.loadtxt(root/'cells.csv',delimiter=','),atol=1e-12,rtol=1e-10)

def test_interval_boundaries_follow_r_findinterval():
    coords=np.array([[0.,0.],[.2,.2],[1.,1.],[1.9,1.9],[2.,2.]])
    grid=np.array([0.,1.,2.]);z=np.arange(9).reshape(3,3)
    np.testing.assert_array_equal(_get_dens(coords,grid,grid,z),[0,0,4,4,8])
