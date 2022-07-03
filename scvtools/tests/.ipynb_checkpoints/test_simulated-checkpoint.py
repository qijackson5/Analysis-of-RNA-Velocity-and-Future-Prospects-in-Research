import scvelo as scv
import numpy as np
import pytest

def test_high_var_subset():
    adata = scv.datasets.simulation(random_seed=1, n_vars=8)
    bdata = adata.copy()
    scv.pp.filter_and_normalize(adata, n_top_genes=5, subset_highly_variable=True)
    scv.pp.filter_and_normalize(bdata, n_top_genes=5, subset_highly_variable=False)
    scv.pp.pca(adata)
    scv.pp.pca(bdata)
    scv.pp.moments(adata, use_rep="pca")
    scv.pp.moments(bdata, use_rep="pca")
    scv.tl.velocity_graph(adata)
    scv.tl.velocity_graph(bdata)
    bdata._inplace_subset_var(bdata.var["highly_variable"])
    assert np.allclose(adata.layers["Ms"][0], bdata.layers["Ms"][0])
    assert np.allclose(adata.layers["velocity"][0], bdata.layers["velocity"][0])
    assert np.allclose(
        adata.uns["velocity_graph"].data[:5], bdata.uns["velocity_graph"].data[:5])

# testing whether the dataset is empty. If datasets are empty, entire notebook will not run!
def test_load_datasets():
    adata = scv.datasets.simulation(n_obs=500, t_max=25, beta=.3, gamma=.15,
                                switches=[.5, .4, .3, .2], noise_level=1)
    scv.tl.velocity(adata, mode='steady_state', vkey='steady_state_velocity', use_raw=True)
    dentategyrus = scv.datasets.dentategyrus()
    assert(adata is not None)
    assert(dentategyrus is not None)

# if recover dynamic runs, then it should not be None
def test_recover_dynamics():
    adata = scv.datasets.simulation(n_obs=500, t_max=25, beta=.3, gamma=.15,
                                switches=[.5, .4, .3, .2], noise_level=1)
    dm = scv.tl.recover_dynamics(adata, var_names=basis, max_iter=2, use_raw=True)
    assert (dm is not None)

# if recover dynamic runs, then it should not be None
def test_recover_dynamics():
    adata = scv.datasets.simulation(n_obs=500, t_max=25, beta=.3, gamma=.15,
                                switches=[.5, .4, .3, .2], noise_level=1)
    basis = adata.var_names[1]
    dm = scv.tl.recover_dynamics(adata, var_names=basis, max_iter=2, use_raw=True)
    assert(dm is not None)

# recovered from the ttils library should not be None. If it is, Figure 2 will not run!
#def test_utils_recover_dynamic():
    #dentategyrus = scv.datasets.dentategyrus()
    #recovered = utils.recover_dynamics(dentategyrus, ['Tmsb10', 'Ppp3ca', 'Hn1', 'Dlg2'])
    #assert(recovered is not None)
