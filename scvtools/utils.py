import scvelo as scv
import matplotlib.pyplot as plt
import numpy as np

def recover_dynamics(data, var_names, max_iter = 100):
    """This function recovers the full splicing kinetics of specified genes.

    The return type must be of scvelo.tools.dynamical_model.DynamicsRecovery

    Parameters
    ----------
    data: :obj: anndata._core.anndata.AnnData
        The data of the specified genes.
    var_names
        The variables names that are accounted for.
    max_iter
        Default = 100
        The max number of iterations
       
    Returns
    -------
    scvelo.tools.dynamical_model.DynamicsRecovery
        recovering dynamics (using 1/8 cores)
            finished (0:00:00)

        outputs model fit of gene: 1

    """
    return scv.tl.recover_dynamics(data, var_names=var_names, max_iter=max_iter)

def plot_genes(data, basis):
    """This function plots a specified gene by passing in the data.

    The return type is None,  a figure will be outputted.
    This function visualizes a sepecified gene given the basis

    Parameters
    ----------
    data: :obj: anndata._core.anndata.AnnData
        The data of the specified genes.
    basis: list
        The variables names that are accounted for in the plot.
       
    Returns
    -------
        None
       
    """
    scv.pl.scatter(data, basis=basis, vkey='dynamics', fontsize=18)

def plot_contour(data, var_name, max_iter, x, y, c, n):
    """This function plots a specified gene by passing in the data.

    The return type is None, a figure will be outputted.
    This function visualizes the contour of a specified gene given the variable name

    Parameters
    ----------
    data: :obj: anndata._core.anndata.AnnData
        The data of the specified genes.
    var_name: str
        The variable name that is accounted for in the contour plot.
    max_iter: int
        The max number of iterations
    x: double
        The x_sight
    y: double
        The y_sight
    c: int
        The contour levels
    n: int
        The number of layers
    Returns
    -------
        None
       
    """
    dm = recover_dynamics(data, var_name, max_iter=max_iter)
    dm.plot_profile_contour(x_sight=x, y_sight=y, contour_levels=c, num=n)

def plot_validation_likelihood(gene):
    """This function plots the validation likelihood of specified genes.

    The return type is None, a figure will be outputted.

    Parameters
    ----------
    gene: :obj: anndata._core.anndata.AnnData
        The data of the specified genes

    Returns
    -------
        None

        Outputs a plot

    """
    top_genes = gene.var.sort_values('fit_likelihood', ascending=False).index
    rnd_genes = np.random.choice(gene.var_names, size=gene.n_vars)
    full_graph = gene.uns['dynamical_velocity_graph'].A
    adj = full_graph / full_graph > 0
    n_genes = np.array([10, 30, 300, 500, 800, 1000])
    rhos = []
    rhos_med = []
    for n in n_genes:
        if n == 0:
            rhos_med.extend([0])
        else:
            idx = gene.var_names.isin(top_genes[:n])

            vgraph = scv.VelocityGraph(gene[:, idx], vkey='dynamical_velocity')
            vgraph.compute_cosines()

            rho = scv.utils.vcorrcoef(full_graph * adj, vgraph.graph.A * adj)

            rhos.extend([rho])
            rhos_med.extend([np.nanmedian(rho)])

    rhos = np.array(rhos)
    rhos_med = np.array(rhos_med)

    ax = scv.pl.scatter(x=n_genes, y=rhos_med, size=1000, show=False, ylim=[.75, 1.04], fontsize=24)
    ax.plot(n_genes, rhos_med, color='grey', linewidth=5)
    #ax.set_yticklabels([.8, .9, .95, 1])
