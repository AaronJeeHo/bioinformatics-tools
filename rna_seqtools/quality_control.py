"""
Author: Aaron Ho
Python Version: 3.8
"""
import scanpy as sc
import numpy as np

from scipy import stats


def get_adata(matrix_path):
    """
    Convert 10x Feature-Barcode Matrix file to AnnData
    :param str matrix_path: path to 10x feature-barcode matrix file
    :return anndata: Anndata object
    """
    sc.settings.verbosity = 3
    adata = sc.read_10x_mtx(matrix_path, var_names='gene_symbols', cache=True)
    adata.var_names_make_unique()
    return adata


def find_cutoffs(data, dev):
    """
    Find minimum and maximum range of data to be within given
     Median Absolute Deviation threshold
    :param iterable data: data to find deviation ranges
    :param int dev: deviation value
    :return (int, int): minimum and maximum values within deviation
    """
    median = np.median(data)
    mad = stats.median_abs_deviation(data)

    min_range = ((-1 * dev) * mad) + median    # cannot be less than cutoff
    max_range = (dev * mad) + median

    if min_range < 0:
        return 0, max_range
    else:
        return min_range, max_range


def preprocess_adata(adata, min_genes, min_cells, min_counts):
    """
    Basic QC for anndata filtering minimum requirements
    :param anndata adata: anndata matrix
    :param int min_genes: Minimum genes expressed in cells
    :param int min_cells: Minimum cells genes expressed in
    :param int min_counts: Minimum counts genes expressed in
    """
    sc.pp.filter_cells(adata, min_genes=min_genes)
    sc.pp.filter_genes(adata, min_cells=min_cells)
    sc.pp.filter_genes(adata, min_cells=min_counts)
    print('Data preprocessed')


def filter_adata(adata):
    """
    Filter data with counts over 5% mitochondrial or
     3x larger than Median Absolute Deviation
    :param anndata adata: anndata Matrix
    :return anndata: Filtered Matrix
    """
    adata.var['mt'] = adata.var_names.str.startswith(('MT-', 'mt-'))

    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'],
                               percent_top=None, log1p=False, inplace=True)

    count_range = find_cutoffs(adata.obs['total_counts'], 3)
    max_counts = int(count_range[1])

    return adata[(adata.obs['pct_counts_mt'] < 5)
                 & (adata.obs['total_counts'] < max_counts), :]


def normalize_adata(adata, norm_factor):
    """
    Normalize counts in matrix
    :param anndata adata: Anndata Matrix
    :param int norm_factor: normalization factor
    """
    sc.pp.normalize_total(adata, target_sum=norm_factor)
    sc.pp.log1p(adata)


def identify_variable(adata):
    """
    Identify and plot highly variable genes
    :param adata: Anndata Matrix
    """
    sc.pp.highly_variable_genes(adata,
                                min_mean=0.0125, max_mean=3, min_disp=0.5)
    sc.pl.highly_variable_genes(adata)
