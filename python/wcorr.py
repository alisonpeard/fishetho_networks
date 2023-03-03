"""
Python adaptation of functions the R wCorr library.

Also contains functions to calculate other similarity metrics such as Cosine similarity and Hamming distance.

References:
-----------
..[1] Bailey P, Emad A (2023). _wCorr: Weighted Correlations_. R package
  version 1.9.6, <https://CRAN.R-project.org/package=wCorr>
"""
from pandas.core import (
    common as com,
    nanops
)
from pandas._libs import lib
import numpy as np


def wcorr(data, weights, method="Spearman", min_periods=1, numeric_only=True):
    """
    Calculate the weighted Spearman/Pearson correlation matrix of a dataframe.
    
    Parameters:
    -----------
    data : pandas.DataFrame
        DataFrame, where correlation will be calculated between each column.
    weights : pandas.DataFrame
        DataFrame (same dimensions as data) containing weights associated with each observation.
    method : string, {"Spearman", "Pearson"}
        A character string indicating which correlation coefficient is to be computed. 
    min_periods : int, optional
        Minimum number of observations required per pair of columns
        to have a valid result. Currently only available for Pearson
        and Spearman correlation.
    numeric_only : bool, default=True
        Include only `float`, `int` or `boolean` data.
        .. versionadded:: 1.5.0
        .. deprecated:: 1.5.0
            The default value of ``numeric_only`` will be ``False`` in a future
            version of pandas.
            
    Returns
    -------
    DataFrame
        Correlation matrix.
        
    """

    numeric_only_bool = com.resolve_numeric_only(numeric_only)
    if numeric_only_bool:
        data = data._get_numeric_data()
        weights = weights._get_numeric_data()

    if numeric_only is lib.no_default and len(data.columns) < len(self.columns):
        com.deprecate_numeric_only_default(type(self), "corr")

    cols = data.columns
    idx = cols.copy()
    mat = data.to_numpy(dtype=float, na_value=np.nan, copy=False)
    wmat = weights.to_numpy(dtype=float, na_value=np.nan, copy=False)

    if min_periods is None:
        min_periods = 1
    mat = mat.T
    wmat = wmat.T

    K = len(cols)
    correl = np.empty((K, K), dtype=float)
    mask = np.isfinite(mat)
    wmask = np.isfinite(wmat)

    for i, mats in enumerate(zip(mat, wmat)):
        ac, wac = mats[0], mats[1]
        for j, mats in enumerate(zip(mat, wmat)):
            bc, wbc = mats[0], mats[1]
            if i > j:
                continue

            wc = wac + wbc + 1
            valid = mask[i] & mask[j] & wmask[i] & wmask[j]
            if valid.sum() < min_periods:
                c = np.nan
            elif i == j:
                c = 1.0
            elif not valid.all():
                c = contCorr(ac[valid], bc[valid], wc[valid], method)
            else:
                c = contCorr(ac, bc, wc, method)

            correl[i, j] = c
            correl[j, i] = c
            
    return data._constructor(correl, index=idx, columns=cols)

    
def wrank(x, w=None):
    """
    Calculate the weighted rank of a vector.
    
    Python version of the private R function wrank from the  wCorr library. Access in R as wCorr:::wrank.
    
    References:
    -----------
    ..[1] Bailey P, Emad A (2023). _wCorr: Weighted Correlations_. R package
      version 1.9.6, <https://CRAN.R-project.org/package=wCorr>
    """
    if w is None:
        w = np.ones(len(x))
    else:
        w = np.array(w)
        
    # push that rank to all final tied units
    ord = np.argsort(x)
    rord = np.argsort(ord)
    xp = x[ord]  # x, permuted
    wp = w[ord]  # weights, permuted
    rnk = np.empty(len(x))
    
    # setup first iteration
    t1 = 0  # total weight of lower ranked elements
    i = 0   # index
    t2 = 0  # total weight of tied elements (including self)
    n = 0   # number of tied elements
    
    while i < len(x) - 1:
        t2 += wp[i]  # tied weight increases by this unit
        n += 1
        
        if xp[i+1] != xp[i]:   # the next one is not a tie
            # find the rank of all tied elements
            rnki = t1 + (1 + (t2 - 1)/2)
            
            # push that rank to all tied units
            for ii in range(n+1):  # (1, n+1)
                rnk[i-ii+1] = rnki
                
            # reset for next iteration
            t1 += t2
            t2 = 0
            n = 0
        
        i += 1
    
    # final row
    t2 += wp[i]
    rnki = t1 + (1 + (t2 - 1) / 2)  # final rank
    
    # push that rank to all final tied units
    for ii in range(n+1):
        rnk[i-ii] = rnki
        
    rnk = rnk[rord]
    
    return list(rnk)



def contCorr(x, y, w, method=['Pearson', 'Spearman']):
    """
    Calculate the weighted Pearson or Spearman rank of a vector.
    
    Python version of the private R function wrank from wCorr library. Access in R as wCorr:::contCorr.
    
    References:
    ----------
    ..[1] Bailey P, Emad A (2023). _wCorr: Weighted Correlations_. R package
      version 1.9.6, <https://CRAN.R-project.org/package=wCorr>
    """

    x = np.array(x, dtype=np.float64)
    y = np.array(y, dtype=np.float64)
    w = np.array(w, dtype=np.float64)

    pm = np.where(np.char.lower(np.array(method)) == np.char.lower(np.array(['Pearson', 'Spearman'])))
    if pm[0][0] == 1:
        # Spearman
        x = wrank(x, w)
        y = wrank(y, w)
        
    xb = np.sum(w * x) / np.sum(w)
    yb = np.sum(w * y) / np.sum(w)
    numerator = np.sum(w * (x - xb) * (y - yb))
    denom = np.sqrt(np.sum(w * (x - xb) ** 2) * np.sum(w * (y - yb) ** 2))
    return numerator / denom


def cosine_similarity(col1, col2):
    
    if (col1.sum() == 0) or (col2.sum() == 0):
        return np.nan
    
    return (col1 * col2).sum() / ( np.sqrt((col1 ** 2).sum()) * np.sqrt((col2 ** 2).sum()))


def hamming_dist(col1, col2):
    return np.sqrt(np.abs(col1 - col2).sum())


def nanhamming_dist(col1, col2):
    valid = np.isfinite(col1) & np.isfinite(col2)
    P = sum(valid)
    T = len(col1)

    return np.sqrt((T / P) * np.abs(col1[valid] - col2[valid]).sum())


def hamming_similarity(col1, col2):
    hdist = hamming_dist(col1, col2)
    return 1 / (1 + hdist)


def nanhamming_similarity(col1, col2):  
    hdist = nanhamming_dist(col1, col2)
    return 1 / (1 + hdist)


def test_wrank():
    # Identity
    x = np.array([1, 2, 3, 4, 5])
    assert (wrank(x) == np.array([1, 2, 3, 4, 5])).all()
    
    # General example
    x = np.array([2, 4, 6, 8, 10])
    assert (wrank(x) == np.array([1, 2, 3, 4, 5])).all()
    
    # Averages ties
    x = np.array([2, 4, 4, 8, 10])
    assert (wrank(x) == np.array([1, 2.5, 2.5, 4, 5])).all()
    
    # Weighted
    x = np.array([2, 4, 6, 8, 10])
    w = np.array([0.5, 1, 1, 1.5, 2])
    assert (wrank(x, w) == np.array([0.75, 1.5, 2.5, 3.75, 5.5])).all()
    
    # Disorderly + weighted
    x = np.array([3, 4, 6, 1, 2])
    w = np.array([0.5, 1, 1, 1.5, 2])
    assert(wrank(x, w) == np.array([4.25, 5., 6., 1.25, 3.])).all()


def test_contCorr():
    # Test with Pearson method
    x = [1, 2, 3, 4, 5]
    y = [2, 4, 6, 8, 10]
    w = [1, 1, 1, 1, 1]
    assert np.isclose(contCorr(x, y, w, method=['Pearson']), 1.0)

    # Test with Spearman method
    x = [1, 2, 3, 4, 5]
    y = [2, 4, 6, 8, 10]
    w = [1, 1, 1, 1, 1]
    assert np.isclose(contCorr(x, y, w, method=['Spearman']), 1.0)

    # Test with weighted data
    x = [1, 2, 3, 4, 5]
    y = [2, 4, 6, 8, 10]
    w = [0.5, 1, 1, 1.5, 2]
    assert np.isclose(contCorr(x, y, w, method=['Spearman']), 1.0)
    
    # Test with negatively correlated data
    x = [1, 2, 3, 4, 5]
    y = [3, 4, 6, 1, 2]
    w = [0.5, 1, 1, 1.5, 2]
    assert np.isclose(contCorr(x, y, w, method=['Pearson']), -0.5430768084288788)
    
    # Test with negatively correlated data
    x = np.array([1, 2, 3, 4, 5])
    y = np.array([3, 4, 6, 1, 2])
    w = np.array([0.5, 1, 1, 1.5, 2])
    assert np.isclose(contCorr(x, y, w, method=['Spearman']), -0.5555556)

    # Test with missing data
    x = [1, 2, 3, 4, np.nan]
    y = [2, 4, 6, 8, 10]
    w = [1, 1, 1, 1, 1]
    assert np.isnan(contCorr(x, y, w, method=['Pearson']))
    
    # Test with tied data
    x = np.array([1, 2, 3, 4, 5])
    y = np.array([3, 2, 1, 2, 4])
    w = np.array([0.5, 1, 1, 1.5, 2])
    assert np.isclose(contCorr(x, y, w, method=['Spearman']), 0.6516946235415335)
