import numba as nb
import numpy as np
import scipy.sparse as sps
import scipy.stats as sts
from tqdm.auto import tqdm

from decoupler._docs import docs
from decoupler._log import _log
from decoupler._Method import Method, MethodMeta


@nb.njit
def _score(rank: np.ndarray, miss: int, n_up: int, lgt: int)-> float:
    sum=miss*n_up
    sum=sum+rank.sum()
    s_min = lgt * (lgt + 1) / 2.0
    s_max = lgt * n_up
    score = 1.0 - (sum - s_min) / (s_max - s_min)
    return score


@nb.njit(parallel=True, cache=True)
def _uc(
    rank: np.ndarray,
    adj: np.ndarray,
    adjj: np.ndarray,
    n_up: int,
    nsrc: int,
    w_neg: float,
) -> np.ndarray:
    # Empty acts
    es = np.zeros(nsrc)
    # For each feature set
    for j in nb.prange(nsrc):
        # Extract feature set
        pos_len=(adj[:,j]>0).sum()
        neg_len=(adj[:,j]<0).sum()       
        adj_col = adjj[:, j]
        pos = rank[adj_col>0]
        neg = rank[adj_col<0]
        pos_miss=pos_len-len(pos)
        neg_miss=neg_len-len(neg)
        pos_score = _score(rank=pos, miss=pos_miss, n_up=n_up, lgt=pos_len) if len(pos) > 0 else 0.0
        neg_score = _score(rank=neg, miss=neg_miss, n_up=n_up, lgt=neg_len) if len(neg) > 0 else 0.0
        es[j] = max(0.0, pos_score - w_neg * neg_score)
    return es




@docs.dedent
def _func_ucell(
    mat: np.ndarray,
    adj: np.ndarray,
    n_up: int = 1500,
    verbose: bool = False,
    pvalue: bool = False,
    ties_method: str = "average",
    missing_genes: str = "skip",
    w_neg: float = 1.0,
) -> tuple[np.ndarray, None]:
    r"""
    Ucell



    %(notest)s

    %(params)s
    n_up
        Number of features to include in the UCell calculation.
        If ``None``, the top 1500 of features based on their magnitude are selected.

    %(returns)s

    Example
    -------
    .. code-block:: python

        import decoupler as dc

        adata, net = dc.ds.toy()
        dc.mt.ucell(adata, net, tmin=3)
    """
    nobs, nvar = mat.shape
    nsrc = adj.shape[1]
    m = f"ucell - calculating {nsrc} Uscore for {nvar} targets across {nobs} observations, categorizing features at rank={n_up}"
    _log(m, level="info", verbose=verbose)
    es = np.zeros(shape=(nobs, nsrc))
    for i in tqdm(range(mat.shape[0]), disable=not verbose):
        if isinstance(mat, sps.csr_matrix):
            row = mat[i].toarray()[0]
        else:
            row = mat[i]

        np.nan_to_num(row, copy=False)
        nz_idx = np.nonzero(row)[0]
        if len(nz_idx) == 0:
            continue
        row=row[nz_idx]
        ranks = sts.rankdata(-row, method=ties_method).astype(np.int32)
        keep_mask = ranks <= n_up
        kept_idx = nz_idx[keep_mask]
        kept_ranks = ranks[keep_mask]
        if len(kept_idx) > n_up:
            kept_idx = kept_idx[:n_up]
            kept_ranks = kept_ranks[:n_up]
        if len(kept_idx) == 0:
            continue
        adjj=adj[kept_idx]
        es[i] = _uc(rank=kept_ranks, adj=adj,adjj=adjj, n_up=n_up, nsrc=nsrc, w_neg=w_neg)
    return es, None


_ucell = MethodMeta(
    name="ucell",
    desc="UCell",
    func=_func_ucell,
    stype="categorical",
    adj=True,
    weight=True,
    test=False,
    limits=(0, 1),
    reference="https://doi.org/10.1093/bioinformatics/btag055",
)
ucell = Method(_method=_ucell)
