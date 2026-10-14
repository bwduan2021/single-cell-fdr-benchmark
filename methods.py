"""Gene-level statistics and FDR procedures for the single-cell benchmark.

Inputs are cell-by-gene matrices. Legacy keys deseq2/mast denote
simplified scores, not calls to the corresponding Bioconductor packages.
The working e-BH rule requires appropriate null-score calibration.
"""


import hashlib
import numpy as np
import scipy.stats as stats

SEED = 2025
ALPHA = 0.05


def stable_seed(*parts):
    """Return a process-independent seed; Python's hash() is intentionally avoided."""
    payload = "|".join(map(str, parts)).encode("utf-8")
    return SEED + int.from_bytes(hashlib.sha256(payload).digest()[:4], "little")


def stat_t(XA, XB):
    nA, nB = XA.shape[0], XB.shape[0]
    vA = XA.var(0, ddof=1) + 1e-10
    vB = XB.var(0, ddof=1) + 1e-10
    se = np.sqrt(vA / nA + vB / nB)
    z  = (XA.mean(0) - XB.mean(0)) / se
    df = (vA/nA + vB/nB)**2 / ((vA/nA)**2/(nA-1) + (vB/nB)**2/(nB-1) + 1e-20)
    pv = 2.0 * stats.t.sf(np.abs(z), df=np.clip(df, 1, None))
    return pv, z


def stat_wilcoxon(XA, XB):
    """Wilcoxon rank-sum statistic with tie correction."""
    nA, nB = XA.shape[0], XB.shape[0]
    N      = nA + nB
    combined = np.vstack([XA, XB])
    ranks    = np.apply_along_axis(stats.rankdata, 0, combined)
    W        = ranks[:nA].sum(0)
    mu_W     = nA * (N + 1) / 2.0
    n_genes  = XA.shape[1]
    z_out    = np.empty(n_genes)
    for g in range(n_genes):
        c = combined[:, g]
        _, cnt = np.unique(c, return_counts=True)
        tc       = np.sum(cnt**3 - cnt) / 12.0
        s2       = nA * nB * ((N+1)/12.0 - tc/(N*(N-1)+1e-20))
        z_out[g] = (W[g] - mu_W) / max(np.sqrt(s2), 1e-10)
    pv = 2.0 * stats.norm.sf(np.abs(z_out))
    return pv, z_out


def stat_deseq2(XA, XB):
    """Negative-binomial score with moment-estimated dispersion."""
    nA, nB = XA.shape[0], XB.shape[0]
    eps    = 0.5
    muA    = XA.mean(0) + eps
    muB    = XB.mean(0) + eps
    mu_p   = np.vstack([XA, XB]).mean(0) + eps
    va_p   = np.vstack([XA, XB]).var(0, ddof=1) + 1e-10
    phi    = np.clip(mu_p**2 / np.maximum(va_p - mu_p, 1e-3), 0.1, 100.0)
    seA    = np.sqrt((1 + muA/phi) / (nA * muA))
    seB    = np.sqrt((1 + muB/phi) / (nB * muB))
    z      = (np.log(muA) - np.log(muB)) / np.sqrt(seA**2 + seB**2 + 1e-20)
    pv     = 2.0 * stats.norm.sf(np.abs(z))
    return pv, z


def stat_mast(XA, XB):
    """Two-part score combining detection and nonzero intensity."""
    nA, nB = XA.shape[0], XB.shape[0]
    # Detection component: signed square root of the chi-square statistic.
    nzA = (XA > 0).sum(0).astype(float); zA = nA - nzA
    nzB = (XB > 0).sum(0).astype(float); zB = nB - nzB
    cont   = np.stack([np.stack([nzA, zA], -1),
                       np.stack([nzB, zB], -1)], axis=-2)
    rs     = cont.sum(-1, keepdims=True)
    cs     = cont.sum(-2, keepdims=True)
    ex     = rs * cs / (nA + nB + 1e-10)
    chi2_v = np.sum((cont - ex)**2 / (ex + 1e-10), axis=(-2, -1))
    z_disc = np.sqrt(chi2_v) * np.sign(nzA/nA - nzB/nB)
    # Intensity component: Welch-type statistic on nonzero log counts.
    logXA  = np.where(XA > 0, np.log(XA + 0.5), np.nan)
    logXB  = np.where(XB > 0, np.log(XB + 0.5), np.nan)
    cA_    = np.sum(~np.isnan(logXA), 0).astype(float)
    cB_    = np.sum(~np.isnan(logXB), 0).astype(float)
    vA_    = np.nanvar(logXA, 0, ddof=1) + 1e-10
    vB_    = np.nanvar(logXB, 0, ddof=1) + 1e-10
    se_    = np.sqrt(vA_/np.maximum(cA_, 2) + vB_/np.maximum(cB_, 2))
    z_cont = np.where((cA_ >= 2) & (cB_ >= 2),
                      (np.nanmean(logXA, 0) - np.nanmean(logXB, 0)) / se_, 0.0)
    z  = (z_disc + z_cont) / np.sqrt(2)
    pv = 2.0 * stats.norm.sf(np.abs(z))
    return pv, z


STAT_FUNC = {"t": stat_t, "wilcoxon": stat_wilcoxon,
             "deseq2": stat_deseq2, "mast": stat_mast}


def bh(pv, alpha=ALPHA):
    m   = len(pv)
    o   = np.argsort(pv)
    thr = np.arange(1, m + 1) * alpha / m
    bl  = pv[o] <= thr
    if not bl.any():
        return np.zeros(m, bool)
    k   = int(np.where(bl)[0][-1])
    rej = np.zeros(m, bool)
    rej[o[:k + 1]] = True
    return rej


def ebh(z_stat, c=2.0, alpha=ALPHA):
    """Working e-BH rule using the normal-score transform."""
    m        = len(z_stat)
    ev       = np.exp(-c**2 / 2) * np.cosh(c * np.abs(z_stat))
    order    = np.argsort(ev)[::-1]
    ev_sort  = ev[order]
    eligible = ev_sort >= m / (np.arange(1, m + 1) * alpha)
    if not eligible.any():
        return np.zeros(m, bool)
    k_star   = int(np.where(eligible)[0][-1]) + 1
    return ev >= ev_sort[k_star - 1]


def _td_scan(W, alpha=ALPHA):
    """Select the largest feasible rejection set using the plus-one rule."""
    best = np.zeros(len(W), bool)
    best_R = 0
    for t in np.sort(np.abs(W))[::-1]:
        R  = int(np.sum(W >= t))
        Dv = int(np.sum(W <= -t))
        if R == 0:
            continue
        if (1 + Dv) / R <= alpha and R > best_R:
            best   = W >= t
            best_R = R
    return best


def independently_permute_columns(Xall, n_a, rng):
    """Permute labels independently per gene while preserving group sizes."""
    n_total, n_genes = Xall.shape
    keys = rng.random((n_total, n_genes))
    order = np.argpartition(keys, n_a - 1, axis=0)
    columns = np.arange(n_genes)[None, :]
    return Xall[order[:n_a], columns], Xall[order[n_a:], columns]


def target_decoy_competition(target_score, decoy_score, rng):
    """Construct signed TD scores, breaking exact target--decoy ties at random."""
    target_score = np.abs(np.asarray(target_score, dtype=float))
    decoy_score = np.abs(np.asarray(decoy_score, dtype=float))
    delta = target_score - decoy_score
    direction = np.sign(delta)
    ties = delta == 0
    if np.any(ties):
        direction[ties] = rng.choice(np.array([-1.0, 1.0]), size=int(ties.sum()))
    return np.maximum(target_score, decoy_score) * direction


def ds_procedure(XA, XB, stat_key, rng, alpha=ALPHA):
    """Single-split mirror procedure with a plus-one threshold."""
    nA, nB = XA.shape[0], XB.shape[0]
    hA, hB = nA // 2, nB // 2
    pA     = rng.permutation(nA); pB = rng.permutation(nB)
    A1, A2 = XA[pA[:hA]], XA[pA[hA:hA*2]]
    B1, B2 = XB[pB[:hB]], XB[pB[hB:hB*2]]
    _, T1  = STAT_FUNC[stat_key](A1, B1)
    _, T2  = STAT_FUNC[stat_key](A2, B2)
    M      = np.sign(T1 * T2) * (np.abs(T1) + np.abs(T2))
    return _td_scan(M, alpha)


def td_procedure(XA, XB, stat_key, rng, alpha=ALPHA):
    """Target-decoy procedure with matched scores and randomized ties."""
    nA = XA.shape[0]
    Xall   = np.vstack([XA, XB])
    _, S_t = STAT_FUNC[stat_key](XA, XB)
    S_t    = np.abs(S_t)
    XA_d, XB_d = independently_permute_columns(Xall, nA, rng)
    _, S_d = STAT_FUNC[stat_key](XA_d, XB_d)
    S_d = np.abs(S_d)
    W = target_decoy_competition(S_t, S_d, rng)
    return _td_scan(W, alpha)


def run_one(XA, XB, stat_key, ctrl_key, rng):
    if ctrl_key == "DS":
        return ds_procedure(XA, XB, stat_key, rng)
    if ctrl_key == "TD":
        return td_procedure(XA, XB, stat_key, rng)
    pv, z = STAT_FUNC[stat_key](XA, XB)
    if ctrl_key == "BH":
        return bh(pv)
    return ebh(z)
