# -*- coding: utf-8 -*-
"""Pure-Python statistics library (no numpy/scipy available on this platform).
Implements: t/Mann-Whitney/chi2/Fisher tests, logistic regression (IRLS),
LASSO logistic (coordinate descent + CV), AUC + bootstrap.
"""
import math, random

# ---------- special functions ----------
def _betacf(a, b, x):
    MAXIT, EPS, FPMIN = 200, 3e-12, 1e-300
    qab, qap, qam = a + b, a + 1.0, a - 1.0
    c, d = 1.0, 1.0 - qab * x / qap
    if abs(d) < FPMIN: d = FPMIN
    d = 1.0 / d
    h = d
    for m in range(1, MAXIT + 1):
        m2 = 2 * m
        aa = m * (b - m) * x / ((qam + m2) * (a + m2))
        d = 1.0 + aa * d
        if abs(d) < FPMIN: d = FPMIN
        c = 1.0 + aa / c
        if abs(c) < FPMIN: c = FPMIN
        d = 1.0 / d
        h *= d * c
        aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2))
        d = 1.0 + aa * d
        if abs(d) < FPMIN: d = FPMIN
        c = 1.0 + aa / c
        if abs(c) < FPMIN: c = FPMIN
        d = 1.0 / d
        de = d * c
        h *= de
        if abs(de - 1.0) < EPS: break
    return h

def betai(a, b, x):
    """Regularized incomplete beta I_x(a,b)."""
    if x <= 0.0: return 0.0
    if x >= 1.0: return 1.0
    lbeta = math.lgamma(a + b) - math.lgamma(a) - math.lgamma(b)
    bt = math.exp(lbeta + a * math.log(x) + b * math.log(1.0 - x))
    if x < (a + 1.0) / (a + b + 2.0):
        return bt * _betacf(a, b, x) / a
    return 1.0 - bt * _betacf(b, a, 1.0 - x) / b

def t_cdf(t, df):
    """Student t CDF."""
    x = df / (df + t * t)
    ib = betai(df / 2.0, 0.5, x)
    return 1.0 - 0.5 * ib if t >= 0 else 0.5 * ib

def norm_cdf(z):
    return 0.5 * (1.0 + math.erf(z / math.sqrt(2.0)))

def chi2_sf(x, df):
    """Chi-square survival function via gamma relation: P = Q(df/2, x/2) = 1 - I_{x/2}(df/2, inf)... use regularized gamma."""
    # Q(a,x) for gamma; use series/continued fraction (Numerical Recipes gammq)
    return _gammq(df / 2.0, x / 2.0)

def _gammq(a, x):
    if x < 0 or a <= 0: return float('nan')
    if x == 0: return 1.0
    if x < a + 1.0:  # series for P, return 1-P
        ap, s, de = a, 1.0 / a, 1.0 / a
        for _ in range(1000):
            ap += 1.0
            de *= x / ap
            s += de
            if abs(de) < abs(s) * 3e-12: break
        P = s * math.exp(-x + a * math.log(x) - math.lgamma(a))
        return 1.0 - P
    # continued fraction for Q
    b, c, d = x + 1.0 - a, 1e300, 1.0 / (x + 1.0 - a)
    h = d
    for i in range(1, 1000):
        an = -i * (i - a)
        b += 2.0
        d = an * d + b
        if abs(d) < 1e-300: d = 1e-300
        c = b + an / c
        if abs(c) < 1e-300: c = 1e-300
        d = 1.0 / d
        de = d * c
        h *= de
        if abs(de - 1.0) < 3e-12: break
    return math.exp(-x + a * math.log(x) - math.lgamma(a)) * h

# ---------- descriptive ----------
def mean(v): return sum(v) / len(v)
def sd(v):
    m = mean(v)
    return math.sqrt(sum((x - m) ** 2 for x in v) / (len(v) - 1)) if len(v) > 1 else 0.0
def median(v):
    s = sorted(v); n = len(s)
    return s[n // 2] if n % 2 else (s[n // 2 - 1] + s[n // 2]) / 2.0
def iqr(v):
    s = sorted(v); n = len(s)
    def q(p):
        idx = p * (n - 1)
        lo, hi = int(idx), min(int(idx) + 1, n - 1)
        return s[lo] + (s[hi] - s[lo]) * (idx - lo)
    return q(0.25), q(0.75)

# ---------- tests ----------
def ttest_welch(a, b):
    na, nb = len(a), len(b)
    ma, mb = mean(a), mean(b)
    va, vb = sd(a) ** 2, sd(b) ** 2
    se = math.sqrt(va / na + vb / nb)
    if se == 0: return float('nan'), float('nan')
    t = (ma - mb) / se
    df = (va / na + vb / nb) ** 2 / ((va / na) ** 2 / (na - 1) + (vb / nb) ** 2 / (nb - 1))
    p = 2.0 * (1.0 - t_cdf(abs(t), df))
    return t, p

def mann_whitney(a, b):
    """Normal approximation with tie & continuity correction, two-sided."""
    vals = [(x, 0) for x in a] + [(x, 1) for x in b]
    vals.sort(key=lambda t: t[0])
    ranks = [0.0] * len(vals)
    i, n = 0, len(vals)
    ties = []
    while i < n:
        j = i
        while j + 1 < n and vals[j + 1][0] == vals[i][0]: j += 1
        r = (i + j) / 2.0 + 1.0
        for k in range(i, j + 1): ranks[k] = r
        if j > i: ties.append(j - i + 1)
        i = j + 1
    R1 = sum(ranks[k] for k in range(n) if vals[k][1] == 0)
    n1, n2 = len(a), len(b)
    U1 = R1 - n1 * (n1 + 1) / 2.0
    mu = n1 * n2 / 2.0
    tie_corr = sum(t ** 3 - t for t in ties)
    sigma = math.sqrt(n1 * n2 / 12.0 * ((n1 + n2 + 1) - tie_corr / ((n1 + n2) * (n1 + n2 - 1))))
    if sigma == 0: return U1, float('nan')
    z = (U1 - mu - (0.5 if U1 > mu else -0.5 if U1 < mu else 0.0)) / sigma
    p = 2.0 * (1.0 - norm_cdf(abs(z)))
    return U1, p

def fisher_exact(a, b, c, d):
    """2x2 Fisher exact, two-sided (sum of tables with prob <= observed)."""
    def hyper_prob(a, b, c, d):
        n = a + b + c + d
        r1, c1 = a + b, a + c
        return math.exp(math.lgamma(r1 + 1) + math.lgamma(n - r1 + 1) + math.lgamma(c1 + 1) + math.lgamma(n - c1 + 1)
                        - math.lgamma(n + 1) - math.lgamma(a + 1) - math.lgamma(b + 1) - math.lgamma(c + 1) - math.lgamma(d + 1))
    p_obs = hyper_prob(a, b, c, d)
    r1, c1, n = a + b, a + c, a + b + c + d
    p = 0.0
    for x in range(max(0, r1 + c1 - n), min(r1, c1) + 1):
        px = hyper_prob(x, r1 - x, c1 - x, n - r1 - c1 + x)
        if px <= p_obs * (1 + 1e-9):
            p += px
    return min(p, 1.0)

def chi2_or_fisher(table2x2):
    (a, b), (c, d) = table2x2
    n = a + b + c + d
    exp = [(a + b) * (a + c) / n, (a + b) * (b + d) / n, (c + d) * (a + c) / n, (c + d) * (b + d) / n]
    if min(exp) < 5:
        return fisher_exact(a, b, c, d), 'Fisher'
    chi = (abs(a * d - b * c) - n / 2.0) ** 2 * n / ((a + b) * (c + d) * (a + c) * (b + d))
    return chi2_sf(chi, 1), 'chi2'

# ---------- logistic regression (IRLS) ----------
def logistic_fit(X, y, l2=1e-6, max_iter=100, tol=1e-8):
    """X: list of rows (list of float), intercept added internally. Returns (beta, se, converged)."""
    n = len(X); p = len(X[0]) + 1
    Xd = [[1.0] + row for row in X]
    beta = [0.0] * p
    converged = False
    for it in range(max_iter):
        eta = [sum(beta[j] * Xd[i][j] for j in range(p)) for i in range(n)]
        mu = [1.0 / (1.0 + math.exp(-max(-30, min(30, e)))) for e in eta]
        w = [max(m * (1 - m), 1e-9) for m in mu]
        # gradient and Hessian
        grad = [sum(Xd[i][j] * (y[i] - mu[i]) for i in range(n)) - l2 * beta[j] for j in range(p)]
        H = [[sum(Xd[i][j] * w[i] * Xd[i][k] for i in range(n)) + (l2 if j == k else 0.0)
              for k in range(p)] for j in range(p)]
        step = _solve(H, grad)
        if step is None: return beta, None, False
        beta = [beta[j] + step[j] for j in range(p)]
        if max(abs(s) for s in step) < tol:
            converged = True
            break
    # standard errors from inverse Hessian
    eta = [sum(beta[j] * Xd[i][j] for j in range(p)) for i in range(n)]
    mu = [1.0 / (1.0 + math.exp(-max(-30, min(30, e)))) for e in eta]
    w = [max(m * (1 - m), 1e-9) for m in mu]
    H = [[sum(Xd[i][j] * w[i] * Xd[i][k] for i in range(n)) + (l2 if j == k else 0.0)
          for k in range(p)] for j in range(p)]
    Hinv = _inverse(H)
    se = [math.sqrt(max(Hinv[j][j], 0.0)) for j in range(p)] if Hinv else None
    return beta, se, converged

def _solve(A, b):
    """Gaussian elimination with partial pivoting. Returns x or None."""
    n = len(A)
    M = [row[:] + [b[i]] for i, row in enumerate(A)]
    for col in range(n):
        piv = max(range(col, n), key=lambda r: abs(M[r][col]))
        if abs(M[piv][col]) < 1e-12: return None
        M[col], M[piv] = M[piv], M[col]
        for r in range(col + 1, n):
            f = M[r][col] / M[col][col]
            for c in range(col, n + 1):
                M[r][c] -= f * M[col][c]
    x = [0.0] * n
    for r in range(n - 1, -1, -1):
        x[r] = (M[r][n] - sum(M[r][c] * x[c] for c in range(r + 1, n))) / M[r][r]
    return x

def _inverse(A):
    n = len(A)
    M = [row[:] + [1.0 if i == j else 0.0 for j in range(n)] for i, row in enumerate(A)]
    for col in range(n):
        piv = max(range(col, n), key=lambda r: abs(M[r][col]))
        if abs(M[piv][col]) < 1e-12: return None
        M[col], M[piv] = M[piv], M[col]
        pv = M[col][col]
        M[col] = [v / pv for v in M[col]]
        for r in range(n):
            if r != col and abs(M[r][col]) > 1e-15:
                f = M[r][col]
                M[r] = [M[r][c] - f * M[col][c] for c in range(2 * n)]
    return [row[n:] for row in M]

def logit_predict(beta, X):
    out = []
    for row in X:
        e = beta[0] + sum(beta[j + 1] * row[j] for j in range(len(row)))
        e = max(-30, min(30, e))
        out.append(1.0 / (1.0 + math.exp(-e)))
    return out

# ---------- ROC / AUC ----------
def auc(scores, y):
    pos = [s for s, t in zip(scores, y) if t == 1]
    neg = [s for s, t in zip(scores, y) if t == 0]
    n1, n0 = len(pos), len(neg)
    if n1 == 0 or n0 == 0: return float('nan')
    cnt = 0.0
    for s in pos:
        for t in neg:
            cnt += 1.0 if s > t else 0.5 if s == t else 0.0
    return cnt / (n1 * n0)

def roc_curve(scores, y):
    pts = sorted(zip(scores, y), key=lambda t: -t[0])
    P = sum(y); N = len(y) - P
    curve = [(0.0, 0.0)]
    tp = fp = 0
    i = 0
    while i < len(pts):
        thr = pts[i][0]
        while i < len(pts) and pts[i][0] == thr:
            if pts[i][1] == 1: tp += 1
            else: fp += 1
            i += 1
        curve.append((fp / N, tp / P))
    curve.append((1.0, 1.0))
    return curve

def bootstrap_auc_ci(scores, y, B=2000, seed=42):
    rnd = random.Random(seed)
    n = len(y)
    stats = []
    for _ in range(B):
        idx = [rnd.randrange(n) for _ in range(n)]
        yy = [y[i] for i in idx]
        if sum(yy) == 0 or sum(yy) == n: continue
        stats.append(auc([scores[i] for i in idx], yy))
    stats.sort()
    lo = stats[int(0.025 * len(stats))]
    hi = stats[int(0.975 * len(stats))]
    return lo, hi

# ---------- LASSO logistic (coordinate descent on penalized log-likelihood) ----------
def soft_thresh(z, g):
    if z > g: return z - g
    if z < -g: return z + g
    return 0.0

def lasso_logistic_fit(X, y, lam, beta0=None, max_iter=200, tol=1e-6):
    """Cyclic coordinate descent, intercept unpenalized. X standardized externally."""
    n = len(X); p = len(X[0])
    beta = [0.0] * p
    b0 = math.log((sum(y) + 0.5) / (n - sum(y) + 0.5))
    if beta0 is not None:
        b0, beta = beta0[0], beta0[1][:]
    eta = [b0 + sum(X[i][j] * beta[j] for j in range(p)) for i in range(n)]
    for it in range(max_iter):
        max_chg = 0.0
        # intercept update (Newton step)
        mu = [1.0 / (1.0 + math.exp(-max(-30, min(30, e)))) for e in eta]
        g0 = sum(y[i] - mu[i] for i in range(n))
        h0 = sum(max(m * (1 - m), 1e-9) for m in mu)
        db0 = g0 / h0
        b0 += db0
        for i in range(n): eta[i] += db0
        max_chg = max(max_chg, abs(db0))
        for j in range(p):
            xj = [X[i][j] for i in range(n)]
            mu = [1.0 / (1.0 + math.exp(-max(-30, min(30, e)))) for e in eta]
            w = [max(m * (1 - m), 1e-9) for m in mu]
            # partial residual gradient for beta_j
            gj = sum(xj[i] * (y[i] - mu[i]) for i in range(n)) / n
            hj = sum(w[i] * xj[i] * xj[i] for i in range(n)) / n
            if hj < 1e-12: continue
            bj_old = beta[j]
            bj_new = soft_thresh(bj_old + gj / hj, lam / hj)
            db = bj_new - bj_old
            if db != 0.0:
                beta[j] = bj_new
                for i in range(n): eta[i] += db * xj[i]
                max_chg = max(max_chg, abs(db))
        if max_chg < tol: break
    return b0, beta

def lasso_cv(X, y, nfolds=10, nlam=100, lam_min_ratio=0.01, seed=42):
    """Single-pass 10-fold CV over lambda path (deviance). Returns lambda path, cv mean, fold devs, lambda.min/1se."""
    n = len(X); p = len(X[0])
    ybar = sum(y) / n
    lam_max = max(abs(sum(X[i][j] * (y[i] - ybar) for i in range(n)) / n) for j in range(p))
    lam_min = lam_max * lam_min_ratio
    lams = [lam_max * (lam_min / lam_max) ** (k / (nlam - 1)) for k in range(nlam)]
    rnd = random.Random(seed)
    pos = [i for i in range(n) if y[i] == 1]
    neg = [i for i in range(n) if y[i] == 0]
    rnd.shuffle(pos); rnd.shuffle(neg)
    folds = [[] for _ in range(nfolds)]
    for k, i in enumerate(pos): folds[k % nfolds].append(i)
    for k, i in enumerate(neg): folds[k % nfolds].append(i)
    fold_devs = [[] for _ in range(nlam)]
    for f in range(nfolds):
        test = set(folds[f])
        Xtr = [X[i] for i in range(n) if i not in test]
        ytr = [y[i] for i in range(n) if i not in test]
        Xte = [X[i] for i in range(n) if i in test]
        yte = [y[i] for i in range(n) if i in test]
        mu = [mean([r[j] for r in Xtr]) for j in range(p)]
        sdd = [sd([r[j] for r in Xtr]) or 1.0 for j in range(p)]
        Xtr_s = [[(r[j] - mu[j]) / sdd[j] for j in range(p)] for r in Xtr]
        Xte_s = [[(r[j] - mu[j]) / sdd[j] for j in range(p)] for r in Xte]
        warm = None
        for li, lam in enumerate(lams):
            b0, beta = lasso_logistic_fit(Xtr_s, ytr, lam, beta0=warm, max_iter=100, tol=1e-5)
            warm = (b0, beta)
            dev = 0.0
            for i in range(len(Xte)):
                e = b0 + sum(beta[j] * Xte_s[i][j] for j in range(p))
                e = max(-30, min(30, e))
                pr = 1.0 / (1.0 + math.exp(-e))
                pr = min(max(pr, 1e-12), 1 - 1e-12)
                dev += -2.0 * (yte[i] * math.log(pr) + (1 - yte[i]) * math.log(1 - pr))
            fold_devs[li].append(dev / len(Xte))
    cv_mean = [mean(fold_devs[li]) for li in range(nlam)]
    i_min = min(range(nlam), key=lambda i: cv_mean[i])
    se_min = sd(fold_devs[i_min]) / math.sqrt(nfolds)
    thr = cv_mean[i_min] + se_min
    cands = [i for i in range(nlam) if cv_mean[i] <= thr]
    i_1se = min(cands)
    return {'lams': lams, 'cv_mean': cv_mean, 'i_min': i_min, 'i_1se': i_1se,
            'lam_min': lams[i_min], 'lam_1se': lams[i_1se]}


def lasso_select(X, y, nfolds=10, nlam=50, lam_min_ratio=0.01, seed=42):
    """One repeat of CV-based selection; returns (set of selected idx at lambda.min, set at lambda.1se)."""
    p = len(X[0])
    cv = lasso_cv(X, y, nfolds=nfolds, nlam=nlam, lam_min_ratio=lam_min_ratio, seed=seed)
    sel_min, sel_1se = set(), set()
    b0, beta = lasso_logistic_fit(X, y, cv['lam_min'], max_iter=300, tol=1e-6)
    sel_min = {j for j in range(p) if abs(beta[j]) > 1e-8}
    b0, beta = lasso_logistic_fit(X, y, cv['lam_1se'], max_iter=300, tol=1e-6)
    sel_1se = {j for j in range(p) if abs(beta[j]) > 1e-8}
    return sel_min, sel_1se


# ---------- Firth logistic (Jeffreys prior) ----------
def firth_logistic_fit(X, y, max_iter=200, tol=1e-7):
    """Firth's penalized logistic regression. Returns (beta, se, converged)."""
    n = len(X); p = len(X[0]) + 1
    Xd = [[1.0] + row for row in X]
    beta = [0.0] * p
    converged = False
    for it in range(max_iter):
        eta = [sum(beta[j] * Xd[i][j] for j in range(p)) for i in range(n)]
        mu = [1.0 / (1.0 + math.exp(-max(-30, min(30, e)))) for e in eta]
        w = [max(m * (1 - m), 1e-9) for m in mu]
        XtWX = [[sum(Xd[i][j] * w[i] * Xd[i][k] for i in range(n)) for k in range(p)] for j in range(p)]
        XtWXinv = _inverse(XtWX)
        if XtWXinv is None: return beta, None, False
        # hat values
        h = []
        for i in range(n):
            s = 0.0
            for j in range(p):
                for k in range(p):
                    s += Xd[i][j] * XtWXinv[j][k] * Xd[i][k]
            h.append(w[i] * s)
        grad = [sum(Xd[i][j] * (y[i] - mu[i] + h[i] * (0.5 - mu[i])) for i in range(n)) for j in range(p)]
        step = [sum(XtWXinv[j][k] * grad[k] for k in range(p)) for j in range(p)]
        beta = [beta[j] + step[j] for j in range(p)]
        if max(abs(s) for s in step) < tol:
            converged = True
            break
    eta = [sum(beta[j] * Xd[i][j] for j in range(p)) for i in range(n)]
    mu = [1.0 / (1.0 + math.exp(-max(-30, min(30, e)))) for e in eta]
    w = [max(m * (1 - m), 1e-9) for m in mu]
    XtWX = [[sum(Xd[i][j] * w[i] * Xd[i][k] for i in range(n)) for k in range(p)] for j in range(p)]
    XtWXinv = _inverse(XtWX)
    se = [math.sqrt(max(XtWXinv[j][j], 0.0)) for j in range(p)] if XtWXinv else None
    return beta, se, converged
