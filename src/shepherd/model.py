"""The Bayesian test between error sequence and true barcode (Section 2.2.2)."""

import math


def log_binom_pmf(k, n, p):
    """Log of the binomial probability mass function, ln C(n, k) p^k (1 - p)^(n - k).

    The same formula as scipy.stats.binom.logpmf, without SciPy's per-call
    overhead (about 60 times faster for scalar arguments).
    """
    log_comb = math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)
    return log_comb + k * math.log(p) + (n - k) * math.log1p(-p)


def get_log_K(f_n, f_c, p_no_err, d, l, total_err_rate, logdenom):
    """Log Bayes factor ln K (Eq. 3) of a sequence with count f_n at distance d
    from a putative barcode with count f_c."""

    n_hat = max(int(f_c / p_no_err), f_c + f_n)
    p_est = (total_err_rate / 3) ** d * (1 - total_err_rate) ** (l - d)

    return log_binom_pmf(f_n, n_hat, p_est) + math.log(p_est) + logdenom
