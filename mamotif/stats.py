"""
mamotif.stats
-------------

Statistical functions used in MAmotif.
"""

import numpy as np
from scipy import stats


def mamotif_t_test(m_values_pos, m_values_neg):
    try:
        t_stat, p_value = stats.ttest_ind(m_values_pos, m_values_neg,
                                          equal_var=False)
        if t_stat < 0:
            p_right = 1 - p_value / 2
        else:
            p_right = p_value / 2
        return t_stat, p_right
    except:
        return np.nan, np.nan


def mamotif_ranksum_test(m_values_pos, m_values_neg):
    try:
        z_stat, p_value = stats.ranksums(m_values_pos, m_values_neg)
        if z_stat < 0:
            p_right = 1 - p_value / 2
        else:
            p_right = p_value / 2
        return z_stat, p_right
    except:
        return np.nan, np.nan


def adjust_p_values(p_values, correction='benjamin'):
    n = len(p_values)
    adjusted_p_values = []
    if correction == 'benjamin':
        order = np.argsort(p_values)
        ranks = np.empty_like(order)
        ranks[order] = np.arange(1, n + 1)
        for p_value, rank in zip(p_values, ranks):
            if np.isnan(p_value):
                adjusted_p_values.append(np.nan)
            else:
                adjusted_p_values.append(min(1, p_value * n / rank))
    elif correction == 'bonferroni':
        for p_value in p_values:
            if np.isnan(p_value):
                adjusted_p_values.append(np.nan)
            else:
                adjusted_p_values.append(min(1, p_value * n))
    else:
        raise ValueError(f"invalid correction type: {correction}")
    return adjusted_p_values
