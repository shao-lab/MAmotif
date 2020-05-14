"""
mamotif.io
----------
I/O module of MAmotif.
"""

import logging

logger = logging.getLogger(__name__)


def write_mamotif_results(path, results, correction):
    logger.info(f"Saving MAmotif results to {path}")
    correction_str = correction.capitalize()
    columns = ["Motif Name", "Target Number", "Average of Target M values",
               "Std. of Target M values", "Non-target Number",
               "Average of Non-target M values", "Std. of Non-target M values",
               "T-test Statistic", "T-test P value (right-tailed)",
               f"T-test P value By {correction_str} correction",
               "RankSum-test Statistic", "RankSum-test P value (right-tailed)",
               f"RankSum-test P value By {correction_str} correction\t",
               "Maximal corrected P value\n"]
    header = '\t'.join(columns)
    results.sort(key=lambda x: x.padj)
    with open(path, 'w') as fout:
        fout.write(header)
        for result in results:
            fout.write(
                f"{result.motif}\t{result.n_pos}\t{result.mean_pos}\t"
                f"{result.std_pos}\t{result.n_neg}\t{result.mean_neg}\t"
                f"{result.std_neg}\t{result.t_stat}\t{result.t_pval}\t"
                f"{result.t_padj}\t{result.r_stat}\t{result.r_pval}\t"
                f"{result.r_padj}\t{result.padj}\n")
