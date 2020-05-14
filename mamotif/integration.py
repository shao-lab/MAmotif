"""
mamotif.integration
-------------------

The final integration module of MAmotif.
"""

import logging
import os

import numpy as np
from motifscan.genome import Genome
from motifscan.region.utils import subset_by_location

from mamotif.io import write_mamotif_results
from mamotif.region import load_mamotif_regions
from mamotif.stats import mamotif_t_test, mamotif_ranksum_test, adjust_p_values

logger = logging.getLogger(__name__)


class MAmotifResult:
    def __init__(self, motif, n_pos, mean_pos, std_pos, n_neg, mean_neg,
                 std_neg, t_stat, t_pval, t_padj, r_stat, r_pval, r_padj,
                 padj):
        self.motif = motif
        self.n_pos = n_pos
        self.mean_pos = mean_pos
        self.std_pos = std_pos
        self.n_neg = n_neg
        self.mean_neg = mean_neg
        self.std_neg = std_neg
        self.t_stat = t_stat
        self.t_pval = t_pval
        self.t_padj = t_padj
        self.r_stat = r_stat
        self.r_pval = r_pval
        self.r_padj = r_padj
        self.padj = padj


def mamotif_test(motifs, regions, negative=False, correction='benjamin'):
    results = []
    for idx, motif in enumerate(motifs):
        m_values_pos = []
        m_values_neg = []
        for region in regions:
            if region.has_motif[idx]:
                m_values_pos.append(region.m_value)
            else:
                m_values_neg.append(region.m_value)
        m_values_pos = np.asarray(m_values_pos)
        m_values_neg = np.asarray(m_values_neg)
        if negative:  # convert M to -M for sample B, log2(A/B)-> log2(B/A)
            m_values_pos = -m_values_pos
            m_values_neg = -m_values_neg
        t_test = mamotif_t_test(m_values_pos, m_values_neg)
        r_test = mamotif_ranksum_test(m_values_pos, m_values_neg)
        result = MAmotifResult(
            motif=motif, n_pos=len(m_values_pos), mean_pos=m_values_pos.mean(),
            std_pos=m_values_pos.std(), n_neg=len(m_values_neg),
            mean_neg=m_values_neg.mean(), std_neg=m_values_neg.std(),
            t_stat=t_test[0], t_pval=t_test[1], t_padj=None,
            r_stat=r_test[0], r_pval=r_test[1], r_padj=None,
            padj=None)
        results.append(result)

    # multiple testing correction
    p_values_t = [result.t_pval for result in results]
    p_values_r = [result.r_pval for result in results]
    adjusted_p_values_t = adjust_p_values(p_values_t, correction=correction)
    adjusted_p_values_r = adjust_p_values(p_values_r, correction=correction)
    for idx, result in enumerate(results):
        result.t_padj = adjusted_p_values_t[idx]
        result.r_padj = adjusted_p_values_r[idx]
        if np.isnan(result.t_padj):
            result.padj = result.r_padj
        elif np.isnan(result.r_padj):
            result.padj = result.t_padj
        else:
            result.padj = max(result.t_padj, result.r_padj)
    return results


def run_integration(f_manorm, f_motifscan, negative=False, genome=None,
                    split=False, upstream=4000, downstream=2000,
                    correction='benjamin', output_dir=None):
    motifs, regions = load_mamotif_regions(f_manorm, f_motifscan)
    sample_name = os.path.basename(f_manorm).replace('_MAvalues.xls', '')
    output_dir = os.path.abspath(output_dir or os.getcwd())
    if not os.path.isdir(output_dir):
        os.makedirs(output_dir)

    results = mamotif_test(motifs=motifs, regions=regions,
                           negative=negative, correction=correction)
    path = os.path.join(output_dir, sample_name + '_MAmotif_output.xls')
    write_mamotif_results(path=path, results=results, correction=correction)

    if split:
        logger.info("Split into promoter/distal regions")
        genome = Genome(genome)
        regions_promoter = subset_by_location(
            regions=regions, genes=genome.genes, location='promoter',
            upstream=upstream, downstream=downstream)
        regions_distal = subset_by_location(
            regions=regions, genes=genome.genes, location='distal',
            upstream=upstream, downstream=downstream)

        logger.info("Performing MAmotif on promoter regions")
        results = mamotif_test(motifs=motifs, regions=regions_promoter,
                               negative=negative, correction=correction)
        path = os.path.join(output_dir,
                            sample_name + '_promoter_MAmotif_output.xls')
        write_mamotif_results(path=path, results=results,
                              correction=correction)

        logger.info("Performing MAmotif on distal regions")
        results = mamotif_test(motifs=motifs, regions=regions_distal,
                               negative=negative, correction=correction)
        path = os.path.join(output_dir,
                            sample_name + '_distal_MAmotif_output.xls')
        write_mamotif_results(path=path, results=results,
                              correction=correction)
