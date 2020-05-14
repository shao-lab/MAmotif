"""
mamotif.region
--------------
Genomic regions use in MAmotif.
"""

import logging

import numpy as np
from motifscan.region import load_motifscan_regions

logger = logging.getLogger(__name__)


class MamotifRegion:
    """Class for a MAmotif genomic region.

    Parameters
    ----------
    chrom : str
        The chromosome name of the region.
    start : int
        The start coordinate of the region.
    end : int
        The end coordinate of the region.
    m_value: float, optional
        The m_value of the region.
    n_sites: list of int
        The motif sites numbers of the region.

    Attributes
    ----------
    chrom : str
        The chromosome name of the region.
    start : int
        The start coordinate of the region.
    end : int
        The end coordinate of the region.
    m_value : float or None
        The m_value of the region or None if not specified.
    has_motif : list of bool
        The target site indicators for motifs.

    Notes
    -----
    The coordinates are 0-based, which means the region range is [start, end).
    """

    def __init__(self, chrom, start, end, n_sites, m_value=None):
        self.chrom = chrom
        self.start = int(start)
        self.end = int(end)
        self.m_value = m_value
        self.has_motif = np.asarray(n_sites) > 0

    def match_manorm(self, manorm_regions):
        for region in manorm_regions:
            if (region.chrom == self.chrom) and (
                    region.start == self.start) and (region.end == self.end):
                self.m_value = region.score
                break
        if self.m_value is None:
            raise ValueError(f"no matched MAnorm region found for: {self!r}")

    def __repr__(self):
        return f"GenomicRegion({self.chrom}:{self.start}-{self.end})"


def load_mamotif_regions(f_manorm, f_motifscan):
    logger.info("Loading MAnorm result")
    manorm_regions = load_motifscan_regions(f_manorm, 'manorm')
    logger.info("Loading MotifScan result")
    logger.info(f"Loading genomic regions from {f_motifscan} [motifscan]")
    regions = []
    with open(f_motifscan, 'r') as fin:
        line = fin.readline()
        header = line.strip().split('\t')
        if header[:3] != ['chr', 'start', 'end']:
            raise ValueError(
                "not a valid MotifScan motif_sites_number.xls file")
        motifs = header[3:]
        for line in fin:
            fields = line.strip().split('\t')
            chrom = fields[0]
            start = int(fields[1]) - 1
            end = int(fields[2])
            n_sites = list(map(int, fields[3:]))
            region = MamotifRegion(chrom=chrom, start=start, end=end,
                                   n_sites=n_sites)
            regions.append(region)
    logger.info(f"Loaded {len(regions)} genomic regions")

    logger.info("Matching MAnorm and MotifScan results")
    if len(manorm_regions) != len(regions):
        logger.warning("the number of genomic regions are unmatched!")
    # group manorm regions by chrom to match with motifscan
    manorm_regions_by_chrom = {}
    for region in manorm_regions:
        manorm_regions_by_chrom.setdefault(region.chrom, [])
        manorm_regions_by_chrom[region.chrom].append(region)

    # find the matched manorm region and set the M value
    for region in regions:
        region.match_manorm(manorm_regions_by_chrom[region.chrom])
    return motifs, regions
