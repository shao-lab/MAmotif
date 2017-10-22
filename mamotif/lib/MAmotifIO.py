import os
import sys
import numpy as np
import pandas as pd
import logging

from MAmotifPeaks import MAnormPeak, MotifScanPeak, MAnormPeakSet
from MAnormPeaksClassifier import MAnormPeaksClassifier, FeaturePromoter, Feature
from sequence import Gene


def coarse(num, keli):
    return int(num / keli) * keli


def read_MAnorm_peaks(file_path, neg=False):
    logging.info("Reading MAnorm results...")
    pk_num = 0
    handle = open(file_path)
    pks = []
    for info in handle:
        if not info.startswith('#'):
            cdt = info.strip().split()
            try:
                pk = MAnormPeak(cdt[0], int(cdt[1]), int(cdt[2]), summit=int(cdt[3]))
            except:
                try:
                    pk = MAnormPeak(cdt[0], int(cdt[1]), int(cdt[2]))
                except:
                    continue
            # pk.set_mvalue(float(cdt[4]))
            try:
                if neg:
                    pk.set_mvalue(-round(float(cdt[4]), 1))  # Keep a decimal point
                else:
                    pk.set_mvalue(round(float(cdt[4]), 1))
            except:
                pass
            # pk.set_mvalue(coarse(float(cdt[4]), 0.5))
            try:
                pk.set_avalue(float(cdt[5]))
            except:
                pass
            try:
                pk.set_pvalue(float(cdt[6]))
            except:
                pass
            pks.append(pk)
            pk_num += 1
    handle.close()
    return pks


def read_refgenes(gene_file):
    """
    read ref genes
    # -------------------------------------------------------------------------------------------
    # field	example	SQL type	info	description

    # bin	637	smallint(5) unsigned	range	Indexing field to speed chromosome range queries.
    # name	NM_021010	varchar(255)	values	Name of gene (usually transcript_id from GTF)
    # chrom	chr8	varchar(255)	values	Reference sequence chromosome or scaffold
    # strand	-	char(1)	values	+ or - for strand
    # txStart	6912828	int(10) unsigned	range	Transcription start position
    # txEnd	6914259	int(10) unsigned	range	Transcription end position
    # cdsStart	6912952	int(10) unsigned	range	Coding region start
    # cdsEnd	6914219	int(10) unsigned	range	Coding region end
    # exonCount	2	int(10) unsigned	range	Number of exons
    # exonStarts	6912828,6914047,	longblob	 	Exon start positions
    # exonEnds	6913065,6914259,	longblob	 	Exon end positions
    # score	0	int(11)	range	score
    # name2	DEFA5	varchar(255)	values	Alternate name (e.g. gene_id from GTF)
    # cdsStartStat	cmpl	enum('none', 'unk', 'incmpl', 'cmpl')	values	enum('none','unk','incmpl','cmpl')
    # cdsEndStat	cmpl	enum('none', 'unk', 'incmpl', 'cmpl')	values	enum('none','unk','incmpl','cmpl')
    # exonFrames	1,0,	longblob	 	Exon frame {0,1,2}, or -1 if no frame for exon
    # --------------------------------------------------------------------------------------------------------
    @param gene_file: gene info file, refgene download from ENCODE
    """
    gene_list = []
    handle = open(gene_file)
    for gene_info in handle:
        ctt = gene_info.split()
        gene = Gene(ctt[1], ctt[12], ctt[2], int(ctt[4]), int(ctt[5]), ctt[3])
        gene_list.append(gene)
    handle.close()
    return gene_list


def read_motifscan_result(motifscan_file_path):
    logging.info('Reading Motifscan result...')
    motifscan_peak_list = []
    if motifscan_file_path.endswith('.mat'):  # from Shao MotifScan result
        from scipy import io as sio
        mat = sio.loadmat(motifscan_file_path)
        mat = mat['peak']
        motifs = [motif[0] for motif in mat['motif_name'][0][0][0]]

        for i in xrange(len(mat['chr_id'][0][0])):
            chrm = mat['chr_id'][0][0][i][0][0]
            start = mat['bpstart'][0][0][i][0]
            end = mat['bpend'][0][0][i][0]
            summit = mat['bpsummit'][0][0][i][0]
            tarnum = mat['motif_tarnum'][0][0][i]
            pk = MotifScanPeak(chrm, start, end, summit)
            pk.set_motif_info(motifs, tarnum)
            motifscan_peak_list.append(pk)
        return motifscan_peak_list, motifs

    if motifscan_file_path.endswith('.csv'):
        motifscan = pd.read_csv(motifscan_file_path)
    else:
        motifscan = pd.read_pickle(motifscan_file_path)  # from Wang Jiawei MotifScan result
    motifs = get_motif_names(motifscan)
    i = 0
    row_num = motifscan.shape[0]
    for k, v in motifscan.iterrows():
        i += 1
        motifscan_peak = MotifScanPeak(v['chr'], int(v['start']), int(v['end']), int(v['summit']))
        target_number_list = [int(v[name + '.number']) for name in motifs]
        motifscan_peak.set_motif_info(motifs, target_number_list)
        motifscan_peak_list.append(motifscan_peak)
    return motifscan_peak_list, motifs


def get_motif_names(motifscan):
    motifs = []
    for c in motifscan.columns:
        if c.endswith('number'):
            motifs.append(c[:-7])
    return motifs


def match_manorm_with_motifscan(pk_file, motifscan_result, neg_mvalue=False):
    """
    read MAnorm peaks and motifscan result, than match the two result.
    """
    # read MAnorm peaks
    manorm_pks = read_MAnorm_peaks(pk_file, neg_mvalue)

    # read motifscan result
    motifscan_pks, motifs = read_motifscan_result(motifscan_result)
    matched_manorm_pks = []
    matched_tarnum_list = []
    for pk in manorm_pks:
        match = False
        for mp in motifscan_pks:
            if pk == mp:
                matched_manorm_pks.append(pk)
                matched_tarnum_list.append(mp.target_number)
                match = True
                break
        if not match:
            logging.warning('MAnorm and Motifscan are not match!')
            pk.prints()

    tarnum_dict = {}
    for i, motif in enumerate(motifs):
        tarnum_dict[motif] = np.array([e[i] for e in matched_tarnum_list])
    return matched_manorm_pks, tarnum_dict, motifs


def classify_MAnorm_pks_by_promoter(pk_fp, refgene_fp):
    # get peaks
    pk_file = os.path.split(pk_fp)[1]
    pklist = read_MAnorm_peaks(pk_fp)
    pkset = MAnormPeakSet(pklist)
    # get feature
    feature = FeaturePromoter(read_refgenes(refgene_fp))
    # make a classifier
    classifier = MAnormPeaksClassifier(pkset, feature)
    classifier.classify_by_feature()
    promoter_pkset = classifier.feature_yes
    distal_pkset = classifier.feature_no

    # saving classified result
    header = \
        '\t'.join(['#chr', 'start', 'end', 'summit', 'MAnorm_Mvalue', 'MAnorm_Avalue', 'pvalue\n'])
    promoter_file_name = pk_file.replace('peak', 'promoter_peak').replace('Peaks', 'PromoterPeaks')
    file_promoter = open(promoter_file_name, 'w')
    file_promoter.write(header)
    [file_promoter.write(pk.tostring()) for pk in promoter_pkset]
    file_promoter.close()
    distal_file_name = pk_file.replace('peak', 'distal_peak').replace('Peaks', 'DistalPeaks')
    file_distal = open(distal_file_name, 'w')
    file_distal.write(header)
    [file_distal.write(pk.tostring()) for pk in distal_pkset]
    file_distal.close()
    return promoter_file_name, distal_file_name


def _get_header(correction_type):
    if correction_type == 'benjamin':
        header = \
            'Motif Name\t' \
            'Target Number\tAverage of Target M-value\tDeviation of Target M-value\t' \
            'Non-target Number\tAverage of Non-target M-value\tDeviation of Non-target M-value\t' \
            'T-test Statistics\tT-test P-value(right-tail)\tT-test P-value By Benjamin correction\t' \
            'RanSum-test Statistics\tRankSum-test P-value(right-tail)\tRankSum-test P-value By Benjamin correction\t' \
            'Maximal P-value\n'
    elif correction_type == 'bonferroni':
        header = \
            'Motif Name\t' \
            'Target Number\tAverage of Target M-value\tDeviation of Target M-value\t' \
            'Non-target Number\tAverage of Non-target M-value\tDeviation of Non-target M-value\t' \
            'T-test Statistics\tT-test P-value(right-tail)\tT-test P-value By Bonferroni correction\t' \
            'RanSum-test Statistics\tRankSum-test P-value(right-tail)\tRankSum-test P-value By Bonferroni correction\t' \
            'Maximal P-value\n'
    else:
        header = \
            'Motif Name\t' \
            'Target Number\tAverage of Target M-value\tDeviation of Target M-value\t' \
            'Non-target Number\tAverage of Non-target M-value\tDeviation of Non-target M-value\t' \
            'T-test Statistics\tT-test P-value(right-tail)\tT-test P-value By correction\t' \
            'RanSum-test Statistics\tRankSum-test P-value(right-tail)\tRankSum-test P-value By correction\t' \
            'Maximal P-value\n'
    return header


def motif_classified_pks_ttest(
        pk_file_path, motifscan_result_path, negative=False, correction_type='benjamin'):
    """
    match peaks with motifscan result once. then output the test result of all jaspar motifs.
    """
    pk_list, tarnum_dict, motifs = \
        match_manorm_with_motifscan(pk_file_path, motifscan_result_path, negative)
    pk_array = np.array(pk_list)

    # classify pks and do t-test&ranksum-test for each motif
    lines = []
    t_stat, ttest_pvalue, r_stat, rtest_pvalue = [], [], [], []
    for moti in motifs:
        classifier = MAnormPeaksClassifier(MAnormPeakSet(), Feature())
        yes = pk_array[np.where(tarnum_dict[moti] > 0)[0]]
        no = pk_array[np.where(tarnum_dict[moti] == 0)[0]]
        pkset_yes = MAnormPeakSet()
        if yes.size > 0:
            pkset_yes.set_sequences(yes)
        pkset_no = MAnormPeakSet()
        if no.size > 0:
            pkset_no.set_sequences(no)
        classifier.set_feature(pkset_yes, pkset_no)
        line = '%s\t' % moti
        line += '%d\t%f\t%f\t' % (
            classifier.feature_yes.size, classifier.feature_yes.mean, classifier.feature_yes.std)
        line += '%d\t%f\t%f\t' % (
            classifier.feature_no.size, classifier.feature_no.mean, classifier.feature_no.std)
        lines.append(line)
        ttest = classifier.ttest_feature_classified_peaks()
        t_stat.append(ttest[0])
        ttest_pvalue.append(ttest[1])
        rtest = classifier.ranksum_feature_classified_peaks()
        # rtest = classifier.kstest_feature_classified_peaks()
        r_stat.append(rtest[0])
        rtest_pvalue.append(rtest[1])

    corrected_ttest_pvalue = correct_pvalues(ttest_pvalue, correction_type)
    corrected_rtest_pvalue = correct_pvalues(rtest_pvalue, correction_type)
    max_pvalue = [max(ctp, crp) for ctp, crp in zip(corrected_ttest_pvalue, corrected_rtest_pvalue)]

    # saving test result
    pk_file_name = os.path.split(pk_file_path)[1]
    test_result = \
        open(pk_file_name[:-4].replace('_MAvalues', '') + '_MAmotif_output.xls', 'w')
    header = _get_header(correction_type)
    test_result.write(header)

    idx_sorted = np.array(max_pvalue).argsort().tolist()
    for i in idx_sorted:
        line = lines[i]
        if ttest_pvalue[i] is not None:
            line += \
                str(t_stat[i]) + '\t' + str(ttest_pvalue[i]) + '\t' + \
                str(corrected_ttest_pvalue[i]) + '\t'
            line += \
                str(r_stat[i]) + '\t' + str(rtest_pvalue[i]) + '\t' + \
                str(corrected_rtest_pvalue[i]) + '\t'
            line += \
                str(max_pvalue[i]) + '\n'
        else:
            line += '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (None, None, None, None, None, None, None)
        test_result.write(line)


def correct_pvalues(pvalues, correction_type='benjamin'):
    sorted_indices = np.array(pvalues).argsort().tolist()
    corrected_pvalues = [None] * len(pvalues)
    if correction_type == 'benjamin':
        for rank, i in enumerate(sorted_indices):
            if pvalues[i] is None:
                continue
            corrected_pvalues[i] = len(pvalues) * pvalues[i] / (rank + 1)
            if corrected_pvalues[i] > 1:
                corrected_pvalues[i] = 1
    elif correction_type == 'bonferroni':
        for rank, i in enumerate(sorted_indices):
            if pvalues[i] is None:
                continue
            corrected_pvalues[i] = len(pvalues) * pvalues[i]
            if corrected_pvalues[i] > 1:
                corrected_pvalues[i] = 1
    else:
        corrected_pvalues = pvalues
    return corrected_pvalues
