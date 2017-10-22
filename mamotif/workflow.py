import os
import logging
from mamotif.lib.MAmotifIO import motif_classified_pks_ttest
from mamotif.lib.MAmotifIO import classify_MAnorm_pks_by_promoter


def main(comparison_pk, motifscan_result, refgene_file='', correction_type='benjamin', neg=False, output_dir=None):
    comparison_pk = os.path.abspath(comparison_pk)
    motifscan_result = os.path.abspath(motifscan_result)
    if not os.path.isdir(output_dir):
        os.mkdir(output_dir)
    os.chdir(output_dir)
    # classify pk into promoter_pk and distal-pk
    if refgene_file == '':
        logging.info("Performing MAmotif...")
        motif_classified_pks_ttest(
            comparison_pk, motifscan_result, correction_type=correction_type, negative=neg)
    else:
        logging.info("Classifying peaks into Promoter/Distal peaks...")
        promoter_pk, distal_pk = classify_MAnorm_pks_by_promoter(comparison_pk, refgene_file)
        # motifs classified peaks test for pk, promoter_pk and distal_pk
        logging.info("Performing MAmotif...")
        motif_classified_pks_ttest(
            comparison_pk, motifscan_result, correction_type=correction_type, negative=neg)
        logging.info("Performing MAmotif on Promoter peaks...")
        motif_classified_pks_ttest(
            promoter_pk, motifscan_result, correction_type=correction_type, negative=neg)
        logging.info("Performing MAmotif on Distal peaks...")
        motif_classified_pks_ttest(
            distal_pk, motifscan_result, correction_type=correction_type, negative=neg)
