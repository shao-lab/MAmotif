"""MAmotif main script for running from the command line."""
import os
import sys
import argparse
import logging
from mamotif import __version__
from mamotif import workflow as mamotif_workflow


def argparser_config():
    """Configure the arguments parser.
    """
    description = """MAmotif -- An integrative toolkit for searching cell type-specific co-factors associated with differential binding."""
    parser = argparse.ArgumentParser(description=description, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--version', action='version', version="%(prog)s {}".format(__version__))
    group_input = parser.add_argument_group("Input File Arguments")
    group_input.add_argument('-p', dest='pk', help='MAnorm comparison peak file path')
    group_input.add_argument('-m', dest='motif', help='motifscan result file path')
    group_input.add_argument('-a', dest='refgene', default='', help='refgene file, default is none.')

    group_advanced = parser.add_argument_group("Advanced arguments")
    group_advanced.add_argument('-n', dest='negative', action='store_true', default=False,
                                help='Using negative test of this pk')
    group_advanced.add_argument('-c', dest='correction', default='benjamin',
                                help='correction type of pvalues, no correction or benjamin or bonferroni,'
                                     'default=benjamin')
    group_output = parser.add_argument_group("Output arguments")
    group_output.add_argument("-o", dest="output_dir", type=str, required=True,
                              help="Output directory.")

    return parser


def main():
    logging.basicConfig(level=logging.INFO, format="%(levelname)-8s @%(asctime)s: %(message)s", stream=sys.stderr,
                        datefmt="%m/%d/%Y %H:%M", filemode="w")
    parser = argparser_config()
    args = parser.parse_args()
    pk = args.pk
    refgene = os.path.abspath(args.refgene)
    motif = args.motif
    negative = args.negative
    correction = args.correction
    output_dir = args.output_dir
    mamotif_workflow.main(comparison_pk=pk, motifscan_result=motif, refgene_file=refgene, correction_type=correction,
                          neg=negative, output_dir=output_dir)


if __name__ == '__main__':
    main()
