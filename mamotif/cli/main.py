"""
mamotif.cli.main
----------------

Main command line interface of MAmotif.
"""

import argparse
import os
from textwrap import dedent

from manorm.read import READ_FORMATS
from manorm.region import REGION_FORMATS

from mamotif import __version__
from mamotif.cli import intergrate, run
from mamotif.logging import setup_logger


def _existed_file(path):
    """Check whether a passed argument is an existed file."""
    if not os.path.isfile(path):
        raise argparse.ArgumentTypeError(f"file not found: {path}")
    return path


def _pos_int(value):
    """Check whether a passed argument is a positive integer."""
    try:
        value_int = int(value)
        if value_int <= 0:
            raise ValueError
    except (ValueError, TypeError):
        raise argparse.ArgumentTypeError(
            f"invalid positive int value: {value!r}")
    return value_int


def _add_verbose_argument(parser):
    parser.add_argument(
        "--verbose", dest="verbose", action="store_true", default=False,
        help="Enable verbose log messages.")
    return parser


def add_mamotif_arguments(parser):
    parser_integrate = parser.add_argument_group("Integration Options")
    parser_integrate.add_argument(
        "--split", dest="split", action="store_true", default=False,
        help="Split genomic regions into promoter/distal regions and run "
             "separately.")
    parser_integrate.add_argument(
        "--upstream", metavar="DISTANCE", dest="upstream",
        type=_pos_int, default=4000,
        help="TSS upstream distance for promoters. Default: 4000")
    parser_integrate.add_argument(
        "--downstream", metavar="DISTANCE", dest="downstream",
        type=_pos_int, default=2000,
        help="TSS downstream distance for promoters. Default: 2000")
    parser_integrate.add_argument(
        "--correction", dest="correction", choices=["benjamin", "bonferroni"],
        default="benjamin",
        help="Method for multiple testing correction. Default: benjamin")
    parser_output = parser.add_argument_group("Output Options")
    parser_output.add_argument(
        "-o", "--output-dir", metavar="DIR", dest="output_dir", required=True,
        help="Directory to write output files.")
    return parser


def configure_parser_run(subparsers):
    help_msg = "Run complete workflow (MAnorm + MotifScan + Integration)."
    desc_msg = help_msg + dedent("""

    Run the complete MAmotif workflow with basic MAnorm/MotifScan options.
    If you want to control other advanced options (MAnorm normalization 
    options or MotifScan scanning options), please run them independently and
    call MAmotif integration module with the `mamotif integrate` sub-command.
    """)
    epilog_msg = dedent("""
    Notes:
    ------
    Before running MAmotif, the MotifScan genome/motif data files should be 
    configured in advance. Please refer to the documentation for more 
    information.      
    """)
    parser = subparsers.add_parser(
        "run", description=desc_msg, help=help_msg, epilog=epilog_msg,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser_input = parser.add_argument_group("Input Options")
    parser_input.add_argument(
        "--p1", "--peak1", metavar="FILE", dest="peak_file1", required=True,
        type=_existed_file, help="Peak file of sample A.")
    parser_input.add_argument(
        "--p2", "--peak2", metavar="FILE", dest="peak_file2", required=True,
        type=_existed_file, help="Peak file of sample B.")
    parser_input.add_argument(
        "--pf", "--peak-format", metavar="FORMAT", dest="peak_format",
        choices=REGION_FORMATS, default="bed",
        help=f"Format of the peak files. Support {REGION_FORMATS}. "
             f"Default: bed")
    parser_input.add_argument(
        "--r1", "--read1", metavar="FILE", dest="read_file1", required=True,
        type=_existed_file, help="Read file of sample A.")
    parser_input.add_argument(
        "--r2", "--read2", metavar="FILE", dest="read_file2", required=True,
        type=_existed_file, help="Read file of sample B.")
    parser_input.add_argument(
        "--rf", "--read-format", metavar="FORMAT", dest="read_format",
        choices=READ_FORMATS, default="bed",
        help=f"Format of the read files. Support {READ_FORMATS}. Default: bed")
    parser_input.add_argument(
        "--n1", "--name1", metavar="NAME", dest="name1",
        help="Name of sample A. If not specified, the peak file name will be "
             "used.")
    parser_input.add_argument(
        "--n2", "--name2", metavar="NAME", dest="name2",
        help="Name of sample B. If not specified, the peak file name will be "
             "used.")

    parser_input.add_argument(
        "-m", "--motif", metavar="NAME", dest="motif", required=True,
        help="Motif set name to scan for.")
    parser_input.add_argument(
        "-g", "--genome", metavar="GENOME", dest="genome", required=True,
        help="Genome assembly name.")

    parser_reads = parser.add_argument_group("MAnorm Options")
    parser_reads.add_argument(
        "--s1", "--shiftsize1", metavar="N", dest="shift_size1",
        type=int, default=100,
        help="Single-end reads shift size for sample A. Reads are shifted by "
             "`N` bp towards 3' direction and the 5' end of each shifted read "
             "is used to represent the genomic locus of the DNA fragment. "
             "Set to 0.5 * fragment size of the ChIP-seq library. "
             "Default: 100")
    parser_reads.add_argument(
        "--s2", "--shiftsize2", metavar="N", dest="shift_size2",
        type=int, default=100,
        help="Single-end reads shift size for sample B. Default: 100")
    parser_reads.add_argument(
        "--pe", "--paired-end", dest="paired", action='store_true',
        default=False,
        help="Paired-end mode. The middle point of each read pair is used to "
             "represent the genomic locus of the DNA fragment. If specified, "
             "`--s1` and `--s2` will be ignored.")

    parser_motifscan = parser.add_argument_group("MotifScan Options")
    parser_motifscan.add_argument(
        "-p", dest="p_value", default="1e-4",
        choices=["1e-2", "1e-3", "1e-4", "1e-5", "1e-6"],
        help="P value cutoff for motif scores. Default: 1e-4")
    parser_motifscan.add_argument(
        "-t", "--threads", metavar="N", dest="n_threads", type=int, default=1,
        help="Number of processes used to run in parallel.")
    parser_mamotif = parser.add_argument_group("MAmotif Options")
    parser_mamotif.add_argument(
        "--mode", dest="mode", choices=['both', 'A', 'B'], default='both',
        help="Which sample to perform MAmotif on. Default: both")
    parser = add_mamotif_arguments(parser)
    parser = _add_verbose_argument(parser)
    parser.set_defaults(func=run.run)


def configure_parser_integrate(subparsers):
    help_msg = "Run the integration module with MAnorm and MotifScan results."
    desc_msg = help_msg + dedent("""

    This command is used when users have already got the MAnorm and MotifScan 
    results, and only run the final integration procedure.
    """)
    epilog_msg = dedent("""
    Examples:
    ---------
    Suppose you have the MAnorm result (sample A vs sample B), and the 
    MotifScan results for both samples:
    
    1) Find cell type-specific co-factors for sample A:
      
        mamotif integrate -i A_MAvalues.xls -m A_motifscan/motif_sites_numbers.xls -o <path>
        
    2) Convert M=log2(A/B) to -M=log2(B/A) and find co-factors for sample B:
      
        mamotif integrate -i B_MAvalues.xls -m B_motifscan/motif_sites_numbers.xls -n -o <path>
        
    """)
    parser = subparsers.add_parser(
        "integrate", description=desc_msg, help=help_msg, epilog=epilog_msg,
        formatter_class=argparse.RawDescriptionHelpFormatter)

    parser_input = parser.add_argument_group("Input Options")
    parser_input.add_argument(
        "-i", metavar="FILE", dest="f_manorm", required=True,
        help="MAnorm result for sample A or B (A/B_MAvalues.xls).")
    parser_input.add_argument(
        "-m", metavar="FILE", dest="f_motifscan", required=True,
        help="MotifScan result for sample A or B (motif_sites_number.xls).")
    parser_input.add_argument(
        "-n", "--negative", dest="negative", action="store_true",
        default=False,
        help="Convert M=log2(A/B) to -M=log2(B/A). Required when finding "
             "co-factors for sample B.")
    parser_input.add_argument(
        "-g", dest="genome", default=None,
        help="Genome name. Required if `--split` is enabled.")
    parser = add_mamotif_arguments(parser)
    parser = _add_verbose_argument(parser)
    parser.set_defaults(func=intergrate.run)


def configure_parser_main():
    """Configure the arguments parsers for MAmotif."""
    description = dedent("""
    MAmotif: An integrative toolkit for detecting cell type-specific regulators
    
    MAmotif is used to compare two ChIP-seq samples of the same protein from 
    different cell types (or conditions, e.g. wild-type vs mutant) and 
    identify transcriptional factors (TFs) associated with the cell type-biased 
    binding of this protein as its co-factors, by using TF binding information 
    obtained from motif analysis.
    
    Citation:
    Sun, H., Wang, J., Gong, Z. et al. Quantitative integration of epigenomic
    variation and transcription factor binding using MAmotif toolkit identifies
    an important role of IRF2 as transcription activator at gene promoters.
    Cell Discov 4, 38 (2018). https://doi.org/10.1038/s41421-018-0045-y
    """)

    epilog_msg = dedent("""
    Please run `mamotif COMMAND -h` to see the subcommand options.

    See also:
      Documentation: https://mamotif.readthedocs.io
      Source code: https://github.com/shao-lab/MAmotif
      Bug reports: https://github.com/shao-lab/MAmotif/issues
    """)

    parser = argparse.ArgumentParser(
        description=description, epilog=epilog_msg,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-v", "--version", action="version",
                        version=f"MAmotif {__version__}")

    subparsers = parser.add_subparsers(title="MAmotif Subcommands",
                                       metavar="command", dest="cmd")
    configure_parser_run(subparsers)
    configure_parser_integrate(subparsers)
    return parser


def main():
    parser = configure_parser_main()
    args = parser.parse_args()
    setup_logger(args.verbose)
    args.func(args)


if __name__ == '__main__':
    main()
