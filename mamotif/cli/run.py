"""
mamotif.cli.run
---------------

Run complete MAmotif workflow (MAnorm + MotifScan + Integration).
"""

import logging
import os

import manorm.cli as cli_manorm
import motifscan.cli.main as cli_motifscan

from mamotif.integration import run_integration

logger = logging.getLogger(__name__)


def run_manorm_from_mamotif(args):
    args.name1 = args.name1 or os.path.splitext(
        os.path.basename(args.peak_file1))[0]
    args.name2 = args.name2 or os.path.splitext(
        os.path.basename(args.peak_file2))[0]
    manorm_dir = os.path.abspath(
        os.path.join(args.output_dir,
                     f'{args.name1}_vs_{args.name2}_manorm_output'))
    args_manorm = [
        "--p1", args.peak_file1, "--p2", args.peak_file2,
        "--pf", args.peak_format, "--r1", args.read_file1,
        "--r2", args.read_file2, "--rf", args.read_format,
        "--n1", args.name1, "--n2", args.name2,
        "--s1", args.shift_size1, "--s2", args.shift_size2,
        "--wa", "-o", manorm_dir]
    if args.paired:
        args_manorm.append("--pe")

    parser_manorm = cli_manorm.configure_parser()
    args_manorm = parser_manorm.parse_args(map(str, args_manorm))
    args_manorm = cli_manorm.preprocess_args(args_manorm)
    cli_manorm.run(args_manorm)
    return manorm_dir


def run_motifscan_from_mamotif(args, f_manorm):
    sample_name = os.path.basename(f_manorm).replace('_MAvalues.xls', '')
    motifscan_dir = os.path.abspath(
        os.path.join(args.output_dir, f'{sample_name}_motifscan_output'))
    args_motifscan = [
        "scan", "-i", f_manorm, "-f", "manorm", "--motif", args.motif,
        "--genome", args.genome, "-p", args.p_value, "-t", args.n_threads,
        "--no-enrich", "-o", motifscan_dir]

    parser_motifscan = cli_motifscan.configure_parser_main()
    args_motifscan = parser_motifscan.parse_args(map(str, args_motifscan))
    args_motifscan.func(args_motifscan)
    return os.path.join(motifscan_dir, "motif_sites_number.xls")


def run(args):
    cli_manorm.setup_logger(args.verbose)
    cli_motifscan.setup_logger(args.verbose)
    manorm_dir = run_manorm_from_mamotif(args)
    if args.mode in ['both', 'A']:
        logger.info("\nScanning motifs for sample A")
        f_manorm = os.path.join(manorm_dir, f"{args.name1}_MAvalues.xls")
        f_motifscan = run_motifscan_from_mamotif(args, f_manorm)
        logger.info("Running MAmotif for sample A")
        run_integration(f_manorm=f_manorm, f_motifscan=f_motifscan,
                        negative=False, genome=args.genome, split=args.split,
                        upstream=args.upstream, downstream=args.downstream,
                        correction=args.correction, output_dir=args.output_dir)
    if args.mode in ['both', 'B']:
        logger.info("\nScanning motifs for sample B")
        f_manorm = os.path.join(manorm_dir, f"{args.name2}_MAvalues.xls")
        f_motifscan = run_motifscan_from_mamotif(args, f_manorm)
        logger.info("Running MAmotif for sample B")
        run_integration(f_manorm=f_manorm, f_motifscan=f_motifscan,
                        negative=True, genome=args.genome, split=args.split,
                        upstream=args.upstream, downstream=args.downstream,
                        correction=args.correction, output_dir=args.output_dir)
