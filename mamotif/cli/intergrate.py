"""
mamotif.cli.integrate
---------------------

Integrate MAnorm and MotifScan results to run MAmoitf.
"""

from manorm.logging import setup_logger as setup_manorm_logger
from motifscan.logging import setup_logger as setup_motifscan_logger

from mamotif.integration import run_integration


def run(args):
    setup_manorm_logger(args.verbose)
    setup_motifscan_logger(args.verbose)
    run_integration(
        f_manorm=args.f_manorm, f_motifscan=args.f_motifscan,
        negative=args.negative, genome=args.genome, split=args.split,
        upstream=args.upstream, downstream=args.downstream,
        correction=args.correction, output_dir=args.output_dir)
