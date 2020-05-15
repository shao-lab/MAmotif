.. _tutorial:

========
Tutorial
========

.. contents::
   :local:

Installation
============

Like many other Python packages and bioinformatics softwares, MAmotif can be
obtained easily from PyPI_ or Bioconda_.

Prerequisites
-------------

* Python >= 3.6
* MAnorm >= 1.3.0
* motifscan >= 1.2.0
* numpy >= 1.15
* scipy >= 1.0


Install with pip
----------------
The latest release of MAmotif is available at PyPI_, you can install via ``pip``::

    $ pip install mamotif

.. _PyPI: https://pypi.org/project/MAmotif/

Install with conda
------------------

You can also install MAmotif with conda_ through Bioconda_ channel::

   $ conda install -c bioconda mamotif

.. _conda: https://conda.io/docs/
.. _Bioconda: https://bioconda.github.io/


Usage of MAmotif
================

To check whether MAmotif is properly installed, you can inspect the version of
MAmotif by the ``-v/--version`` option::

  $ mamotif --version

Configuration
-------------

Before running MAmotif, you need to configure the genome and motif data files
for `MotifScan`:

Please refer to the QuickStart_ section of MotifScan for the details.

.. _QuickStart: https://motifscan.readthedocs.io/en/latest/quickstart.html

Run complete MAmotif workflow
-----------------------------

MAmotif provide a console script ``mamotif`` for running the program, the
``mamotif run`` sub-command is used to run complete MAmotif workflow
(MAnorm + MotifScan + Integration).

.. code-block:: shell

    $ mamotif run --p1 sampleA_peaks.bed --p2 sampleB_peaks.bed --r1 sampleA_reads.bed --r2 sampleB_reads.bed -g <genome>
    –m <motif_set> -o <output_dir>

.. tip::

    The ``run`` sub-command only provides basic MAnorm/MotifScan options.
    If you want to control other advanced options (MAnorm normalization
    options or MotifScan scanning options), please run them independently and
    call MAmotif integration module with the ``mamotif integrate`` sub-command.

Options
^^^^^^^

-h, --help           Show help message and exit.
--verbose            Enable verbose log output.
--p1, --peak1        **[Required]** Peak file of sample A.
--p2, --peak2        **[Required]** Peak file of sample B.
--pf, --peak-format  Format of the peak files. Default: bed
--r1, --read1        **[Required]** Read file of sample A.
--r2, --read2        **[Required]** Read file of sample B.
--rf, --read-format  Format of the read files. Default: bed
--n1, --name1        Name of sample A.
--n2, --name2        Name of sample B.
--s1, --shiftsize1   Single-end reads shiftsize of sample A. Default: 100
--s2, --shiftsize2   Single-end reads shiftsize of sample B. Default: 100
--pe, --paired-end   Paired-end mode.
-m                   **[Required]** Motif set to scan for.
-g                   **[Required]** Genome name.
-p                   P value cutoff for motif scores. Default: 1e-4
-t, --threads        Number of processes used to run in parallel.
--mode               Which sample to perform MAmotif on {both,A,B}. Default: both
--split              Split genomic regions into promoter/distal regions and
                     run separately.
--upstream           TSS upstream distance for promoters. Default: 4000
--downstream         TSS downstream distance for promoters. Default: 2000
--correction         Method for multiple testing correction {benjamin,bonferroni}.
                     Default: benjamin
-o, --output-dir     Directory to write output files.


Integrate MAnorm and MotifScan results
--------------------------------------

The ``mamotif integrate`` sub-command is used when users have already got the
MAnorm and MotifScan results, and only run the final integration procedure.

Suppose you have the MAnorm result (sample A vs sample B), and the MotifScan
results for both samples:

To find cell type-specific co-factors for sample A:

.. code-block:: shell

    $ mamotif integrate -i A_MAvalues.xls -m A_motifscan/motif_sites_numbers.xls -o <path>

Convert M=log2(A/B) to -M=log2(B/A) and find co-factors for sample B:

.. code-block:: shell

    $ mamotif integrate -i B_MAvalues.xls -m B_motifscan/motif_sites_numbers.xls -n -o <path>

Options
^^^^^^^

-h, --help        Show help message and exit.
--verbose         Enable verbose log output.
-i                MAnorm result for sample A or B (A/B_MAvalues.xls).
-m                MotifScan result for sample A or B (motif_sites_number.xls).
-n, --negative    Convert M=log2(A/B) to -M=log2(B/A). Required when finding
                  co-factors for sample B.
-g                Genome name. Required if `--split` is enabled.
--split           Split genomic regions into promoter/distal regions and  run separately.
--upstream        TSS upstream distance for promoters. Default: 4000
--downstream      TSS downstream distance for promoters. Default: 2000
--correction      Method for multiple testing correction {benjamin,bonferroni}.
                  Default: benjamin
-o, --output-dir  Directory to write output files.

MAmotif Output
==============

After finished running MAmotif, all output files will be written to the directory
you specified with "-o" argument.s

The MAmotif output table includes the following columns:

::

    1. Motif Name
    2. Target Number: Number of motif-present peaks
    3. Average of Target M values: Average M-value of motif-present peaks
    4. Std. of Target M values: M-value Std. of motif-present peaks
    5. Non-target Number: Number of motif-absent peaks
    6. Average of Non-target M-value: Average M-value of motif-absent peaks
    7. Std. of Non-target M-value: M-value Std. of motif-absent peaks
    8. T-test Statistics: T-Statistic for M-values of motif-present peaks against motif-absent peaks
    9. T-test P-value: Right-tailed P-value of T-test
    10. T-test P-value By Benjamin/Bonferrroni correction
    11. RanSum-test Statistic
    12. RankSum-test P-value
    13. RankSum-test P-value By Benjamin/Bonferroni correction
    14. Maximal P-value: Maximal corrected P-value of T-test and RankSum-test
