MAmotif
=======

|travis-ci| |Documentation Status| |pypi| |license|

.. |travis-ci| image:: https://travis-ci.org/shao-lab/MAmotif.svg?branch=master
   :target: https://travis-ci.org/shao-lab/MAmotif
.. |Documentation Status| image:: https://readthedocs.org/projects/mamotif/badge/?version=latest
   :target: http://mamotif.readthedocs.io/en/latest/?badge=latest
.. |pypi| image:: https://img.shields.io/pypi/v/mamotif.svg
   :target: https://pypi.python.org/pypi/mamotif
.. |license| image:: https://img.shields.io/pypi/l/MAmotif.svg
   :target: https://github.com/shao-lab/MAmotif/blob/master/LICENSE

Introduction
------------

**MAmotif** is used to compare two ChIP-seq samples of the same protein from different cell types or conditions
(e.g. Mutant vs Wild-type) and **identify transcriptional factors (TFs) associated with the cell-type biased binding**
of this protein as its **co-factors**, by using TF binding information obtained from motif analysis
(or from other ChIP-seq data).

MAmotif automatically combines **MAnorm** model to perform quantitative comparison on given ChIP-seq samples together
with Motif-Scan toolkit to scan ChIP-seq peaks for **TF binding motifs**, and uses a systematic integrative analysis to
search for TFs whose binding sites are significantly associated with the cell-type biased peaks between two ChIP-seq samples.

When applying to ChIP-seq data of histone marks of regulatory elements (such as H3K4me3 for active promoters and
H3K9/27ac for active promoter/enhancers), or DNase/ATAC-seq data, MAmotif can be used to detect **cell-type specific regulators**.

Workflow
--------

.. image:: https://github.com/shao-lab/MAmotif/blob/master/docs/source/image/MAmotif_workflow.png

Documentation
-------------

To see the full documentation of MAmoitf, please refer to: http://bioinfo.sibs.ac.cn/shaolab/mamotif/index.php

Installation
------------

The latest release of MAmotif is available at `PyPI <https://pypi.python.org/pypi/mamotif>`__:

::

    $ pip install mamotif

Or you can install MAmotif via conda:

**WIP!**

::

    $ conda install -c bioconda mamotif

MAmotif uses `setuptools <https://setuptools.readthedocs.io/en/latest/>`__ for installation from source code.
The source code of MAmotif is hosted on GitHub: https://github.com/shao-lab/MAmotif

You can clone the repo and execute the following command under source directory:

::

    $ python setup.py install

Galaxy Installation
-------------------

**WIP!**


Usage
-----

You need to build some prerequisites before running MAmotif:

Build genomes
^^^^^^^^^^^^^

Preprocess sequences and genome-wide nucleotide frequency for the corresponding genome assembly.

::

    $ genomecompile [-h] [-v] -G hg19.fa -o hg19_genome

**Note:** You only need to run this command once for each genome

Build motifs (Optional)
^^^^^^^^^^^^^^^^^^^^^^^

**Note:** MAmotif provides some preprocessed motif PWM files under **data/motif** of the MotifScan package.

Build motif PWM/motif-score cutoff for custom motifs that are not included in our pre-complied motif collection:

::

    $ motifcompile [-h] [-v] –M motif_pwm_demo.txt –g hg19_genome -o hg19_motif

run MAmotif
^^^^^^^^^^^

::

    $ mamotif --p1 sample1_peaks.bed --p2 sample2_peaks.bed --r1 sample1_reads.bed --r2 sample2_reads.bed -g hg19_genome
    –m hg19_motif_p1e-4.txt -o sample1_vs_sample2

**Note:** Using -h/--help for the details of all arguments.

Output of MAmotif
-----------------

After finished running MAmotif, all output files will be written to the directory you specified with "-o" argument.

Main output
^^^^^^^^^^^

::

    1.Motif Name
    2.Target Number: Number of motif-present peaks
    3.Average of Target M-value: Average M-value of motif-present peaks
    4.Deviation of Target M-value: M-value Std of motif-present peaks
    5.Non-target Number: Number of motif-absent peaks
    6.Average of Non-target M-value: Average M-value of motif-absent peaks
    7.Deviation of Non-target M-value: M-value Std of motif-absent peaks
    8.T-test Statistics: T-Statistics for M-values of motif-present peaks against motif-absent peaks
    9.T-test P-value: Right-tailed P-value of T-test
    10.T-test P-value By Benjamin correction
    11.RanSum-test Statistics
    12.RankSum-test P-value
    13.RankSum-test P-value By Benjamin correction
    14.Maximal P-value: Maximal corrected P-value of T-test and RankSum-test

MAnorm output
^^^^^^^^^^^^^

MAmotif will invoke MAnorm and output the normalization results and MA-plot for samples under comparison.


Motif output
^^^^^^^^^^^^

MAmotif will also output tables to summarize the enrichment of motifs and the motif target number and motif-score
of each peak region.

If you specified "-s" with MAmotif, it will also output the genome coordinates of every motif target site.


License
-------

`BSD 3-Clause
License <https://github.com/shao-lab/MAmotif/blob/master/LICENSE>`__


