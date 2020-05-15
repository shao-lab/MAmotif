MAmotif
=======

|Documentation Status| |pypi| |license|

.. |Documentation Status| image:: https://readthedocs.org/projects/mamotif/badge/?version=latest
   :target: https://mamotif.readthedocs.io/en/latest/?badge=latest
.. |pypi| image:: https://img.shields.io/pypi/v/mamotif.svg
   :target: https://pypi.org/project/MAmotif/
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

Citation
--------

`Sun H, Wang J, Gong Z, Yao J, Wang Y, Xu J, Yuan GC, Zhang Y, Shao Z. Quantitative integration of epigenomic variation
and transcription factor binding using MAmotif toolkit identifies an important role of IRF2 as transcription activator
at gene promoters. Cell discovery. 2018 Jul 10;4(1):38. <https://www.nature.com/articles/s41421-018-0045-y>`__

Installation
------------

The latest release of MAmotif is available at `PyPI <https://pypi.org/project/MAmotif/>`__:

.. code-block:: shell

    $ pip install mamotif

Or you can install MAmotif via conda:

.. code-block:: shell

    $ conda install -c bioconda mamotif


Documentation
-------------

To see the full documentation of MAmotif, please refer to: http://mamotif.readthedocs.io/en/latest/

License
-------

`BSD 3-Clause
License <https://github.com/shao-lab/MAmotif/blob/master/LICENSE>`__


