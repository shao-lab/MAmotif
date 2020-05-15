.. _demo:

Example Usage
=============

Here we provide a step-by-step instruction on how to use MAmotif to find candidate
cell-type specific regulators associated with certain histone modifications.

We take the H3K4me3 analysis between adult and fetal ProES in the MAmotif
paper as an example:

1. Install MAmotif::

    $pip install mamotif
    $conda install -c bioconda mamotif

2. Download example data files::

    $mkdir MAmotif_demo
    $cd MAmotif_demo
    $wget ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM908nnn/GSM908038/suppl/GSM908038_H3K4me3-F_peaks.bed.gz
    $wget ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM908nnn/GSM908039/suppl/GSM908039_H3K4me3-A_peaks.bed.gz
    $wget ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM908nnn/GSM908038/suppl/GSM908038_H3K4me3-F.bed.gz
    $wget ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM908nnn/GSM908039/suppl/GSM908039_H3K4me3-A.bed.gz
    $gzip -d *gz

    Remove the header line and ribosomal reads
    $sed -i '1d' GSM908038_H3K4me3-F.bed
    $sed -i '1d' GSM908039_H3K4me3-A.bed
    $sed -i '8986927,$d' GSM908038_H3K4me3-F.bed
    $sed -i '14916308,$d' GSM908039_H3K4me3-A.bed

    Substitute space into tab for bed files
    $sed -i "s/ /\t/g" GSM908038_H3K4me3-F.bed
    $sed -i "s/ /\t/g" GSM908039_H3K4me3-A.bed

.. note::

    The modification steps of data files above is specific to the example,
    since the format does not follow the standard. You do not have to do this
    for your own data.


3. Install genome `hg18` from UCSC database::

    $ motifscan genome --install -n hg18 -r hg18

4. Install motif PFMs set from JASPAR database::

    $ motifscan motif --install -n vertebrates -r vertebrates_non-redundant -g hg19

5. Run MAmotif::

   $mamotif run --p1 GSM908039_H3K4me3-A_peaks.bed --p2 GSM908038_H3K4me3-F_peaks.bed --r1 GSM908039_H3K4me3-A.bed --r2 GSM908038_H3K4me3-F.bed -g hg18 -m vertebrates -o AvsF_H3K4me3_MAmotif

6. Check the output of MAmotif
