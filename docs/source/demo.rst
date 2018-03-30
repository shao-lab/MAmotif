.. _demo:

Example Usage
=============

Here we provide a step-by-step instruction on how to use MAmotif to find candidate cell-type specific regulators
associated with certain histone modifications.

We take the H3K4me3 analysis between adult and fetal ProES in MAmotif paper as an example:

1. Install MAmotif::

    $pip install mamotif
    $conda install -c bioconda mamotif

2. Download all data MAmotif needs::

    $mkdir MAmotif_demo
    $cd MAmotif_demo
    $wget ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM908nnn/GSM908038/suppl/GSM908038_H3K4me3-F_peaks.bed.gz
    $wget ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM908nnn/GSM908039/suppl/GSM908039_H3K4me3-A_peaks.bed.gz
    $wget ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM908nnn/GSM908038/suppl/GSM908038_H3K4me3-F.bed.gz
    $wget ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM908nnn/GSM908039/suppl/GSM908039_H3K4me3-A.bed.gz
    $gzip -d *gz

    Remove the header line and ribosomal reads (You do not need to do this for modern ChIP-seq mapping softwares)
    $sed -i '1d' GSM908038_H3K4me3-F.bed
    $sed -i '1d' GSM908039_H3K4me3-A.bed
    $sed -i '8986927,$d' GSM908039_H3K4me3-F.bed
    $sed -i '14916308,$d' GSM908039_H3K4me3-A.bed

    Substitute space into tab for bed files (You do not need to do this for your own bed files are tab-separated)
    $sed -i "s/ /\t/g" GSM908038_H3K4me3-F.bed
    $sed -i "s/ /\t/g" GSM908039_H3K4me3-A.bed


3. Build for genome sequences::

    $mkdir genome
    $cd genome
    $wget http://hgdownload.cse.ucsc.edu/goldenPath/hg18/bigZips/chromFa.zip
    $unzip chromFa.zip
    $cat *fa > hg18.fa
    $genomecompile -G hg18.fa -o hg18
    $cd ..

4. Build for motif PWM (Optional)

The motif matrix file which containing the motif score cutoff is already packaged under /data directory under MotifScan package.

If you want you compile for your custom motifs, please run the following commands::

    $mkdir motif
    $cd motif
    $wget http://jaspar2016.genereg.net/html/DOWNLOAD/JASPAR_CORE/pfm/nonredundant.tar.gz
    $tar -xzvf nonredundant.tar.gz
    $motifcompile -M nonredundant/pfm_vertebrates.txt -g ../genome/hg18 -o hg18_jaspar2016_nonredundant_vertebrates
    $cd ..

5. Run MAmotif::

   $mamotif --p1 GSM908039_H3K4me3-A_peaks.bed --p2 GSM908038_H3K4me3-F_peaks.bed --r1 GSM908039_H3K4me3-A.bed --r2 GSM908038_H3K4me3-F.bed -g genome/hg18 -m motif/hg18_jaspar2016_nonredundant_vertebrates_1e-4.txt -o AvsF_H3K4me3_MAmotif

6. Check the output of MAmotif