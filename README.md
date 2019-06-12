# MINTIE #

A pipeline for detecting differentially expressed cryptic variants in RNA-seq data. MINTIE can detect canonical and non-canonical fusions, as well as fusions with novel sequence between the boundary. The software is also able to detect insertions, deletions and small or large rearrangements in RNA. 

## How do I get set up? ##

```
git clone git@github.com:mcmero/MINTIE.git
cd MINTIE
chmod u+x install_linux64.sh
./install_linux64.sh
```

This will install the dependencies required to run the pipeline. 

Make sure you have [python 3+](https://www.python.org/downloads/) installed, then run the following:

```
pip install -r requirements.txt
```

Make sure you also have [R](https://www.r-project.org/) v3.2+ installed with the following dependencies:

```
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("edgeR")
BiocManager::install("GenomicRanges")

install.packages("data.table")
install.packages("reshape2")
install.packages("dplyr")
install.packages("stringr")
install.packages("seqinr")
```

### References ###

You can quickly set up the references required if you're running hg38:

```
mkdir -p ref && cd ref
wget ftp://ftp.ensembl.org/pub/release-84/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz
wget ftp://ftp.ensembl.org/pub/release-84/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz
wget http://ccb.jhu.edu/chess/data/chess2.2.gtf.gz

for i in `ls *gz`; do 
	gunzip $i
done
```

MINTIE needs to set up a specific exon reference, so go ahead and run this script in the MINTIE directory:

```
python util/make_exon_reference.py ref/chess2.2.gtf
```

The last step is to create a GMAP reference:

```
tools/bin/gmap_build -s chrom -k 15 -d hg38 -D ref/gmap_genome ref/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa
```

Finally, check that all references under your `tools.groovy` and `references.groovy` files are correct.

## How do I run it? ##

Make the following directories:

```
mkdir -p cases
mkdir -p controls
```

Now either copy or symlink your desired fastq files into the respective directories for your case and control samples. Cases are your tumour/cancer samples to look for crytic variants in, and the controls are the samples you are testing against. Ideally, controls will be benign tissue of the same type as the tumour primary -- however, as in blood-based cancers where normal of the same tissue type is difficult to acquire, remission samples, or tumours from other samples (ideally the same cancer type) can be used. More controls means more power, so aim for at least 5 but ideally more controls. 

Make sure your fastq files match the pattern found in your `params.txt` file: 

```
-p fastqCaseFormat=cases/%_R*.fastq.gz
-p fastqControlFormat=controls/%_R*.fastq.gz
```

The above is the default and indicates the directories of cases and controls, and that the fastq files are in `<sample>_R*.fastq.gz` format, where the `*` refers to the paired read number (i.e. R1 and R2). 

MINTIE uses [bpipe](https://github.com/ssadedin/bpipe). Please see documentation for bpipe under [http://docs.bpipe.org/](http://docs.bpipe.org/) for more information. 

You should double-check the parameters in the `params.txt` file, and adjust as needed. You can also copy this file and create a custom parameters file (make sure to use `bpipe run @<your param file>` in this case). If you are running MINTIE on a cluster, also check the `bpipe.config` file and check that the executor, queue, processor and memory settings are set according to your preferences.

Now make a shell script like the following:

```
#!/bin/sh

cases=`ls cases/*fastq.gz`
controls=`ls controls/*fastq.gz`

bpipe run @params.txt MINTIE.groovy $cases $controls
```

Cases are run on a 1 vs. all controls basis. Several cases can be specified and they will be run in parallel. Be careful when running >5 simultaneous cases as bpipe might start throwing errors about too many open files. 


Users might like to separate the pipeline into two parts: detection and visualisation. As the visualisation part is more time-consuming and aligns all the controls to the supertranscipt, you may want to analyse the results for significant contigs, then visualise only the samples you care about. Reducing the number of controls that you run the visualisation with is also a really good idea.

```
bpipe run @params.txt -u make_super_supertranscript MINTIE.groovy $cases $controls
```

Running the above in place of the run step will run everything for cases up to collating and visualising the results. Once you've done this, you can run without the -u make_super_supertranscript to finish up the pipeline.
