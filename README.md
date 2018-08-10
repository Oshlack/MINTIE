# CrypticVariantFinder #

A pipeline for detecting differentially expressed cryptic variants in RNA-seq data. CrypticVariantFinder can detect canonical and non-canonical fusions, as well as fusions with novel sequence between the boundary. The software is also able to detect insertions, deletions and small or large rearrangements in RNA. 

## How do I get set up? ##

Install the following dependencies:

* [trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)
* [BBmap](https://github.com/BioInfoTools/BBMap)
* [fastuniq](https://sourceforge.net/projects/fastuniq/)
* [fastx-toolkit](http://hannonlab.cshl.edu/fastx_toolkit/)
* [bedops](https://bedops.readthedocs.io/en/latest/)
* [bedtools](http://bedtools.readthedocs.io/en/latest/)
* [gmap](https://github.com/juliangehring/GMAP-GSNAP)
* [salmon](https://github.com/COMBINE-lab/salmon)
* [SOAPdenovo](http://soap.genomics.org.cn/soapdenovo.html)
* [samtools](http://samtools.sourceforge.net/)
* [star](https://github.com/alexdobin/STAR)
* [bpipe](https://github.com/ssadedin/bpipe) (>=0.9.9.5)

Then set all the paths correctly for `software' in crypticVariantFinder.groovy.

Now run:

```
git clone git@github.com:mcmero/CrypticVariant.git
cv_code=$PWD/CrypticVariant
```

That's it.

## How do I run it? ##

Make the following directories:

```
mkdir -p cases
mkdir -p controls
```

Now either copy or symlink your desired fastq files into the respective directories for your case and control samples. Cases are your tumour/cancer samples to look for crytic variants in, and the controls are the samples you are testing against. Ideally, controls will be benign tissue of the same type as the tumour primary -- however, as in blood-based cancers where normal of the same tissue type is difficult to acquire, remission samples, or tumours from other samples (ideally the same cancer type) can be used. More controls means more power, so aim for at least 5 but ideally more controls. 

Make sure your fastq files are in the following format:

    <sample>_<fastqinfo>_R1.fastq.gz

Ensure no underscores are present in your sample name, as the first part of the fastq filename will be assumed to be the sample name. 

Now make a shell script like the following:

```
#!/bin/sh

cases="cases/<sample>_<fastqinfo>_R1.fastq.gz cases/<sample>_<fastqinfo>_R2.fastq.gz"
controls="controls/<sample>_<fastqinfo>_R1.fastq.gz controls/<sample>_<fastqinfo>_R2.fastq.gz
          controls/<sample>_<fastqinfo>_R1.fastq.gz controls/<sample>_<fastqinfo>_R2.fastq.gz"

bpipe run -n 120 ${cv_code}/crypticVariantFinder.groovy $cases $controls
```

Cases are run on a 1 vs. all controls basis. Several cases can be specified and they will be run in parallel. Be careful when running >5 simultaneous cases as bpipe might start throwing errors about too many open files. 

If you are running the pipeline outside of the `${cv_code}` directory, you'll want to copy the bpipe.config file to the directory from which you are running. Make sure you update this with the correct modules and tweak any run parameters for your cluster scenario. Alternatively, you can run everything under the current environment (won't launch cluster jobs), by setting `executor=''`. 

Users might like to separate the pipeline into two parts: detection and visualisation. As the visualisation part is more time-consuming and aligns all the controls to the supertranscipt, you may want to analyse the results for significant contigs, then visualise only the samples you care about. Reducing the number of controls that you run the visualisation with is also a really good idea.

```
bpipe run -n 120 -u make_super_supertranscript /group/bioi1/marekc/20170918_cryptic_variant/CrypticVariant/crypticVariantFinder.groovy $cases $controls
```

Running the above in place of the run step will run everything for cases up to collating and visualising the results. Once you've done this, you can run without the -u make_super_supertranscript to finish up the pipeline.
