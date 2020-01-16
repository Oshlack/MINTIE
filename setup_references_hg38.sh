#!/bin/bash

## This script will setup required hg38 references for running MINTIE

mkdir -p ref
cd ref

gmap_refdir=$PWD
commands="genome_fasta trans_fasta tx_annotation ann_info tx2gene gmap_refdir gmap_genome"

function genome_fasta_setup {
    file=hg38.fa
    wget http://ccb.jhu.edu/chess/data/hg38_p8.fa.gz
    gunzip hg38_p8.fa.gz
    samtools=`which samtools 2>/dev/null`
    if [ -z $samtools ] ; then
        samtools=$PWD/../tools/bin/samtools
    fi
    chr=$(seq 1 1 22)" X Y M"
    for i in $chr; do
        $samtools faidx hg38_p8.fa chr$i >> $file ;
    done
    if [ -s $file ]; then
        rm hg38_p8.fa*
        echo -e "$PWD/$file" > genome_fasta.success
    fi
}

function trans_fasta_setup {
    file=chess2.2.fa
    wget --no-check-certificate http://ccb.jhu.edu/software/stringtie/dl/gffread-0.11.6.Linux_x86_64.tar.gz
    tar -xvzf gffread-0.11.6.Linux_x86_64.tar.gz && rm gffread-0.11.6.Linux_x86_64.tar.gz
    gffread-0.11.6.Linux_x86_64/gffread chess2.2.gtf -g hg38.fa -w $file
    if [ -s $file ]; then
        echo -e "$PWD/$file" > trans_fasta.success
    fi
}

function tx_annotation_setup {
    file=chess2.2.gtf
    wget http://ccb.jhu.edu/chess/data/${file}.gz
    gunzip ${file}.gz
    grep -vE "^K|^chrUn|alt|random" $file > ${file}.tmp && mv ${file}.tmp $file

    if [ -s $file ]; then
        echo -e "$PWD/$file" > tx_annotation.success
    fi
}

function ann_info_setup {
    file=chess2.2.info
    python ../util/make_exon_reference.py chess2.2.gtf
    if [ -s $file ]; then
        echo -e "$PWD/$file" > ann_info.success
    fi
}

function tx2gene_setup {
    file=tx2gene.txt
    python ../util/make_tx2gene_lookup.py chess2.2.gtf > $file
    if [ -s $file ]; then
        echo -e "$PWD/$file" > tx2gene.success
    fi
}

function gmap_refdir_setup {
    mkdir -p $gmap_refdir && echo -e "$gmap_refdir" > gmap_refdir.success
}

function gmap_genome_setup {
    gmap_build=`which gmap_build 2>/dev/null`
    if [ -z $gmap_build ] ; then
        gmap_build=$PWD/../tools/bin/gmap_build
    fi
    $gmap_build -s chrom -k 15 -d gmap_genome -D $gmap_refdir hg38.fa
    if [ -d gmap_genome ]; then
        echo "gmap_genome" > gmap_genome.success
    fi
}

echo "// Path to references used by the MINTIE pipeline" > ../references.groovy
echo "gmap_refdir=\"$PWD/\"" >> ../references.groovy

for c in $commands ; do
    if [ ! -f ${c}.success ] ; then
        echo "$c not found, setting this up..."
        ${c}_setup
    fi
    c_path=`cat ${c}.success`
    echo "$c=\"$c_path\"" >> ../references.groovy
done

#loop through commands to check they are all installed
echo "Checking that all required references were setup:"
Final_message="All commands installed successfully!"
for c in $commands ; do
    if [ ! -f ${c}.success ] ; then
        echo -n "WARNING: $c could not be found!!!! "
        echo "You will need to setup $c manually, then add its path to references.groovy"
    Final_message="WARNING: One or more command did not install successfully. See warning messages above. \
                       You will need to correct this before running MINTIE."
    else
        echo "$c looks like it has been setup"
    fi
done
echo "**********************************************************"
echo $Final_message
