#!/bin/bash

## This script will setup required hg38 references for running MINTIE

mkdir -p ref
cd ref

gmap_refdir=$PWD
commands="genome_fasta trans_fasta tx_annotation ann_info tx2gene gmap_refdir gmap_genome"

function genome_fasta_setup {
    file=Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa
    wget ftp://ftp.ensembl.org/pub/release-84/fasta/homo_sapiens/dna/${file}.gz
    gunzip ${file}.gz
    if [ -f $file ]; then
        echo -e "$PWD/$file" > genome_fasta.success
    fi
}

function trans_fasta_setup {
    file=Homo_sapiens.GRCh38.cdna.all.fa
    wget ftp://ftp.ensembl.org/pub/release-84/fasta/homo_sapiens/cdna/${file}.gz
    gunzip ${file}.gz
    if [ -f $file ]; then
        echo -e "$PWD/$file" > trans_fasta.success
    fi
}

function tx_annotation_setup {
    file=chess2.2.gtf
    wget http://ccb.jhu.edu/chess/data/${file}.gz
    gunzip ${file}.gz
    if [ -f $file ]; then
        echo -e "$PWD/$file" > tx_annotation.success
    fi
}

function ann_info_setup {
    file=chess2.2.info
    python ../util/make_exon_reference.py chess2.2.gtf
    if [ -f $file ]; then
        echo -e "$PWD/$file" > ann_info.success
    fi
}

function tx2gene_setup {
    file=tx2gene.txt
    python ../util/make_tx2gene_lookup.py chess2.2.gtf > $file
    if [ -f $file ]; then
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
    $gmap_build -s chrom -k 15 -d gmap_genome -D $gmap_refdir Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa
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
