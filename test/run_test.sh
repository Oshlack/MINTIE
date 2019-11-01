#!/bin/sh

# modify reference to control fasta
grep -v allvars-control references.groovy > backup_references.groovy
echo "trans_fasta=\"$PWD/test/data/allvars-control.fasta\"" > tmp.txt
cat backup_references.groovy tmp.txt > references.groovy ; rm tmp.txt

cases=`ls test/data/cases/*gz`
controls=`ls test/data/controls/*gz`
tools/bin/bpipe @test/test_params.txt MINTIE.groovy $cases $controls
tools/bin/bpipe @test/test_params.txt -p tsv="allvars-case/novel_contigs_info.tsv" -p vcf="allvars-case/novel_contigs.vcf" MINTIEVIZ.groovy $cases $controls

# restore original references
mv backup_references.groovy references.groovy
