#!/bin/sh

cases=`ls test/data/cases/*gz`
controls=`ls test/data/controls/*gz`

if [ -s allvars-case/novel_contigs_info.tsv ] && [ -s allvars-case/novel_contigs.vcf ]; then
    tools/bin/bpipe @test/test_params.txt -p tsv="allvars-case/novel_contigs_info.tsv" -p vcf="allvars-case/novel_contigs.vcf" MINTIEVIZ.groovy $cases $controls
else
    echo "ERROR: Could not run test. Input files are missing. Check that MINTIE.groovy ran successfully."
fi
