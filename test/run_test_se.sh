#!/bin/sh

# modify reference to control fasta
grep -v allvars-control references.groovy > backup_references.groovy
echo "trans_fasta=\"$PWD/test/data/allvars-control.fasta\"" > tmp.txt
cat backup_references.groovy tmp.txt > references.groovy ; rm tmp.txt

cases=`ls test/data/cases/*gz | grep R1`
controls=`ls test/data/controls/*gz | grep R1`
bpipe @test/test_params.txt MINTIE_SE.groovy $cases $controls

# restore original references
mv backup_references.groovy references.groovy
