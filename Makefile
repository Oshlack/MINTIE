
all: parse_transcriptome_blat parse_genome_blat make_superReference psl2gtf psl2sjdbFileChrStartEnd

parse_transcriptome_blat: parse_transcriptome_blat.c
	g++ -O3 $^ -o $@

parse_genome_blat: parse_genome_blat.c
	g++ -O3 $^ -o $@

make_superReference: make_superReference.c
	g++ -O3 $^ -o $@

psl2gtf: psl2gtf.c
	g++ -O3 $^ -o $@

psl2sjdbFileChrStartEnd: psl2sjdbFileChrStartEnd.c
	g++ -O3 $^ -o $@

clean:
	-rm parse_transcriptome_blat parse_genome_blat make_superReference


