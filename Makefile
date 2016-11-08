
all: parse_transcriptome_blat parse_genome_blat make_superReference

parse_transcriptome_blat: parse_transcriptome_blat.c
	g++ -O3 $^ -o $@

parse_genome_blat: parse_genome_blat.c
	g++ -O3 $^ -o $@

make_superReference: make_superReference.c
	g++ -O3 $^ -o $@

clean:
	-rm parse_transcriptome_blat parse_genome_blat make_superReference


