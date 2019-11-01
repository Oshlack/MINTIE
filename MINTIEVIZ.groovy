/*
Visualisation for the MINTIE pipeline

Author: Marek Cmero
*/

code_base = file(bpipe.Config.config.script).parentFile.absolutePath
load code_base + "/tools.groovy"
load code_base + "/references.groovy"

// initialise defaults if not provided
if(!binding.variables.containsKey("fastqCaseFormat")){
    fastqCaseFormat="cases/%_R*.fastq.gz"
}
if(!binding.variables.containsKey("fastqControlFormat")){
    fastqControlFormat="controls/%_R*.fastq.gz"
}
if(!binding.variables.containsKey("vcf") || !binding.variables.containsKey("tsv")){
    print "Please provide both the VCF and TSV file output from MINTIE's refine_contigs step."
    System.exit(0)
}

make_supertranscript_reference = {
    def sample_name = branch.name
    output.dir = sample_name
    produce(sample_name + "_supertranscript.fasta", sample_name + "_blocks_supertranscript.bed", sample_name + "_genes_supertranscript.bed"){
        exec """
        $python ${code_base}/collate/make_supertranscript.py $tsv $vcf \
            $tx_annotation $genome_fasta $output.dir $sample_name --log $output.dir/makest.log
        """, "make_supertranscript_reference"
    }
}

make_super_supertranscript = {
    def files = inputs.fasta.collect{ it as String }.collect{ it as File }
    def colpath = files.collect{ it.getName().split('_').first() }.join('_')
    def workingDir = System.getProperty("user.dir");
    input_beds = inputs.bed.collect{ it as String }
    block_beds = input_beds.findAll{ it.contains('_blocks_') }.join(' ')
    gene_beds = input_beds.findAll{ it.contains('_genes_') }.join(' ')
    colpath = workingDir + '/' + colpath + '_collated'
    output.dir = colpath
    produce("supersupertranscript.fasta", "supersupertranscript_blocks.bed", "supersupertranscript_genes.bed"){
        exec """
        cat $inputs.fasta | $python ${code_base}/util/remove_redundant_records.py - > $output.fasta ;
        cat $block_beds > $output2 ; cat $gene_beds > $output3
        """
    }
}

make_supertranscript_gmap_reference = {
    def collated_dir = new File(input.fasta).getParentFile().getName()
    output.dir = collated_dir + "/st_gmap_ref"
    produce("st_gmap_ref.chromosome"){
        exec """
        ${gmap}_build -s chrom -k 15 -d st_gmap_ref -D $collated_dir $input.fasta ;
        """
    }
}

align_contigs_to_supertranscript = {
    def colpath = new File(input.fasta).getParentFile().getName()
    def sample_name = branch.name.split('\\.').first() //remove branch dot suffix
    def index_dir = colpath
    output.dir = colpath + "/alignment"
    produce(sample_name + "_novel_contigs_st_aligned.sam"){
        exec """
        $gmap -D $index_dir -d st_gmap_ref -f samse -t $threads \
            -n 0 ${sample_name}/de_contigs.fasta > $output.sam ;
        """, "align_contigs_to_supertranscript"
    }
}

hisat_index = {
    def colpath = new File(input).getParentFile().getName()
    output.dir = colpath + "/genome"
    def idx_prefix = output.dir + "/hisat_index"
    produce("hisat_index.1.ht2") {
        exec """
        ${hisat}-build $input.fasta $idx_prefix
        """
    }
}

hisat_align = {
    def workingDir = System.getProperty("user.dir");
    def (rf1, rf2) = inputs.split().collect { workingDir+"/$it" }
    def sample_name = branch.name.split('\\.').first() //remove branch dot suffix
    def colpath = new File(input.fasta).getParentFile().getName()
    output.dir = colpath + "/alignment"
    if(type == "controls"){
        output.dir = colpath + "/alignment/controls/"
    }
    produce(sample_name + "_hisatAligned.sam"){
        exec """
        hisat_idx=$input.ht2; idx=\${hisat_idx%.1.ht2} ;
        $hisat --threads $threads -x \$idx -1 $rf1 -2 $rf2 > $output ;
        """, "hisat_align"
    }
}

sort_and_index_bam = {
    output.dir = new File(input.sam).getParentFile()
    transform('sam') to ('bam') {
        exec """
        $samtools sort -@ $threads -m ${sort_ram} $input.sam -o $output ;
        $samtools index $output
        """, "sort_and_index_bam"
    }
}

collate = {
    def colpath = new File(input.fasta).getParentFile().getName()
    def sample_name = branch.name.split('\\.').first() //remove branch dot suffix
    output.dir = colpath + "/results"
    produce(sample_name + '_results.tsv'){
        exec """
        $python ${code_base}/collate/collate.py \
            $sample_name \
            $sample_name/${sample_name}_blocks_supertranscript.bed \
            $colpath/alignment/${sample_name}_novel_contigs_st_aligned.bam \
            $colpath/alignment/${sample_name}_hisatAligned.bam \
            $input.tsv > $output
        """
    }
}

run { fastqCaseFormat * [ make_supertranscript_reference ] +
        make_super_supertranscript +
        make_supertranscript_gmap_reference +
        hisat_index +
            [ fastqCaseFormat * [ align_contigs_to_supertranscript + sort_and_index_bam,
                                  hisat_align + sort_and_index_bam],
              fastqControlFormat * [ hisat_align.using(type:"controls") + sort_and_index_bam ] ]
}
