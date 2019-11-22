/*
  __  __ ___ _   _ _____ ___ _____
 |  \/  |_ _| \ | |_   _|_ _| ____|
 | |\/| || ||  \| | | |  | ||  _|
 | |  | || || |\  | | |  | || |___
 |_|  |_|___|_| \_| |_| |___|_____|

Method for Inferring Novel Transcripts and Isoforms using Equivalences classes

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

if(!binding.variables.containsKey("assemblyFasta")){
    assemblyFasta=""
}

fastq_dedupe = {
    from("*.gz"){
        def sample_name = branch.name
        output.dir = sample_name
        produce(sample_name+'.1.fastq.gz',sample_name+'.2.fastq.gz'){
            exec """
            gunzip -c $input1.gz > $sample_name/temp1.fastq ;
            gunzip -c $input2.gz > $sample_name/temp2.fastq ;
            echo $sample_name/temp1.fastq > $sample_name/fastq.list ;
            echo $sample_name/temp2.fastq >> $sample_name/fastq.list ;
            echo "Reads before:" ; wc -l $sample_name/temp1.fastq ;
            $fastuniq -i $sample_name/fastq.list -o $output1.prefix -p $output2.prefix ;
            echo "Reads after:" ; wc -l $output1.prefix ;
            gzip $output1.prefix $output2.prefix ;
            rm $sample_name/fastq.list $sample_name/temp1.fastq $sample_name/temp2.fastq
            """, "fastq_dedupe"
        }
    }
}

trim = {
    output.dir = branch.name
    produce('trim1.fastq.gz', 'trim2.fastq.gz') {
        if (assemblyFasta != '') {
            // no need to trim if assembly provided
            exec """
            touch $output1 ; touch $output2
            """
        } else {
            exec """
            $trimmomatic PE -threads $threads -phred$scores $input1 $input2
                $output1.prefix /dev/null $output2.prefix /dev/null
                LEADING:$minQScore TRAILING:$minQScore MINLEN:$min_read_length ;
            gzip $output1.prefix $output2.prefix ;
            """, "trim"
        }
    }
}

assemble = {
    def sample_name = branch.name
    def Ks_for_soap = Ks.split(',').join(' ')
    output.dir = sample_name
    produce(sample_name+'_denovo_filt.fasta', sample_name+'.fasta'){
        if (assemblyFasta != '') {
            exec """
            ln -s $assemblyFasta $output1 ;
            cat $assemblyFasta $trans_fasta > $output2
            """
        } else if (assembler.toLowerCase() == 'trinity') {
            exec """
            $Trinity --seqType fq --max_memory ${assembly_mem}G --output $sample_name/trinity_assembly \
                --left $input1 --right $input2 --CPU $threads ;
            ln -s trinity_assembly/Trinity.fasta $output1 ;
            cat $output1 $trans_fasta > $output2
            """, "assemble"
        } else if (assembler.toLowerCase() == 'spades') {
            exec """
            $rnaspades -1 $input1 -2 $input2 -k $Ks -t $threads -m $assembly_mem -o $sample_name/SPAdes_assembly ;
            ln -s SPAdes_assembly/contigs.fasta $output1 ;
            cat $output1 $trans_fasta > $output2
            """, "assemble"
        } else {
            exec """
            if [ ! -d $output.dir/SOAPassembly ]; then mkdir $output.dir/SOAPassembly ; fi ;
            cd $sample_name/SOAPassembly ;

            echo \"max_rd_len=$max_read_length\" > config.config ;
            echo -e \"[LIB]\\nq1=../../$input1\\nq2=../../$input2\" >> config.config ;
            if [ -e SOAP.fasta ]; then rm SOAP.fasta ; fi ;
            for k in $Ks_for_soap ; do
                $soapdenovotrans pregraph -s config.config -o outputGraph_\$k -K \$k -p $threads ;
                $soapdenovotrans contig -g outputGraph_\$k ;
                cat outputGraph_\$k.contig | sed "s/^>/>k\${k}_/g" >> SOAP.fasta ;
            done ;
            cd ../../ ;
            $dedupe in=$sample_name/SOAPassembly/SOAP.fasta out=stdout.fa threads=$threads overwrite=true |
            $fasta_formatter |
            awk '!/^>/ { next } { getline seq } length(seq) > $max_read_length { print \$0 "\\n" seq }'
            > $output1 ;
            cat $output1 $trans_fasta > $output2 ;
            """, "assemble"
        }
    }
}

create_salmon_index = {
    def sample_name = branch.name
    def salmon_index = sample_name + "/all_fasta_index"
    output.dir = sample_name + "/all_fasta_index"
    produce('sa.bin', 'hash.bin', 'rsd.bin'){
        exec """
         $salmon index -t $input2 -i $salmon_index ;
         """
    }
}

run_salmon = {
    def workingDir = System.getProperty("user.dir");
    def (rf1, rf2) = inputs.split().collect { workingDir+"/$it" }
    def salmon_index="all_fasta_index"
    def base_outdir = "salmon_out"
    def controls_dir = fastqControlFormat.split("/")[-2]
    def sample_name = branch.name

    if(type == "controls"){
        sample_name = branch.parent.parent.name
        def control_name = branch.name

        output.dir = sample_name + "/" + controls_dir + "/" + control_name + "_salmon_out/aux_info"
        base_outdir = control_name + "_salmon_out"
        salmon_index = "../all_fasta_index"
    } else {
        output.dir = sample_name + "/salmon_out/aux_info"
    }

    produce("eq_classes.txt"){
        exec """
        cd $output.dir/../.. ;
        $salmon quant --dumpEq --seqBias -i $salmon_index -l A -r $rf1 $rf2 -p $threads -o $base_outdir
        """, "run_salmon"
    }
}

run_de = {
    def sample_name = branch.name
    output.dir = sample_name
    def test_flag = test_mode.toBoolean() ? "--test" : ""
    produce("eq_classes_de.txt"){
        exec """
        ${R}script $code_base/DE/compare_eq_classes.R $sample_name $input $trans_fasta $output --FDR=$fdr --minCPM=$min_cpm --minLogFC=$min_logfc $test_flag
        """, "run_de"
    }
}

create_ec_count_matrix = {
    def sample_name = branch.name
    def sample_names = inputs.split().collect { it.split('/')[-3].split('_salmon_out')[0] }
    sample_names.set(0, sample_name) // case sample, rest are controls
    sample_names = sample_names.join(',')
    output.dir = sample_name
    produce("ec_count_matrix.txt"){
        exec """
        $python $code_base/DE/create_ec_count_matrix.py $inputs $sample_names $output1 ;
        """, "create_ec_count_matrix"
    }
}

filter_on_significant_ecs = {
    def sample_name = branch.name
    output.dir = sample_name
    produce("de_contigs.fasta"){
        exec """
        $python $code_base/util/filter_fasta.py $input.fasta $input.txt --col_id contig > $output1 ;
        """
    }
}

align_contigs_against_genome = {
    def sample_name = branch.name
    output.dir = sample_name
    produce('aligned_contigs_against_genome.sam'){
        exec """
        $gmap -D $gmap_refdir -d $gmap_genome -f samse -t $threads -n 0 $input.fasta > $output
        """, "align_contigs_against_genome"
    }
}

annotate_contigs = {
    def sample_name = branch.name
    output.dir = sample_name
    produce("annotated_contigs.vcf", "annotated_contigs_info.tsv", "annotated_contigs.bam"){
        exec """
        $python ${code_base}/annotate/annotate_contigs.py \
            $sample_name $input.bam \
            $ann_info $tx_annotation \
            $output.bam $output.tsv \
            --minClip $min_clip \
            --minGap $min_gap \
            --minMatch $min_match \
            --log $output.dir/annotate.log > $output.vcf
        """
    }
}

refine_contigs = {
    def sample_name = branch.name
    output.dir = sample_name
    def motif_check = test_mode.toBoolean() ? "--skipMotifCheck" : ""
    produce("novel_contigs.vcf", "novel_contigs_info.tsv", "novel_contigs.bam", "novel_contigs.fasta"){
        exec """
        $python ${code_base}/annotate/refine_annotations.py \
            $input.tsv $input.vcf $input.bam $tx_annotation \
            $genome_fasta $output.prefix \
            --minClip $min_clip \
            --minGap $min_gap \
            $motif_check \
            --log $output.dir/refine.log > $output.vcf ;
        $samtools index $output.bam ;
        $python $code_base/util/filter_fasta.py $input.fasta $output.tsv --col_id contig_id > $output.fasta ;
        """
    }
}

salmon_quant = {
    def workingDir = System.getProperty("user.dir");
    def (rf1, rf2) = inputs.split().collect { workingDir+"/$it" }
    def sample_name = branch.name
    def salmon_index = sample_name + "/salmon_quant_index"
    output.dir = sample_name + "/salmon_quant_out"

    produce("quant.sf"){
        exec """
        cat $trans_fasta $input.fasta > $sample_name/salmon_quant_index.fasta ;
        $salmon index -t $sample_name/salmon_quant_index.fasta -i $salmon_index ;
        $salmon quant --seqBias --validateMappings -i $salmon_index -l A -r $rf1 $rf2 -p $threads -o $output.dir
        """, "salmon_quant"
    }
}

calculate_VAF = {
    output.dir = branch.name
    produce("vaf_estimates.txt"){
        exec """
        ${R}script ${code_base}/annotate/estimate_VAF.R $input.sf $input.tsv $tx2gene $output
        """
    }
}

post_process = {
    def sample_name = branch.name
    def var_filter = var_filter.split(',').join(' ')
    def gf_arg = gene_filter == '' ? '' : '--gene_filter ' + gene_filter
    def vf_arg = var_filter == '' ? '' : '--var_filter ' + var_filter
    output.dir = sample_name
    produce(sample_name + '_results.tsv'){
        exec """
        $python ${code_base}/annotate/post_process.py \
            $sample_name \
            $input.tsv \
            $sample_name/eq_classes_de.txt \
            $sample_name/vaf_estimates.txt \
            $gf_arg \
            $vf_arg \
            --log $output.dir/postprocess.log > $output
        """
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

run { fastqCaseFormat * [ fastq_dedupe +
                          trim +
                          assemble +
                          create_salmon_index +
                          [ fastqCaseFormat * [ run_salmon ],
                            fastqControlFormat * [ run_salmon.using(type:"controls") ] ] +
                          create_ec_count_matrix +
                          run_de +
                          filter_on_significant_ecs +
                          align_contigs_against_genome +
                          sort_and_index_bam +
                          annotate_contigs +
                          refine_contigs +
                          [ fastqCaseFormat * [ salmon_quant ] ] +
                          calculate_VAF +
                          post_process ]
}
