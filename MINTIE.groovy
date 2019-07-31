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
            """
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
            """
        } else if (assembler.toLowerCase() == 'spades') {
            exec """
            $rnaspades -1 $input1 -2 $input2 -k $Ks -t $threads -m $assembly_mem -o $sample_name/SPAdes_assembly ;
            ln -s SPAdes_assembly/contigs.fasta $output1 ;
            cat $output1 $trans_fasta > $output2
            """
        } else {
            exec """
            if [ ! -d $output.dir/SOAPassembly ]; then mkdir $output.dir/SOAPassembly ; fi ;
            cd $sample_name/SOAPassembly ;

            echo \"max_rd_len=$max_read_length\" > config.config ;
            echo -e \"[LIB]\\nq1=../../$input1\\nq2=../../$input1\" >> config.config ;
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
            ""","assemble"
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
    def sample_name = branch.parent.parent.name.split('_').first()

    if(type == "controls"){
        sample_name = branch.parent.parent.name.split('_').first()
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
    test_flag = test_mode.toBoolean() ? "--test" : ""
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
    motif_check = test_mode.toBoolean() ? "--skipMotifCheck" : ""
    produce("novel_contigs.vcf", "novel_contigs_info.tsv", "novel_contigs.bam"){
        exec """
        $python ${code_base}/annotate/refine_annotations.py \
            $input.tsv $input.vcf $input.bam $tx_annotation \
            $genome_fasta $output.tsv $output.bam \
            --minClip $min_clip \
            --minGap $min_gap \
            $motif_check \
            --log $output.dir/refine.log > $output.vcf ;
        $samtools index $output.bam
        """
    }
}

create_supertranscript_reference = {
    def sample_name = branch.name
    output.dir = sample_name
    produce(sample_name + "_supertranscript.fasta"){
        exec """
        $python ${code_base}/annotate/make_supertranscript.py $input.tsv $input.vcf \
            $tx_annotation $genome_fasta $output.dir $sample_name --log $output.dir/makest.log
        """, "create_supertranscript_reference"
    }
}

make_super_supertranscript = {
    def workingDir = System.getProperty("user.dir");
    colpath = inputs.collect { it.split('/').last().split('_').first() }.join('_')
    colpath = workingDir + '/' + colpath + '_collated'
    output.dir = colpath
    produce('supersupertranscript.fasta'){
        exec """
        cat $inputs.fasta | $python ${code_base}/util/remove_redundant_records.py - > $output ;
        """
    }
}

make_supertranscript_gmap_reference = {
    output.dir = colpath
    produce("st_gmap_ref"){
        exec """
        ${gmap}_build -s chrom -k 15 -d st_gmap_ref -D $output.dir $input.fasta ;
        """
    }
}

align_contigs_to_supertranscript = {
    output.dir = colpath+"/alignment"
    index_dir = colpath
    def sample_name = branch.name.split('\\.').first() //remove branch dot suffix
    produce(sample_name+"_novel_contigs_st_aligned.sam"){
        exec """
        $gmap -D $index_dir -d st_gmap_ref -f samse -t $threads \
            -n 0 ${sample_name}/de_contigs.fasta > $output.sam ;
        """, "align_contigs_to_supertranscript"
    }
}

hisat_index = {
    output.dir=colpath+"/genome"
    def idx_prefix = output.dir+"/hisat_index"
    produce("hisat_index.1.ht2") {
        exec """
        ${hisat}-build $input.fasta $idx_prefix
        """
    }
}

hisat_align = {
    output.dir = colpath+"/alignment"
    def workingDir = System.getProperty("user.dir");
    def (rf1, rf2) = inputs.split().collect { workingDir+"/$it" }
    def sample_name = branch.name.split('\\.').first() //remove branch dot suffix
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
        $samtools sort -@ $threads -m ${sort_ram}G $input.sam -o $output ;
        $samtools index $output
        """, "sort_and_index_bam"
    }
}

post_process = {
    output.dir = colpath + '/results'
    def sample_name = branch.name.split('\\.').first() //remove branch dot suffix
    def var_filter = var_filter.split(',').join(' ')
    gf_arg = gene_filter == '' ? '' : '--gene_filter ' + gene_filter
    produce(sample_name + '_results.tsv'){
        exec """
        $python ${code_base}/annotate/post_process.py $sample_name \
            $sample_name/novel_contigs_info.tsv \
            $sample_name/eq_classes_de.txt \
            $sample_name/${sample_name}_blocks_supertranscript.bed \
            $colpath/alignment/${sample_name}_novel_contigs_st_aligned.bam \
            $colpath/alignment/${sample_name}_hisatAligned.bam \
            $gf_arg \
            --var_filter $var_filter > $output
        """
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
                          create_supertranscript_reference ] +
        make_super_supertranscript +
        make_supertranscript_gmap_reference +
        hisat_index +
            [ fastqCaseFormat * [ align_contigs_to_supertranscript + sort_and_index_bam,
                                  hisat_align + sort_and_index_bam],
              fastqControlFormat * [ hisat_align.using(type:"controls") + sort_and_index_bam ] ] +
              fastqCaseFormat * [ post_process ]
}
