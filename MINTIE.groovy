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
if(!binding.variables.containsKey("run_de_step")){
    run_de_step="true"
}
if(!binding.variables.containsKey("splice_motif_mismatch")){
    splice_motif_mismatch=0
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
    def Ks_for_soap = Ks.toString().contains(',') ? Ks.split(',').join(' ') : Ks
    output.dir = sample_name
    produce(sample_name + '_denovo_filt.fasta'){
        if (assemblyFasta != '') {
            exec """
            ln -s $assemblyFasta $output ;
            """
        } else if (assembler.toLowerCase() == 'trinity') {
            exec """
            $Trinity --seqType fq --max_memory ${assembly_mem}G --output $sample_name/trinity_assembly \
                --left $input1 --right $input2 --CPU $threads ;
            ln -s trinity_assembly/Trinity.fasta $output ;
            """, "assemble"
        } else if (assembler.toLowerCase() == 'spades') {
            exec """
            $rnaspades -1 $input1 -2 $input2 -k $Ks -t $threads -m $assembly_mem -o $sample_name/SPAdes_assembly ;
            ln -s SPAdes_assembly/contigs.fasta $output ;
            """, "assemble"
        } else {
            exec """
            rlens=`zcat $input1 $input2 \
                       | awk -v mrl=$min_read_length 'BEGIN {minlen = mrl; maxlen = 0} {
                            if (NR % 4 == 2) {
                                rlen = length(\$1) ;
                                if (rlen > maxlen) {maxlen = rlen}
                                if (rlen < minlen) {minlen = rlen}
                            }} END {print minlen" "maxlen}'` ;
            min_rlen=\${rlens% *} ;
            max_rlen=\${rlens#* } ;

            if [ ! -d $output.dir/SOAPassembly ]; then
                mkdir $output.dir/SOAPassembly ;
            fi ;
            cd $output.dir/SOAPassembly ;

            echo \"max_rd_len=\$max_rlen\" > config.config ;
            echo -e \"[LIB]\\nq1=../../$input1\\nq2=../../$input2\" >> config.config ;
            if [ -e SOAP.fasta ]; then rm SOAP.fasta ; fi ;
            for k in $Ks_for_soap ; do
                if [ \$k -gt \$min_rlen ]; then
                    echo "WARNING: Kmer size \$k exceeds minimum read length \${min_rlen}. Please double check parameters." ;
                else
                    $soapdenovotrans pregraph -s config.config -o outputGraph_\$k -K \$k -p $threads ;
                    $soapdenovotrans contig -g outputGraph_\$k ;
                    cat outputGraph_\$k.contig | sed "s/^>/>k\${k}_/g" >> SOAP.fasta ;
                fi ;
            done ;

            cd ../../ ;
            $dedupe in=$sample_name/SOAPassembly/SOAP.fasta out=stdout.fa threads=$threads overwrite=true | \
                $fasta_formatter | \
                awk '!/^>/ { next } { getline seq } length(seq) > $min_contig_len { print \$0 "\\n" seq }' > $output ;
            if [ ! -s $output ] ; then
                rm $output ;
                echo "ERROR: de novo assembled contigs fasta file is empty." ;
                echo "Please check paths for SOAPdenovoTrans, dedupe and fasta" ;
                echo "formatter are correct, and their dependencies are installed." ;
            fi ;
            """, "assemble"
        }
    }
}

create_salmon_index = {
    def sample_name = branch.name
    def salmon_index = sample_name + "/all_fasta_index"
    if (type == "quant") {
        salmon_index = sample_name + "/salmon_quant_index"
    }
    output.dir = salmon_index
    def index_fasta = output.dir + "/" + sample_name + ".fasta"
    produce(index_fasta, '*.bin'){
        exec """
        cat $trans_fasta $input.fasta > $output1 ;
        $salmon index -t $output1 -i $salmon_index -p $threads ;
        """, "create_salmon_index"
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

    produce("eq_classes.txt*"){
        exec """
        cd $output.dir/../.. ;
        $salmon quant --dumpEq --seqBias --validateMappings --hardFilter -i $salmon_index -l A -r $rf1 $rf2 -p $threads -o $base_outdir
        """, "run_salmon"
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

run_de = {
    def run_de_bool = run_de_step.toBoolean()
    def sample_name = branch.name
    output.dir = sample_name
    produce("eq_classes_de.txt"){
        if(run_de_bool) {
            exec """
            ${R}script $code_base/DE/compare_eq_classes.R $sample_name $input $trans_fasta $output --FDR=$fdr --minCPM=$min_cpm --minLogFC=$min_logfc
            """, "run_de"
        } else {
            exec """
            $python $code_base/DE/get_novel_contigs.py $input $trans_fasta $output.dir/${sample_name}_denovo_filt.fasta
            """, "run_de"
        }
    }
}

filter_on_significant_ecs = {
    def sample_name = branch.name
    output.dir = sample_name
    produce("de_contigs.fasta"){
        exec """
        $python $code_base/util/filter_fasta.py $output.dir/${sample_name}_denovo_filt.fasta $input.txt --col_id contig > $output1 ;
        """
    }
}

align_contigs_against_genome = {
    def sample_name = branch.name
    output.dir = sample_name
    produce('aligned_contigs_against_genome.sam'){
        exec """
        $gmap -D $gmap_refdir -d $gmap_genome -f samse -t $threads -x $min_gap --max-intronlength-ends=500000 -n 0 $input.fasta > $output
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
        """, "annotate_contigs"
    }
}

refine_contigs = {
    def sample_name = branch.name
    output.dir = sample_name
    produce("novel_contigs.vcf", "novel_contigs_info.tsv", "novel_contigs.bam", "novel_contigs.fasta"){
        exec """
        $python ${code_base}/annotate/refine_annotations.py \
            $input.tsv $input.vcf $input.bam $tx_annotation \
            $genome_fasta $output.prefix \
            --minClip $min_clip \
            --minGap $min_gap \
            --mismatches $splice_motif_mismatch \
            --log $output.dir/refine.log > $output.vcf ;
        $samtools index $output.bam ;
        $python $code_base/util/filter_fasta.py $input.fasta $output.tsv --col_id contig_id > $output.fasta ;
        """
    }
}

calculate_VAF = {
    output.dir = branch.name
    produce("vaf_estimates.txt"){
        exec """
        ${R}script ${code_base}/annotate/estimate_VAF.R $branch.name/ec_count_matrix.txt $branch.name/salmon_out/quant.sf $input.tsv $trans_fasta $tx2gene $output
        """, "calculate_VAF"
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
                          calculate_VAF +
                          post_process ]
}
