code_base="/group/bioi1/marekc/20170918_cryptic_fusions/CrypticVariant/"

//Trim options
minQScore=20 //trimmomatic quality cut off
threads=8
scores=33
min_read_length=50;

//Assembly options
Ks="79 49 19" //"31 25 19"
max_read_length=150;

//software
trimmomatic="trimmomatic"
soap="/usr/bin/time -v /group/bioi1/nadiad/software/SOAP/SOAPdenovo-Trans-127mer" ;
trinity="Trinity"
fasta_dedupe="/group/bioi1/nadiad/software/bbmap/dedupe.sh" ;
fastq_dedupe="fastuniq" ;
bowtie2="/usr/bin/time -v bowtie2"
gtf2bed="gtf2bed"
bedops="bedops"
gmap="gmap"

//reference
genome_fasta="/group/bioi1/shared/genomes/hg38/fasta/hg38.fa"
trans_fasta="/group/bioi1/shared/transcriptomes/hg38/Homo_sapiens.GRCh38.cdna.all.fa"
ann_info=code_base+"/gen24_hg38.info"
ann_superTranscriptome=code_base+"/gen24_hg38.super_transcriptome.fasta"
gmap_index="/group/bioi1/shared/genomes/hg38/gmapdb"
gmap_tx_index="/group/bioi1/shared/transcriptomes/hg38/indexes/gmapdb"

controls_dir="controls"
sample_n_controls=29
bootstrap_iters=1

//Make a directory for each sample
make_sample_dir= {
   from("*.gz"){
      output.dir=branch.name
      produce(branch.name+".ignore"){
         exec """
            if [ ! -d $output.dir ]; then mkdir $output.dir ; fi ;
                  touch $output #this is just to get around the dir. being passed.
         """
       }
   }
}

dedupe = {
   from("*.gz"){
      output.dir=branch.name ;
      produce(branch.name+'.1.fastq.gz',branch.name+'.2.fastq.gz'){
         exec """
             gunzip -c $input1.gz > $branch.name/temp1.fastq ;
             gunzip -c $input2.gz > $branch.name/temp2.fastq ;
         echo $branch.name/temp1.fastq > $branch.name/fastq.list ;
         echo $branch.name/temp2.fastq >> $branch.name/fastq.list ;
         echo "Reads before:" ; wc -l $branch.name/temp1.fastq ;
         $fastq_dedupe -i $branch.name/fastq.list -o $output1.prefix -p $output2.prefix ;
         echo "Reads after:" ; wc -l $output1.prefix ;
         gzip $output1.prefix $output2.prefix ;
         rm $branch.name/fastq.list $branch.name/temp1.fastq $branch.name/temp2.fastq
      """
      }
   }
}

SOAPassemble = {
     output.dir=branch.name ;
     produce(branch.name+'_denovo_filt.fasta', branch.name+'.fasta'){
         exec """
             $trimmomatic PE -threads $threads -phred$scores $input1.gz $input2.gz
                 $branch.name/trim1.fastq /dev/null $branch.name/trim2.fastq /dev/null
                 LEADING:$minQScore TRAILING:$minQScore MINLEN:$min_read_length ;

         if [ ! -d $output.dir/SOAPassembly ]; then mkdir $output.dir/SOAPassembly ; fi ;
         cd $branch.name/SOAPassembly ;

         echo \"max_rd_len=\$max_read_length\" > config.config ;
         echo -e \"[LIB]\\nq1=../trim1.fastq\\nq2=../trim2.fastq\" >> config.config ;
         if [ -e SOAP.fasta ]; then rm SOAP.fasta ; fi ;
         for k in $Ks ; do
               $soap pregraph -s config.config -o outputGraph_\$k -K \$k -p $threads ;
           $soap contig -g outputGraph_\$k -p $threads ;
           cat outputGraph_\$k.contig | sed "s/^>/>k\${k}_/g" >> SOAP.fasta ;
         done ;
         cd ../../ ;
         $fasta_dedupe in=$branch.name/SOAPassembly/SOAP.fasta out=stdout.fa threads=$threads overwrite=true |
         fasta_formatter |
         awk '!/^>/ { next } { getline seq } length(seq) > $max_read_length { print \$0 "\\n" seq }'
         > $output1 ;
         cat $output1 $trans_fasta > $output2 ;
         rm $branch.name/trim1.fastq $branch.name/trim2.fastq ;
     ""","SOAPassemble"
     }
}

align_contigs_against_genome = {
   output.dir=branch.name
   produce('aligned_contigs_against_genome.sam'){
      exec "$gmap -D $gmap_index -d hg38 -f samse -t $threads -n 0 $input.fasta > $output.sam"
   }
}

filter_contigs_against_genome = {
   output.dir=branch.name
   produce('filtered_contigs.bam', 'interesting_contigs.txt', 'genome_filtered.fasta'){
      exec """
      python ${code_base}/filter_contigs.py $input.sam $output.bam --splice_juncs $ann_info;
      python ${code_base}/filter_fasta.py $input.fasta $output.txt > $output.fasta ;
      """
   }
}
      //cat $output.dir/genome_filtered.fasta $trans_fasta > $output.fasta

align_contigs_against_transcriptome = {
   output.dir=branch.name
   produce('aligned_contigs_against_txome.sam'){
      exec "$gmap -D $gmap_tx_index -d hg38 -f samse -t $threads -n 0 $input.fasta > $output.sam"
   }
}

filter_contigs_against_transcriptome = {
   output.dir=branch.name
   produce('filtered_contigs_against_txome.bam', 'all.groupings', 'all.fasta'){
      exec """
      python ${code_base}/filter_contigs.py $input.sam $output.bam --groupings $trans_fasta ;
      python ${code_base}/filter_fasta.py $input.fasta $output.groupings > $output.dir/transcriptome_filtered.fasta ;
      cat $output.dir/transcriptome_filtered.fasta $trans_fasta > $output.fasta ;
      """
   }
}

//blat_against_genome = {
//   output.dir=branch.name
//   produce('against_genome.psl'){
//      exec "blat $genome_fasta $input1.fasta -minIdentity=98 -minScore=30 $output"
//   }
//}
//
//filter_blat_against_genome = {
//   output.dir=branch.name
//   produce('genome_filtered.fasta'){
//      exec "${code_base}/parse_genome_blat $ann_info $input.psl $input.fasta > $output"
//   }
//}

//blat_against_transcriptome = {
//   output.dir=branch.name
//   produce('against_transcriptome.psl'){
//      exec "blat $trans_fasta $input.fasta -minIdentity=98 -minScore=30 -mask=out $output"
//   }
//}
//
//filter_blat_against_transcriptome = {
//   output.dir=branch.name
//   produce('all.groupings','all.fasta'){
//      exec """
//         ${code_base}/parse_transcriptome_blat $ann_info $input.psl $input.fasta > $output.groupings ;
//         python ${code_base}/filter_fasta.py $input.fasta $output.groupings > $output.dir/transcriptome_filtered.fasta ;
//         cat $output.dir/transcriptome_filtered.fasta $trans_fasta > $output.fasta ;
//     """
//   }
//}

create_salmon_index = {
   def salmon_index=branch.name+"/all_fasta_index"
   output.dir=branch.name+"/all_fasta_index"
   produce('bwaidx.sa'){
      exec """
         module load salmon ;
         salmon index -t $input.fasta -i $salmon_index --type fmd ;
     """
   }
}

run_salmon = {
   def workingDir = System.getProperty("user.dir");
   def (rf1, rf2)=inputs.fastq.gz.split().collect { workingDir+"/$it" }
   def salmon_index="all_fasta_index"
   def base_outdir = "salmon_out"

   if(type=="controls"){
        output.dir=branch.parent.parent.name+"/"+controls_dir+"/"+branch.name+"_salmon_out/aux_info"
        base_outdir=branch.name+"_salmon_out"
        salmon_index="../all_fasta_index"
   } else {
        output.dir=branch.parent.name+"/salmon_out/aux_info"
   }

   produce("eq_classes.txt"){
      exec """
        cd $output.dir/../.. ;
        module load salmon ;
        salmon quant --dumpEq -i $salmon_index -l A -1 $rf1 -2 $rf2 -p $threads -o $base_outdir
     """
   }
}

filter_on_significant_ecs = {
   output.dir=branch.name
   output_prefix = branch.name+"/eq_class_comp"
   salmon_dir = branch.name+"/salmon_out/aux_info"
   def sample_names=inputs.split().collect { it.split('/')[-3].split('_salmon_out')[0] }
   sample_names.set(0, branch.name) // case sample, rest are controls
   sample_names = sample_names.join(',')
   produce("ec_count_matrix.txt", "eq_class_comp_diffsplice.txt", "diffspliced_contigs.fasta", "all_filt.fasta"){
      exec """
        module load R/3.3.2 ;
        python $code_base/create_ec_count_matrix.py $inputs $sample_names $output1 ;
        Rscript $code_base/compare_eq_classes.R $output1 $output.dir/all.groupings $salmon_dir $output2 \
            --sample=$sample_n_controls --iters=$bootstrap_iters ;
        python $code_base/filter_fasta.py $input.fasta $output2.txt --col_id transcript > $output3 ;
        cat $trans_fasta $output3 > $output4
      """
   }
}

run_lace = {
   output.dir=branch.name
//       source /group/bioi1/nadiad/software/anaconda2/bin/activate lace ;
   produce('SuperDuper.fasta'){
    exec """
       module load blat ;
       if [ ! -d $output.dir/lace ]; then mkdir $output.dir/lace ; fi ;
       python /group/bioi1/nadiad/software/Lace/Lace.py
          --cores $threads
          --outputDir $output.dir/lace
          $input.fasta $input.groupings ;
       mv $output.dir/lace/SuperDuper.fasta $output ;
       rm -rf $output.dir/lace
    """
   }
}

annotate_superTranscript = {
   output.dir=branch.name
   produce("SuperDuper-Ann.gtf","SuperDuper-Ass.gtf","SuperDuper-Ann.bed",
   "SuperDuper-Ass.bed", "SuperDuper-part.bed","SuperDuper-Ann.juncs","SuperDuper.bed"){
      exec """
         module load blat ;
         module load bedops ;
         module load bedtools ;
         module load fastx-toolkit ;
         blat $input.fasta $trans_fasta -minScore=30 -minIdentity=98 $output.dir/SuperDuper-Ann.psl ;
         $code_base/psl2gtf $output.dir/SuperDuper-Ann.psl | bedtools sort | awk '{if(length(\$1)<128){print \$0}}' > $output1 ;
         blat $input.fasta $output.dir/transcriptome_filtered.fasta -minScore=30 -minIdentity=98 $output.dir/SuperDuper-Ass.psl ;
         $code_base/psl2gtf $output.dir/SuperDuper-Ass.psl | bedtools sort > $output2 ;
         $code_base/psl2sjdbFileChrStartEnd $output.dir/SuperDuper-Ann.psl > $output6 ;
         $gtf2bed < $output1 > $output3 ; $gtf2bed < $output2 > $output4 ;
         $bedops -p $output3 $output4 | cut -f1-3 | uniq - | awk '{if(\$3 - \$2 > 2){print}}' - > $output5 ;
         python $code_base/get_intersect_bed.py $output3 $output4 $output5 > $output7 ;
      """
   }
}

build_STAR_reference = {
    output.dir=branch.name+"/STARRef"
    produce("Genome"){
        exec """
            module load star ;
            cd $branch.name ;
            if [ ! -d STARRef ]; then mkdir STARRef ; fi ;
            STAR --runMode genomeGenerate --genomeDir STARRef --sjdbFileChrStartEnd ../$input.juncs
                --sjdbOverhang 99 --genomeFastaFiles ../$input.fasta  --runThreadN $threads  --genomeSAindexNbases 5 ;
        """
    }
}

map_reads = {
   def out_prefix="STAR"
   def ref="STARRef"
   def workingDir = System.getProperty("user.dir");
   def read_files=inputs.fastq.gz.split().collect { workingDir+"/$it" }.join(' ')
   if(type=="controls"){
        output.dir=branch.parent.name+"/"+controls_dir
        out_prefix=branch.name
        ref="../"+ref
   } else {
         output.dir=branch.name
   }
   produce(out_prefix+"Aligned.sortedByCoord.out.bam",out_prefix+"SJ.out.tab"){
    exec """
        module load star ;
        module load samtools ;
        cd $output.dir ;
        time STAR --genomeDir $ref --readFilesCommand zcat
           --readFilesIn $read_files
           --outSAMtype BAM SortedByCoordinate --outFileNamePrefix $out_prefix
           --runThreadN $threads --limitBAMsortRAM 5050000000 ;
        samtools index ${workingDir}/${output1} ;
        rm -rf ${out_prefix}_STARtmp ;
    """
   }
}

map_reads_controls = {
    def out_prefix="STAR"
    def ref="../STARRef"
    def workingDir = System.getProperty("user.dir");

    def reads = inputs.fastq.gz.split().collect { workingDir+"/$it" }
    def (reads_r1, reads_r2) = [[],[]]
    reads.eachWithIndex { readfile, idx ->
        if (idx % 2 == 0)
            reads_r1.add(readfile)
        else
            reads_r2.add(readfile)
    }
    def readPairs = [reads_r1, reads_r2].transpose().collect { it.join(" ") }

    int n = reads.size / 2
    def outPrefixes = (1..n).collect { "control${it}." }
    def outBams = outPrefixes.collect { "${it}Aligned.sortedByCoord.out.bam" }
    def outTabs = outPrefixes.collect { "${it}SJ.out.tab" }

    output.dir = branch.parent.name+"/"+controls_dir
    def args = [readPairs, outPrefixes, outBams].transpose()
    def commands = args.collect {
        "cd $output.dir; time STAR --genomeDir $ref --readFilesCommand zcat" +
        " --readFilesIn " + it[0] +
        " --outSAMtype BAM SortedByCoordinate --outFileNamePrefix " + it[1] +
        " --runThreadN $threads --limitBAMsortRAM 3050000000 ; " +
        "samtools index " + it[2] + " ; " +
        "rm -rf " + it[1] + "_STARtmp ;"
    }

    produce(outBams, outTabs){
        multiExec commands
    }
}

get_info_on_novel_events = {
    def out_prefix=""
    output.dir=branch.name
    bed_dir=branch.name
    if(type=="controls"){
        output.dir=branch.parent.name+"/"+controls_dir
        out_prefix=branch.name+"."
        bed_dir=branch.parent.name
    }
    produce(out_prefix+"novel.junctions",out_prefix+"novel.blocks"){
    exec """
        module load samtools ;
        awk '\$6 == \"0\" { print \$0 }' $input.tab > $output1 ;
        rm -rf $output2 ;
        cat ${bed_dir}/SuperDuper.bed | while read line ; do
            region=`echo $line | awk '{ printf "%s:%s-%s\\n", \$1, \$2 + 1, \$3}'` ;
            ave=`samtools depth $input.bam -r \$region | awk '{sum+=\$3;} END { if (NR > 0){print sum/NR} else { print "0"} }'`;
            echo -e "\$line\\t\$ave" >> $output2 ;
        done ;
    """
    }
}

get_filtered_variants = {
    output.dir=branch.name
    produce("novel.summary"){
    exec """
        cd $output.dir ;
        $code_base/parse_superTranscript_results novel.blocks novel.junctions all.groupings controls/*.novel.*  > ../$output
    """
    }
}

//fastqInputFormat="%_L001_R*.fastq.gz"
fastqInputFormat="%_R*.fastq.gz"

run { fastqInputFormat * [ make_sample_dir +
                        dedupe +
                        SOAPassemble +
              align_contigs_against_genome +
              filter_contigs_against_genome +
//              blat_against_genome +
//              filter_blat_against_genome +
              align_contigs_against_transcriptome +
              filter_contigs_against_transcriptome +
              create_salmon_index +
              [run_salmon, "controls/%.*.fastq.gz" * [ run_salmon.using(type:"controls") ]] +
              filter_on_significant_ecs +
              run_lace +
              annotate_superTranscript +
              build_STAR_reference +
              map_reads
//              get_info_on_novel_events +
//              "controls/%.*.fastq.gz" *  [ map_reads.using(type:"controls") +
//                            get_info_on_novel_events.using(type:"controls") ] +
//              get_filtered_variants
              ]
}
