code_base="/group/bioi1/marekc/cryptic_fusions/CrypticVariant"

//Trim options
minQScore=20 //trimmomatic quality cut off
threads=16
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

//reference
genome_fasta="/group/bioi1/shared/genomes/hg38/fasta/hg38.fa"
trans_fasta=code_base+"/Homo_sapiens.GRCh38.cdna.all.fa"
ann_info=code_base+"/gen24_hg38.info"
ann_superTranscriptome=code_base+"/gen24_hg38.super_transcriptome.fasta"

controls_dir="controls"

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

blat_against_genome = {
   output.dir=branch.name
   produce('against_genome.psl'){
      exec "blat $genome_fasta $input1.fasta -minIdentity=98 -minScore=100 $output"
   }
}

filter_blat_against_genome = {
   output.dir=branch.name
   produce('genome_filtered.fasta'){
      exec "${code_base}/parse_genome_blat $ann_info $input.psl $input.fasta > $output"
   }
}

blat_against_transcriptome = {
   output.dir=branch.name
   produce('against_transcriptome.psl'){
      exec "blat $trans_fasta $input -minIdentity=98 -minScore=100 $output"
   }
}

filter_blat_against_transcriptome = {
   output.dir=branch.name
   produce('all.groupings','all.fasta'){
      exec """
         ${code_base}/parse_transcriptome_blat $ann_info $input.psl $input.fasta > $output.groupings ;
         cat $input.fasta $ann_superTranscriptome > $output.fasta ;
     """
   }
}

create_salmon_index = {
   def salmon_index=branch.name+"/all_fasta_index"
   output.dir=branch.name+"/all_fasta_index"
   produce('bwaidx.sa'){
      exec """
         salmon index -t $input.fasta -i $salmon_index --type fmd ;
     """
   }
}

run_salmon = {
   def workingDir = System.getProperty("user.dir");
   def (rf1, rf2)=inputs.fastq.gz.split().collect { workingDir+"/$it" }
   def salmon_index="all_fasta_index"
   output.dir=branch.parent.name+"/salmon_out/aux_info"
   produce("eq_classes.txt"){
      exec """
        cd $output.dir/../.. ;
        salmon quant --dumpEq -i $salmon_index -l A -1 $rf1 -2 $rf2 -o salmon_out
     """
   }
}

run_salmon_controls = {
    def workingDir = System.getProperty("user.dir");
    def reads = inputs.fastq.gz.split().collect { workingDir+"/$it" }

    //contruct output dirs
    int n = reads.size / 2
    def outdirs = (1..n).collect { branch.parent.parent.name + "/" + controls_dir + "${it}_salmon_out" }
    def outfiles = outdirs.collect { "${it}/aux_info/eq_classes.txt" }

    //pair up reads
    def (reads_r1, reads_r2) = [[],[]]
    reads.eachWithIndex { readfile, idx ->
        if (idx % 2 == 0)
            reads_r1.add(readfile)
        else
            reads_r2.add(readfile)
    }
    def readPairs = [reads_r1, reads_r2, outdirs].transpose()

    //construct salmon command args
    def salmon_index = branch.parent.parent.name + "/all_fasta_index"
    def run_salmon = "salmon quant --dumpEq -i $salmon_index -l A "
    def commands_salmon = readPairs.collect { run_salmon+"-1 "+it[0]+" -2 "+it[1]+" -o "+it[2] }

    produce(outfiles) {
        multiExec commands_salmon
    }
}

filter_on_significant_ecs = {
   output.dir=branch.name
   output_prefix = branch.name+"/eq_class_comp"
   produce("eq_class_comp_de.txt", "eq_class_comp_diffsplice.txt", "filtered_all_fasta.fasta"){
      exec """
        Rscript $code_base/compare_eq_classes.R $inputs $output.dir/all.groupings $output_prefix ;
        python $code_base/filter_fasta.py $input.fasta $output.txt | sed --expression='/^\$/d' - > $output.fasta ;
      """
   }
}
        //python $code_base/filter_fasta.py $input1.fasta $output.txt | sed --expression='/^\$/d' - > $output2 ;
        //python $code_base/filter_fasta.py $input2.fasta $output.txt | sed --expression='/^\$/d' - > $output3 ;

run_lace = {
   output.dir=branch.name
//       source /group/bioi1/nadiad/software/anaconda2/bin/activate lace ;
   produce('SuperDuper.fasta'){
    exec """
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
         blat $input.fasta $trans_fasta -minScore=100 -minIdentity=98 $output.dir/SuperDuper-Ann.psl ;
         $code_base/psl2gtf $output.dir/SuperDuper-Ann.psl | bedtools sort > $output1 ;
         blat $input.fasta $output.dir/genome_filtered.fasta -minScore=100 -minIdentity=98 $output.dir/SuperDuper-Ass.psl ;
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
        cd $output.dir ;
        time STAR --genomeDir $ref --readFilesCommand zcat
           --readFilesIn $read_files
           --outSAMtype BAM SortedByCoordinate --outFileNamePrefix $out_prefix
           --runThreadN $threads --limitBAMsortRAM 3050000000 ;
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

fastqInputFormat="%_L001_R*.fastq.gz"

run { fastqInputFormat * [ make_sample_dir +
                        dedupe +
                        SOAPassemble +
              blat_against_genome +
              filter_blat_against_genome +
              blat_against_transcriptome +
              filter_blat_against_transcriptome +
              create_salmon_index +
              [run_salmon, "controls/%.*.fastq.gz" * [ run_salmon_controls ]] +
              filter_on_significant_ecs
//              run_lace +
//              annotate_superTranscript +
//              build_STAR_reference +
//              map_reads +
//              //get_info_on_novel_events +
//              "controls/%.*.fastq.gz" *  [ map_reads_controls ]
//              //              get_info_on_novel_events.using(type:"controls") ]
//              //get_filtered_variants
              ]
}
