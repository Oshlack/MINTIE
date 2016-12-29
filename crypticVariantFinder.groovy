code_base="/mnt/storage/nadiad/work_area/20160203_ALL/new_new_code"

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
soap="/usr/bin/time -v /mnt/storage/nadiad/software/SOAP/SOAPdenovo-Trans-127mer" ;
fasta_dedupe="/mnt/storage/nadiad/software/bbmap/dedupe.sh" ;
fastq_dedupe="fastuniq" ;
bowtie2="/usr/bin/time -v bowtie2"

//reference
genome_fasta="/mnt/storage/shared/genomes/hg38/fasta/hg38.fa"
trans_fasta=code_base+"/Homo_sapiens.GRCh38.cdna.all.fa"
ann_info=code_base+"/gen24_hg38.info"
ann_superTranscriptome=code_base+"/gen24_hg38.super_transcriptome.fasta"

/**ref_fasta="/mnt/storage/nadiad/work_area/20160203_ALL/new_code/super_reference.fasta"
ref_coding_blocks="/mnt/storage/nadiad/work_area/20160203_ALL/new_code/super_reference.coding_blocks"
ref_tab="/mnt/storage/nadiad/work_area/20160203_ALL/code/hg38_genCode20.tab"
gene_list=code_base+"/ALL_genes_of_interest" **/

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
     produce(branch.name+'.fasta'){
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
	     > $output ;
	     rm $branch.name/trim1.fastq $branch.name/trim2.fastq ;
	 """ 
     }
}

blat_against_genome = {
   output.dir=branch.name
   produce('against_genome.psl'){
      exec "blat $genome_fasta $input -minIdentity=98 -minScore=100 $output"
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

run_lace = {
   output.dir=branch.name
   produce('SuperDuper.fasta'){
	exec """
	   if [ ! -d $output.dir/lace ]; then mkdir $output.dir/lace ; fi ; 
	   source activate lace ;
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
   produce("SuperDuper-Ann.gtf","SuperDuper-Ass.gtf",
	   "SuperDuper-Ann.juncs","SuperDuper.bed"){
      exec """
         blat $input.fasta $trans_fasta -minScore=100 -minIdentity=98 $output.dir/SuperDuper-Ann.psl ;
	 $code_base/psl2gtf $output.dir/SuperDuper-Ann.psl | bedtools sort > $output1 ;
         blat $input.fasta $output.dir/genome_filtered.fasta -minScore=100 -minIdentity=98 $output.dir/SuperDuper-Ass.psl ;
	 $code_base/psl2gtf $output.dir/SuperDuper-Ass.psl | bedtools sort > $output2 ;
	 $code_base/psl2sjdbFileChrStartEnd $output.dir/SuperDuper-Ann.psl > $output3 ;
	 bedtools multiinter -names Annotation Assembly -i $output1 $output2 | cut -f1-3,5 > $output4 ;
      """
   }
}

map_reads = {
   output.dir=branch.name
   produce("STARAligned.sortedByCoord.out.bam"){
      exec """
      cd $branch.name ;
      if [ ! -d STARRef ]; then mkdir STARRef ; fi ;
      STAR --runMode genomeGenerate --genomeDir STARRef --sjdbFileChrStartEnd ../$input.juncs 
      	   --sjdbOverhang 99 --genomeFastaFiles ../$input.fasta  --runThreadN $threads  --genomeSAindexNbases 5 ;
      F1=`echo $inputs.fastq.gz | cut -f1 -d \" \"` ;
      F2=`echo $inputs.fastq.gz | cut -f2 -d \" \"` ;
      STAR --genomeDir STARRef --readFilesCommand zcat 
         --readFilesIn ../\$F1 ../\$F2
         --outSAMtype BAM SortedByCoordinate --outFileNamePrefix STAR
	 --runThreadN $threads ;
      samtools index ../$output ;
      rm -rf STAR_STARtmp ;
      """
   }
}

map_controls = {
   output.dir=branch.parent.name+"/controls"
   produce(branch.name+"Aligned.sortedByCoord.out.bam"){
      exec """
      if [ ! -d $output.dir ]; then mkdir $output.dir ; fi ;
      cd $output.dir ;
      STAR --genomeDir ../STARRef --readFilesCommand zcat
         --readFilesIn ../../$input1 ../../$input2
         --outSAMtype BAM SortedByCoordinate --outFileNamePrefix $branch.name
         --runThreadN $threads ;
      cd ../../ ;
      samtools index $output ;
      """
      }
}

get_info_on_novel_events = {
      output.dir=branch.name
      produce("novel.junctions","novel.blocks","novel.summary"){
      exec """
      	 awk '\$6 == \"0\" { print \$0 }' $output.dir/STARSJ.out.tab > $output1 ;
	 rm -rf $output2 ;
	 cat $input.bed | while read line ; do 
	    region=`echo $line | awk '{ printf "%s:%s-%s\\n", \$1, \$2 + 1, \$3}'` ;
	    ave=`samtools depth $input.bam -r \$region | awk '{sum+=\$3;} END { if (NR > 0){print sum/NR} else { print "0"} }'`;
	    echo -e "\$line\\t\$ave" >> $output2 ;
	 done ;
	 $code_base/parse_superTranscript_results $output1 $output2 $output.dir/all.groupings > $output3
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
			  run_lace + 
			  annotate_superTranscript +
			  map_reads +
			  "controls/%_*.fastq.gz" *  [ map_controls ] +
			  get_info_on_novel_events 
			  ]
}
