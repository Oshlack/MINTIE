code_base="/group/bioi1/marekc/20170918_cryptic_variant/CrypticVariant/"

//Trim options
minQScore=20 //trimmomatic quality cut off
threads=8
scores=33
min_read_length=50;
genome_mem=31000000000
sort_ram=4 //ram per code for bam sorting

//Assembly options
Ks="79 49 19" //"31 25 19"
max_read_length=150;

//software
trimmomatic="trimmomatic"
soap="/usr/bin/time -v /group/bioi1/nadiad/software/SOAP/SOAPdenovo-Trans-127mer" ;
trinity="Trinity"
fasta_dedupe="/group/bioi1/nadiad/software/bbmap/dedupe.sh" ;
fastq_dedupe="/group/bioi1/marekc/apps/FastUniq/source/fastuniq" ;
bowtie2="/usr/bin/time -v bowtie2"
gtf2bed="gtf2bed"
bedops="bedops"
gmap="/group/bioi1/marekc/apps/GMAP-GSNAP/src/gmap"
gmap_build="perl /group/bioi1/marekc/apps/GMAP-GSNAP/util/gmap_build.pl -B=/group/bioi1/marekc/apps/GMAP-GSNAP/src"
salmon="/group/bioi1/marekc/apps/Salmon-latest_linux_x86_64/bin/salmon"
hisat="/group/bioi1/marekc/apps/hisat-0.1.6-beta/hisat"
hisat_build="/group/bioi1/marekc/apps/hisat-0.1.6-beta/hisat-build"

//reference
genome_fasta="/group/bioi1/shared/genomes/hg38/fasta/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
trans_fasta="/group/bioi1/shared/transcriptomes/hg38/Homo_sapiens.GRCh38.cdna.all.fa"
//ann_info=code_base+"/gen24_hg38.info"
ann_superTranscriptome=code_base+"/gen24_hg38.super_transcriptome.fasta"
gmap_index="/group/bioi1/shared/genomes/hg38/gmapdb"
gmap_tx_index="/group/bioi1/shared/transcriptomes/hg38/indexes/gmapdb"
ann_info="/group/bioi1/shared/genomes/hg38/gtf/gencode.v24.annotation.gtf.info"
tx_annotation="/group/bioi1/shared/genomes/hg38/gtf/gencode.v24.annotation.gtf.gz"

controls_dir="controls"
sample_n_controls=29
bootstrap_iters=1

//Make a directory for each sample
make_sample_dir= {
    output.dir=branch.name
    from("*.gz"){
      output.dir=branch.name
      produce(branch.name+".ignore"){
         exec """
            mkdir -p $output.dir ;
            touch $output
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
      """, "dedupe"
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
           $soap contig -g outputGraph_\$k ;
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

filter_contigs_against_genome = {
   output.dir=branch.name
   produce('filtered_contigs_against_genome.bam', 'novel_contigs.txt', 'genome_filtered.fasta'){
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
      exec """
        $gmap -D $gmap_tx_index -d hg38 -f samse -t $threads -n 0 $input.fasta > $output.sam
      """, "align_contigs_against_transcriptome"
   }
}

filter_contigs_against_transcriptome = {
   output.dir=branch.name
   produce('filtered_contigs_against_txome.bam', 'all.groupings', 'all.fasta'){
      exec """
      python ${code_base}/filter_contigs.py $input.sam $output.bam --groupings $trans_fasta ;
      python ${code_base}/filter_fasta.py $input.fasta $output.groupings > $output.dir/filtered_contigs_against_txome.fasta ;
      cat $output.dir/filtered_contigs_against_txome.fasta $trans_fasta > $output.fasta ;
      """
   }
}

create_salmon_index = {
   def salmon_index=branch.name+"/all_fasta_index"
   output.dir=branch.name+"/all_fasta_index"
   produce('sa.bin','hash.bin','rsd.bin'){
      exec """
         $salmon index -t $input2 -i $salmon_index ;
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
        $salmon quant --dumpEq --seqBias -i $salmon_index -l A -r $rf1 $rf2 -p $threads -o $base_outdir
     """, "run_salmon"
   }
}

create_ec_count_matrix = {
   output.dir = branch.name
   def sample_names=inputs.split().collect { it.split('/')[-3].split('_salmon_out')[0] }
   sample_names.set(0, branch.name) // case sample, rest are controls
   sample_names = sample_names.join(',')
   produce("ec_count_matrix.txt"){
      exec """
        python $code_base/create_ec_count_matrix.py $inputs $sample_names $output1 ;
      """, "create_ec_count_matrix"
   }
}

run_de = {
   output.dir = branch.name
   def case_name = branch.name
   //def salmon_dir = branch.name+"/salmon_out/aux_info"
   produce("eq_class_comp_diffsplice.txt"){
      exec """
        Rscript $code_base/compare_eq_classes.R $case_name $input $trans_fasta $output ;
      """, "run_de"
   }
}

filter_on_significant_ecs = {
   output.dir=branch.name
   produce("de_contigs.fasta"){
      exec """
        python $code_base/filter_fasta.py $input.fasta $input.txt --col_id contig > $output1 ;
      """
   }
}

align_contigs_against_genome = {
   output.dir=branch.name
   produce('aligned_contigs_against_genome.bam'){
      exec """
        $gmap -D $gmap_index -d hg38 -f samse -t $threads -n 0 $input.fasta > $output.dir/tmp.sam ;
        samtools sort -o $output $output.dir/tmp.sam ;
        samtools index $output ;
        rm $output.dir/tmp.sam ;
      """, "align_contigs_against_genome"
   }
}

annotate_contigs = {
   output.dir=branch.name
   def sample_name = inputs.fastq.gz.split()[0].split('/').last().split('\\.').first().split('_R').first()
   produce("novel_contigs.vcf", "novel_contigs_info.tsv", "novel_contigs.bam"){
      exec """
        time python ${code_base}/annotate/annotate_contigs.py \
            $sample_name $input.bam $output.bam $output.tsv $ann_info \
            $tx_annotation --log $output.dir/annotate.log> $output.vcf
      """
   }
}

create_supertranscript_reference = {
   output.dir=branch.name
   produce("tx_annotation.gtf", "supertranscript.fasta"){
      exec """
          echo "extracting relevant transcripts to gtf reference..." ;
          cat $input1 | cut -f 1 | sed 1d | sort | uniq | awk '{split(\$0, x, "|")}{print x[1]"\\n"x[2]}' | \
                sort | uniq | sed '/^ *\$/d' > $output.dir/gene_list.txt ;
          zcat < $tx_annotation | fgrep -wf $output.dir/gene_list.txt > $output1 ;
          echo "generating supertranscript blocks..." ;
          python ${code_base}/generate_st_blocks.py $output1 ${output.dir}/st_blocks.bed ;
          echo "extracting fasta sequence..." ;
          bedtools getfasta -fi $genome_fasta -bed ${output.dir}/st_blocks.bed -fo ${output.dir}/st_blocks.fasta -s ;
          echo "making supertranscripts..." ;
          python ${code_base}/make_supertranscript_ref.py $input.fasta \
                $input1 ${output.dir}/st_blocks.bed ${output.dir}/st_blocks.fasta $output2 ;
      """
   }
}

make_super_supertranscript = {
    colpath = inputs.fastq.gz.split()
    colpath = colpath[(0..(colpath.length-1)).step(2)]
    colpath = colpath.collect { it.split('/').last().split('_').first() }.join('_')
    colpath = colpath+'_collated_output'
    output.dir = colpath
    produce('novel_contigs_annotated.txt', 'supersupertranscript.fasta'){
        exec """
            cat $inputs.fasta > $output.fasta ;
            file1=`echo "$inputs.txt" | cut -f1 -d' '` ;
            head -1 \$file1 > $output.txt ;
            for file in $inputs.txt; do sed 1d $file >> $output.txt ; done ;
        """
    }
}

annotate_supertranscript = {
   //clinker_out=branch.name+"/clinker_out"
   clinker_out=colpath+'/clinker'
   output.dir=clinker_out+"/reference"
   produce("fst_reference.fasta"){
      exec """
         python ${code_base}/Clinker/main.py -in $input.txt -out $clinker_out -pos 4,5,6,7,21 \
            -del t -header true -competitive true -st $input.fasta ;
         faidx $output ;
      """
   }
}

make_supertranscript_gmap_reference = {
   //clinker_out=branch.name+"/clinker_out"
   output.dir=colpath+'/clinker/reference'
   produce("st_gmap_ref"){
      exec """
         $gmap_build -s chrom -k 15 -d st_gmap_ref -D $output.dir $input.fasta ;
      """
   }
}

align_contigs_to_supertranscript = {
   output.dir=colpath+"/clinker/alignment"
   index_dir=colpath+"/clinker/reference"
   //TODO: this needs to play nicely with the given mask for cases
   def sample_name=inputs.fastq.gz.split()[0].split('/').last().split('\\.').first().split('_R').first()
   produce(sample_name+"_novel_contigs_st_aligned.bam"){
      exec """
         outfile=$output ; basename="\${outfile%.*}" ;
         $gmap -D $index_dir -d st_gmap_ref -f samse -t $threads -n 0 ${sample_name}/diffspliced_contigs.fasta > \${basename}.sam ;
         samtools view -hb \${basename}.sam > \${basename}_unsort.bam ;
         samtools sort \${basename}_unsort.bam > $output ;
         samtools index $output ; rm \${basename}_unsort.bam ; rm \${basename}.sam
      """, "align_contigs_to_supertranscript"
   }
}

star_genome_gen = {
    doc "Generate STAR genome index"
    output.dir=colpath+"/clinker/"
    genome_folder = output.dir+"/genome"
    produce("$genome_folder/Genome") {
        exec """STAR --runMode genomeGenerate
                --runThreadN $threads
                --genomeDir $genome_folder
                --genomeFastaFiles $input.fasta
                --limitGenomeGenerateRAM $genome_mem
                --genomeSAindexNbases 5""", "star_genome_gen"
    }
}

star_align = {
   //output.dir=branch.parent.name+"/clinker_out/alignment/"
   output.dir=colpath+"/clinker/alignment/"
   def workingDir=System.getProperty("user.dir");
   def read_files=inputs.fastq.gz.split().collect { workingDir+"/$it" }.join(' ')
   def sample_name=read_files.split()[0].split('/').last().split('\\.').first()
   def out_prefix=output.dir+'/'+sample_name
   if(type=="controls"){
        //output.dir=branch.parent.parent.name+"/clinker_out/alignment/controls/"
        output.dir=colpath+"/clinker/alignment/controls/"
        out_prefix=output.dir+'/'+sample_name
   }
   produce(out_prefix+"Aligned.sortedByCoord.out.bam",out_prefix+"SJ.out.tab"){
        //mkdir -p $output.dir ;
        exec """
        rm -rf ${out_prefix}_STARtmp ;
        time STAR --genomeDir $genome_folder
           --readFilesCommand zcat
           --readFilesIn $read_files
           --outSAMtype BAM SortedByCoordinate
           --outFileNamePrefix $out_prefix
           --runThreadN $threads
           --limitBAMsortRAM $genome_mem
           --genomeSAindexNbases 5
           --outWigStrand Unstranded
           --outWigType bedGraph
           --outWigNorm RPM ;
        samtools index ${workingDir}/${output1} ;
    """, "star_align"
   }
}

hisat_index = {
    output.dir=colpath+"/clinker/"
    genome_folder = output.dir+"/genome"
    def idx_prefix = genome_folder+"/hisat_index"
    produce(genome_folder+"/hisat_index.rev.1.bt2") {
        exec """
        $hisat_build $input.fasta $idx_prefix
        """
    }
}

hisat_align = {
   output.dir=colpath+"/clinker/alignment/"
   def workingDir=System.getProperty("user.dir");
   def (r1, r2)=inputs.fastq.gz.split().collect { workingDir+"/$it" }
   def sample_name=r1.split('/').last().split('\\.').first()
   def out_prefix=output.dir+'/'+sample_name
   if(type=="controls"){
        output.dir=colpath+"/clinker/alignment/controls/"
        out_prefix=output.dir+'/'+sample_name
   }
   produce(out_prefix+"_hisatAligned.bam"){
        exec """
        outfile=$output.bam ; basename="\${outfile%.*}" ;
        hisat_idx=$input.bt2; idx=\${hisat_idx%.rev.1.bt2} ;
        time $hisat --threads $threads -x \$idx -1 $r1 -2 $r2 > \${basename}.sam ;
        samtools sort -@ $threads -m ${sort_ram}G \${basename}.sam -o $output ;
        samtools index $output ;
        if [[ -f $output ]] ; then
            rm -f \${basename}.sam ;
        fi
    """, "hisat_align"
   }
}

if(!binding.variables.containsKey("fastqCaseFormat")){
    fastqCaseFormat="cases/%_R*.fastq.gz"
}
if(!binding.variables.containsKey("fastqControlFormat")){
    fastqControlFormat="controls/%_R*.fastq.gz"
}

run { fastqCaseFormat * [ make_sample_dir +
                          dedupe +
                          SOAPassemble +
//                          filter_contigs_against_genome +
//                          align_contigs_against_transcriptome +
//                          filter_contigs_against_transcriptome +
                          create_salmon_index +
                          [run_salmon, fastqControlFormat * [ run_salmon.using(type:"controls") ]] +
                          create_ec_count_matrix +
                          run_de +
                          filter_on_significant_ecs +
                          align_contigs_against_genome +
                          annotate_contigs +
                          create_supertranscript_reference ] +
       make_super_supertranscript +
       annotate_supertranscript +
       make_supertranscript_gmap_reference +
       star_genome_gen +
       hisat_index +
       fastqCaseFormat * [ align_contigs_to_supertranscript, hisat_align] +
       fastqControlFormat * [ hisat_align.using(type:"controls") ]
}
