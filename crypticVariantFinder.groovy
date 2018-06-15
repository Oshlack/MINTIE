code_base="/group/bioi1/marekc/20170918_cryptic_variant/CrypticVariant/"

//Trim options
minQScore=20 //trimmomatic quality cut off
threads=8
scores=33
min_read_length=50;
genome_mem=5050000000

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
      """
      }
   }
}

SOAPassemble = {
     output.dir=branch.name ;
     produce(branch.name+'_denovo_filt.fasta', branch.name+'.fasta'){
         exec """
             module load trimmomatic ;
             module load fastx-toolkit ;
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

align_contigs_against_genome = {
   output.dir=branch.name
   produce('aligned_contigs_against_genome.sam'){
      exec "$gmap -D $gmap_index -d hg38 -f samse -t $threads -n 0 $input.fasta > $output.sam"
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
      exec "$gmap -D $gmap_tx_index -d hg38 -f samse -t $threads -n 0 $input.fasta > $output.sam"
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
         $salmon index -t $input.fasta -i $salmon_index ;
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
        $salmon quant --dumpEq --gcBias --seqBias -i $salmon_index -l A -r $rf1 $rf2 -p $threads -o $base_outdir
     """
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
      """
   }
}

run_diffsplice = {
   case_name = branch.name
   output.dir = branch.name
   salmon_dir = branch.name+"/salmon_out/aux_info"
   produce("eq_class_comp_diffsplice.txt"){
      exec """
        module load R/3.3.2 ;
        Rscript $code_base/compare_eq_classes.R $case_name $input $output.dir/all.groupings \
            $trans_fasta $salmon_dir $output ;
      """
   }
}

filter_on_significant_ecs = {
   output.dir=branch.name
   produce("diffspliced_contigs.fasta", "all_filt.fasta", "ds_novel_contigs.txt"){
      exec """
        python $code_base/filter_fasta.py $input.fasta $input.txt --col_id transcript > $output1 ;
        cat $trans_fasta $output1 > $output2 ;
        grep -ohE "k[0-9]+_[0-9]+" $input.txt | sort | uniq > $output3 ;
      """
   }
}

annotate_diffspliced_contigs = {
   output.dir=branch.name
   ds_results = branch.name+"/eq_class_comp_diffsplice.txt"
   produce("novel_contigs_annotated.txt", "novel_contigs.bam"){
      exec """
        module load samtools ;
        samtools view -H $output.dir/filtered_contigs_against_genome.bam > $output.dir/tmp.sam ;
        samtools view $output.dir/filtered_contigs_against_genome.bam | fgrep -w -f $input3 >> $output.dir/tmp.sam ;
        python ${code_base}/filter_contigs.py $output.dir/tmp.sam $output.bam --splice_juncs $ann_info --annotate $ds_results ;
        rm $output.dir/tmp.sam ;
      """
   }
}

create_supertranscript_reference = {
   output.dir=branch.name
   produce("tx_annotation.gtf", "st_blocks.bed", "st_blocks.fasta", "supertranscript.fasta"){
      exec """
          module load bedtools ;
          echo "extracting relevant transcripts to gtf reference..." ;
          cat $input1 | cut -f 1 | sed 1d | sort | uniq | awk '{split(\$0, x, "|")}{print x[1]"\\n"x[2]}' | sort | uniq | sed '/^ *\$/d' > $output.dir/gene_list.txt ;
          zcat < $tx_annotation | fgrep -wf $output.dir/gene_list.txt > $output1 ;
          echo "generating supertranscript blocks..." ;
          python ${code_base}/generate_st_blocks.py $output1 $output2 ;
          echo "extracting fasta sequence..." ;
          bedtools getfasta -fi $genome_fasta -bed $output2 -fo $output3 -s ;
          echo "making supertranscripts..." ;
          python ${code_base}/make_supertranscript_ref.py $input.fasta $input1 $output2 $output3 $output4 ;
      """
   }
}

annotate_supertranscript = {
   clinker_out=branch.name+"/clinker_out"
   output.dir=clinker_out+"/reference"
   produce("fst_reference.fasta"){
      exec """
         python ${code_base}/Clinker/main.py -in $input.txt -out $clinker_out -pos 4,5,6,7 -del t -header true -competitive false -st $input4 ;
      """
   }
}

make_supertranscript_gmap_reference = {
   clinker_out=branch.name+"/clinker_out"
   output.dir=branch.name
   produce("st_gmap_ref"){
      exec """
         $gmap_build -s chrom -k 15 -d st_gmap_ref -D $branch.name $input.fasta ;
      """
   }
}

align_contigs_to_supertranscript = {
   clinker_out=branch.name+"/clinker_out"
   output.dir=branch.name
   produce("novel_contigs_st_aligned.bam"){
      exec """
         outfile=$output ; basename="\${outfile%.*}" ;
         $gmap -D $output.dir -d st_gmap_ref -f samse -t $threads -n 0 ${branch.name}/diffspliced_contigs.fasta > \${basename}.sam ;
         module load samtools;
         samtools view -hb \${basename}.sam > \${basename}_unsort.bam ;
         samtools sort \${basename}_unsort.bam > $output ;
         samtools index $output ; rm \${basename}_unsort.bam ; rm \${basename}.sam
      """
   }
}

star_genome_gen = {
    doc "Generate STAR genome index"

    output.dir = branch.name+"/clinker_out"
    genome_folder = output.dir+"/genome"

    // Generate Fusion SuperTranscriptome Genome for STAR
    produce("$genome_folder/Genome") {
        exec """module load star ;
                STAR --runMode genomeGenerate
                --runThreadN $threads
                --genomeDir $genome_folder
                --genomeFastaFiles $input.fasta
                --limitGenomeGenerateRAM $genome_mem
                --genomeSAindexNbases 5""","stargen"
    }
}

//star_align = {
//    doc "Map paired-end reads using the STAR aligner: 1st pass"
//
//    //Map paired-end reads using the STAR aligner: 1st pass
//    from("*.fastq.gz") {
//        transform("(.*)_R1.fastq.gz","(.*)_R2.fastq.gz"){
//
//            // Setup stage
//            files = inputs.toString()
//            output.dir = branch.name+"/clinker_out/alignment"
//            String bam = "$output.dir/Aligned.sortedByCoord.out.bam"
//
//            // Align to fusion SuperTranscriptome
//            produce("$bam"){
//                exec """module load star ; STAR --genomeDir $genome_folder
//                    --readFilesIn ${files}
//                    --readFilesCommand gunzip -c
//                    --outSAMtype BAM SortedByCoordinate
//                    --runThreadN $threads
//                    --outFileNamePrefix "$output.dir/"
//                    --limitBAMsortRAM $genome_mem
//                    --outWigType bedGraph
//                    --outWigNorm RPM
//                    --genomeSAindexNbases 5
//                    --outWigStrand Unstranded ;
//                    module load samtools ;
//                    samtools index $bam
//                """, "star1pass"
//            }
//        }
//    }
//}
//
//run_lace = {
//   output.dir=branch.name
////       source /group/bioi1/nadiad/software/anaconda2/bin/activate lace ;
//   produce('SuperDuper.fasta'){
//    exec """
//       module load blat ;
//       if [ ! -d $output.dir/lace ]; then mkdir $output.dir/lace ; fi ;
//       python /group/bioi1/nadiad/software/Lace/Lace.py
//          --cores $threads
//          --outputDir $output.dir/lace
//          $input.fasta $input.groupings ;
//       mv $output.dir/lace/SuperDuper.fasta $output ;
//       rm -rf $output.dir/lace
//    """
//   }
//}
//
//annotate_superTranscript = {
//   output.dir=branch.name
//   produce("SuperDuper-Ann.gtf","SuperDuper-Ass.gtf","SuperDuper-Ann.juncs"){
//      exec """
//         module load blat ;
//         module load bedops ;
//         module load bedtools ;
//         module load fastx-toolkit ;
//         blat $input.fasta $trans_fasta -minScore=30 -minIdentity=98 $output.dir/SuperDuper-Ann.psl ;
//         $code_base/psl2gtf $output.dir/SuperDuper-Ann.psl | bedtools sort | awk '{if(length(\$1)<128){print \$0}}' > $output1 ;
//         blat $input.fasta $output.dir/filtered_contigs_against_txome.fasta -minScore=30 -minIdentity=98 $output.dir/SuperDuper-Ass.psl ;
//         $code_base/psl2gtf $output.dir/SuperDuper-Ass.psl | bedtools sort > $output2 ;
//         $code_base/psl2sjdbFileChrStartEnd $output.dir/SuperDuper-Ann.psl > $output3 ;
//      """
//   }
//}
//         //blat $input.fasta $genome_fasta -minScore=30 -minIdentity=98 $output.dir/SuperDuper_against_genome.psl ;
//         //python annotate_supertranscript.py $output.dir/SuperDuper_against_genome.psl $gtf_tx $output4;
//         //$gtf2bed < $output1 > $output3 ; $gtf2bed < $output2 > $output4 ;
//         //$bedops -p $output3 $output4 | cut -f1-3 | uniq - | awk '{if(\$3 - \$2 > 2){print}}' - > $output5 ;
//         //python $code_base/get_intersect_bed.py $output3 $output4 $output5 > $output7 ;
//
//build_STAR_reference = {
//    output.dir=branch.name+"/STARRef"
//    produce("Genome"){
//        exec """
//            module load star ;
//            cd $branch.name ;
//            if [ ! -d STARRef ]; then mkdir STARRef ; fi ;
//            STAR --runMode genomeGenerate --genomeDir STARRef --sjdbFileChrStartEnd ../$input.juncs
//                --sjdbOverhang 99 --genomeFastaFiles ../$input.fasta  --runThreadN $threads  --genomeSAindexNbases 5 ;
//        """
//    }
//}

star_align = {
   def out_prefix=branch.name+"/clinker_out/alignment/"
   def workingDir=System.getProperty("user.dir");
   def read_files=inputs.fastq.gz.split().collect { workingDir+"/$it" }.join(' ')
   output.dir=branch.name+"/clinker_out/alignment"
   if(type=="controls"){
        output.dir=branch.parent.parent.name+"/clinker_out/alignment/controls/"
        sample_name=read_files.split()[0].split('/').last().split('\\.').first()
        out_prefix=output.dir+'/'+sample_name+'_'
   }
   produce(out_prefix+"Aligned.sortedByCoord.out.bam",out_prefix+"SJ.out.tab"){
    exec """
        mkdir -p $output.dir ;
        module load star ;
        module load samtools ;
        time STAR --genomeDir $genome_folder
           --readFilesCommand zcat
           --readFilesIn $read_files
           --outSAMtype BAM SortedByCoordinate
           --outFileNamePrefix $out_prefix
           --runThreadN $threads
           --limitBAMsortRAM 5050000000
           --genomeSAindexNbases 5
           --outWigStrand Unstranded
           --outWigType bedGraph
           --outWigNorm RPM ;
        samtools index ${workingDir}/${output1} ;
        rm -rf ${out_prefix}_STARtmp ;
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
              align_contigs_against_transcriptome +
              filter_contigs_against_transcriptome +
              create_salmon_index +
              [run_salmon, "controls/%.*.fastq.gz" * [ run_salmon.using(type:"controls") ]] +
              create_ec_count_matrix +
              run_diffsplice +
              filter_on_significant_ecs +
              annotate_diffspliced_contigs +
              create_supertranscript_reference +
              annotate_supertranscript +
              make_supertranscript_gmap_reference +
              align_contigs_to_supertranscript +
              star_genome_gen + [star_align, "controls/%.*.fastq.gz" * [ star_align.using(type:"controls") ]]]
//              run_lace +
//              annotate_superTranscript +
//              build_STAR_reference +
//              map_reads
}
