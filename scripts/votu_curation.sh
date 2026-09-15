#########################################
#########################################
########## GeNomad and CheckV ###########
#########################################
#########################################


#geNomad, to predict viruses
sbatch -p base -t 48:00:00 --mem=125000 -J genomad -n 15 --wrap="genomad end-to-end --threads 15 input_fasta genomad_outdir path/to/genomad_db"
#CheckV, to assess quality of predicted viruses
sbatch -p base -t 48:00:00 --mem=125000 -J checkv -n 15 --wrap="checkv end_to_end genomad_outdir/virus.fna checkv_outdir -t 10 -d path/to/checkv_db"
#CheckV AGAIN, to assess quality of trimmed proviruses
sbatch -p base -t 48:00:00 --mem=125000 -J checkv -n 15 --wrap="checkv end_to_end checkv_outdir/proviruses.fna checkv_outdir/proviruses_checkv_outdir -t 10 -d path/to/checkv_db"

# Get only good quality viruses, no proviruses
awk '$3 == "No" && ($0 ~ /Complete|Medium|High-quality/)' checkv_outdir/quality_summary.tsv | cut -f1 | grep -v "contig_id" > comp.med.hq.viruses.txt

# Get good quality proviruses
awk '$3 == "No" && ($0 ~ /Complete|Medium|High-quality/)' checkv_outdir/proviruses_checkv_outdir/quality_summary.tsv | cut -f1 | grep -v "contig_id" > comp.med.hq.proviruses.txt

cat comp.med.hq.viruses.txt comp.med.hq.proviruses.txt > comp.med.hq.all.viruses.txt

cat checkv_outdir/viruses.fna checkv_outdir/proviruses_checkv_outdir/proviruses.fna > all.viruses.fa

seqtk subseq all.viruses.fa comp.med.hq.all.viruses.txt > comp.med.hq.all.viruses.fa

#########################################
#########################################
##########     Clustering     ###########
#########################################
#########################################

# Step 1: Create a blast+ database for each sets of datasets to cluster
  makeblastdb -in comp.med.hq.all.viruses.fa -dbtype nucl -out blast_output.db

# Step 2: Use megablast from blast+ package to perform all-vs-all blastn of sequences
  sbatch -p base -t 12:00:00 --mem=250000 -J blastn -n 30 --wrap="blastn -query comp.med.hq.all.viruses.fa \
  -db blast_output.db \
  -outfmt '6 std qlen slen' -max_target_seqs 10000 \
  -out blast_output.tsv \
  -num_threads 30"

# Step 3: Calculate pairwise ANI by combining local alignments between sequence pairs
  sbatch -p highmem -t 12:00:00 --mem=250000 -J anicalc -n 30 \
  --wrap="python anicalc.py -i blast_output.tsv \
  -o anicalc.tsv"

# Step 4: Cluster!
python aniclust.py --fna comp.med.hq.all.viruses.fa \
  --ani anicalc.tsv \
  --out aniclust.tsv \
  --min_ani 95 --min_tcov 85 --min_qcov 0

#######################################################
#######################################################
################# Mapping begins here #################
#######################################################
#######################################################

# Get sequences of cluster representatives
cut -f1 anicalc.tsv > cluster_reps.txt
seqtk subseq comp.med.hq.all.viruses.fa cluster_reps.txt > cluster_reps.fa

ref=cluster_reps.fa build=1 -Xmx64g -Xms64g

# For short reads
sbatch -p base -t 48:00:00 --mem=75000 -J bbmap --cpus-per-task=15 --ntasks=1 --wrap="bbmap.sh -Xmx50000m nodisk=t local=f threads=15 in=R1.fastq.gz in2=R2.fastq.gz minid=0.95 \
      out=out.sam mappedonly=t local=f overwrite=t && \
      anvi-init-bam -o out.bam out.sam && \
      rm  out.sam"

# For long reads
sbatch -p base -t 12:00:00 --mem=50000 -J minimap_ant  -n 10 \
  --wrap="minimap2 -ax asm5 --sam-hit-only -t 10 \
  cluster_reps.fa longreads_fastq.gz > minimap2_outdir/longreads.sam"

# Convert .sam to .bam
samtools view -bS minimap2_outdir/longreads.sam > minimap2_outdir/longreads.bam

# Sort .bam files
samtools sort -o minimap2_outdir/longreads.sorted.bam minimap2_outdir/longreads.bam

# Get only primary alignments
samtools view -b -F 0x800 -F 0x100 minimap2_outdir/longreads.sorted.bam > minimap2_outdir/longreads.sorted.primary.bam
done

# Index .bam files
samtools index minimap2_outdir/longreads.sorted.primary.bam

# CoverM
coverm contig -m trimmed_mean --min-covered-fraction 0.25 --bam-files input.bam > output_tm.tsv &
      

#######################################################
#######################################################
################# Get sequence lengths ################
#######################################################
#######################################################

sbatch -p base -t 6:00:00 --mem=50000 -J seqkit -o seqkit.out -e seqkit.err -n 1 \
      --wrap="seqkit stats R1.fasq.gz R2.fastq.gz | sed 's/  */ /g' | sed 's/ /\t/g' | cut -f1,4,5,6,7,8 > output_file.txt"

