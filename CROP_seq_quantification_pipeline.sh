#!/bin/bash

# Load required modules
module load cellranger/7.0.0
module load bcl2fastq2/v2.20.0.422

#1.Make fastq files of each lane
cellranger mkfastq --run=/home/zjiayu/CROP-seq/new_data/ --id=CROP \
--samplesheet=/home/zjiayu/CROP-seq/new_data/samplesheet.csv \
--output-dir=/home/zjiayu/CROP-seq/new_data/Count_output/

#2. Make reference genome by adding sgRNAs to human genome.
U6P="GAGGGCCTATTTCCCATGATTCCTTCATATTTGCATATACGATACAAGGCTGTTAGAGAGATAATTAGAATTAATTTGACTGTAAACACAAAGATATTAGTACAAAATACGTGACGTAGAAAGTAATAATTTCTTGGGTAGTTTGCAGTTTTAAAATTATGTTTTAAAATGGACTATCATATGCTTACCGTAACTTGAAAGTATTTCGATTTCTTGGCTTTATATATCTTGTGGAAAGGAC"
#The U6 promoter need to add GAAACACCG at end to link the sgRNA, finally is 250bp.
U6Pplus="GAGGGCCTATTTCCCATGATTCCTTCATATTTGCATATACGATACAAGGCTGTTAGAGAGATAATTAGAATTAATTTGACTGTAAACACAAAGATATTAGTACAAAATACGTGACGTAGAAAGTAATAATTTCTTGGGTAGTTTGCAGTTTTAAAATTATGTTTTAAAATGGACTATCATATGCTTACCGTAACTTGAAAGTATTTCGATTTCTTGGCTTTATATATCTTGTGGAAAGGACGAAACACCG"
#The subset LTR region is short to 230bp.
LTR="GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGCTTTTTTAAGCTTGGCGTAACTAGATCTTGAGACACTGCTTTTTGCTTGTACTGGGTCTCTCTGGTTAGACCAGATCTGAGCCTGGGAGCTCTCTGGCTAACTAGGGAACCCACTGCTTAAGCCTCAATAAAGCTTGCCTTGAGTGCTTCAAGTA"
#
cat /home/zjiayu/CROP-seq/new_data/hg38_GFP/sgRNA.txt | awk -vOFS='\t' -vLTR=$LTR -vU6=$U6Pplus '{print ">"$1;print U6$2LTR}' > /home/zjiayu/CROP-seq/new_data_another_method/hg38_GFP/CROP_sgRNA_extend_500bp.fa
#cut them into 60bp in each line
cat /home/zjiayu/CROP-seq/new_data/hg38_GFP/CROP_sgRNA_extend_500bp.fa | awk -F ' ' '(NR%2==1){print $0}(NR%2==0){for(i=1;i<=10;i++){start=(i-1)*60+1;len=60;if(start<length($0)){if(end>length($0)){print substr($0,start,length($0)-start)}else{print substr($0,start,len)}}}}' > /home/zjiayu/CROP-seq/new_data_another_method/hg38_GFP/CROP_sgRNA_extend_500bp_format.fa

#Add gRNA sequence to the reference genome(.fasta file)
# Define paths to input and output files
GENOME_FA="/home/zjiayu/CROP-seq/new_data/hg38_GFP/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
SGRNA_FA="/home/zjiayu/CROP-seq/new_data/hg38_GFP/CROP_sgRNA_extend_500bp_format.fa"
OUTPUT_FA="/home/zjiayu/CROP-seq/new_data/hg38_GFP/Homo_sapiens.GRCh38.mix_sgRNA.fa.gz"

# Process the genome file to replace metadata, add "chr" prefix, and handle mitochondrial chromosome
zcat $GENOME_FA \
    | sed -E 's/^>(\S+).*/>\1 \1/' \
    > tmp.fa
    
# Append the sgRNA sequence to the processed genome file and gzip the result
cat tmp.fa $SGRNA_FA | gzip -c > ${OUTPUT_FA}

# Clean up temporary files
rm tmp.fa

gunzip ${OUTPUT_FA}

#add sgRNA to the human reference GTF.
#first process the GTF file as metioned in 10x cell ranger notes.
ID="(ENS(MUS)?[GTE][0-9]+)\.([0-9]+)"
# Construct the gene ID allowlist. We filter the list of all transcripts
# based on these criteria:
#   - allowable gene_type (biotype)
#   - allowable transcript_type (biotype)
#   - no "PAR" tag (only present for Y chromosome PAR)
#   - no "readthrough_transcript" tag
BIOTYPE_PATTERN="(protein_coding|lncRNA|\
  IG_C_gene|IG_D_gene|IG_J_gene|IG_LV_gene|IG_V_gene|\
  IG_V_pseudogene|IG_J_pseudogene|IG_C_pseudogene|\
  TR_C_gene|TR_D_gene|TR_J_gene|TR_V_gene|\
  TR_V_pseudogene|TR_J_pseudogene)"
  GENE_PATTERN="gene_biotype \"${BIOTYPE_PATTERN}\""
  TX_PATTERN="transcript_biotype \"${BIOTYPE_PATTERN}\""
  READTHROUGH_PATTERN="tag \"readthrough_transcript\""
  PAR_PATTERN="tag \"PAR\""

#Transcript passing the filters (filter necessary genes, and this step is optional)
zcat /home/zjiayu/CROP-seq/new_data/hg38_GFP/Homo_sapiens.GRCh38.105.gtf.gz \
    | awk '$3 == "transcript"' \
    | grep -E "$GENE_PATTERN" \
    | grep -E "$TX_PATTERN" \
    | grep -Ev "$READTHROUGH_PATTERN" \
    | grep -Ev "$PAR_PATTERN" \
    | sed -E 's/.*(gene_id "[^"]+").*/\1/' \
    | sort \
    | uniq \
    > "/home/zjiayu/CROP-seq/new_data/hg38_GFP/gene_allowlist"


#Filter to the gene allowlist
zcat /home/zjiayu/CROP-seq/new_data/hg38_GFP/Homo_sapiens.GRCh38.105.gtf.gz | \
grep -Ff "/home/zjiayu/CROP-seq/new_data/hg38_GFP/gene_allowlist" - >> /home/zjiayu/CROP-seq/new_data/hg38_GFP/Homo_sapiens.GRCh38.105.gtf.filtered

#merge with other artificial sequence gtf (based on the format of gtf file)
cat /home/zjiayu/CROP-seq/new_data/hg38_GFP/sgRNA.txt | awk -vOFS='\t' '{print $1,"ensembl_havana","gene","1","500",".","+",".","gene_id \""$1"\"; gene_version \"11\"; gene_name \""$1"\"; gene_source \"ensembl_havana\"; gene_biotype \"protein_coding\";";print $1,"ensembl_havana","transcript","1","500",".","+",".","gene_id \""$1"\"; gene_version \"11\"; transcript_id \""$1"\"; transcript_version \"4\" ;gene_name \""$1"\"; gene_source \"ensembl_havana\"; gene_biotype \"protein_coding\"; transcript_name \""$1"\"; transcript_source \"ensembl_havana\";";print $1,"ensembl_havana","exon","1","500",".","+",".","gene_id \""$1"\"; gene_version \"11\"; transcript_id \""$1"\"; transcript_version \"4\" ; exon_number "1"; gene_name \""$1"\"; gene_source \"ensembl_havana\"; gene_biotype \"protein_coding\"; transcript_name \""$1"\"; transcript_source \"ensembl_havana\";";print $1,"ensembl_havana","CDS","1","500",".","+",".","gene_id \""$1"\"; gene_version \"11\"; transcript_id \""$1"\"; transcript_version \"4\" ; exon_number "1"; gene_name \""$1"\"; gene_source \"ensembl_havana\"; gene_biotype \"protein_coding\"; transcript_name \""$1"\"; transcript_source \"ensembl_havana\";"}' > /home/zjiayu/CROP-seq/new_data_another_method/hg38_GFP/CROP_sgRNA_extend_500bp.gtf

cat /home/zjiayu/CROP-seq/new_data/hg38_GFP/Homo_sapiens.GRCh38.105.gtf.filtered /home/zjiayu/CROP-seq/new_data/hg38_GFP/CROP_sgRNA_extend_500bp.gtf | gzip -c > /home/zjiayu/CROP-seq/new_data/hg38_GFP/Homo_sapiens.GRCh38.mix_sgRNA.gtf.gz
gunzip /home/zjiayu/CROP-seq/new_data/hg38_GFP/Homo_sapiens.GRCh38.mix_sgRNA.gtf.gz

module load cellranger/7.0.0
# Create reference package
cellranger mkref \
--genome=CROP_genome \
--fasta=/home/zjiayu/CROP-seq/new_data/hg38_GFP/Homo_sapiens.GRCh38.mix_sgRNA.fa \
--genes=/home/zjiayu/CROP-seq/new_data/hg38_GFP/Homo_sapiens.GRCh38.mix_sgRNA.gtf

#map to reference with sgRNA
cellranger count --localcores=10 --id=ENRICH \
--sample=library_2 --fastqs=/home/zjiayu/CROP-seq/new_data/Fastq_output/HM7L2AFX3// \
--transcriptome=/home/zjiayu/CROP-seq/new_data/hg38_GFP/CROP_genome/ \
--output-dir=/home/zjiayu/CROP-seq/new_data/Count_output/

#process the scRNA-seq part
cellranger count --localcores=10 --id=CROP_scRNA \
--sample=library_1 --fastqs=/home/zjiayu/CROP-seq/new_data/Fastq_output/HM7L2AFX3/Library_1/ \
--transcriptome=/home/zjiayu/CROP-seq/new_data/hg38_GFP/CROP_genome/ \
--output-dir=/home/zjiayu/CROP-seq/new_data/Count_output/

#Assign to cells
module load python/3.7

#this is the python file which available from the https://github.com/shendurelab/single-cell-ko-screens/blob/master/get_barcodes.py, and whitelist is the single guide RNA sequence for each target gene. 
python /home/zjiayu/CROP-seq/new_data/Count_output/get_barcodes.py --all_reads --input_bams /home/zjiayu/CROP-seq/new_data/Count_output/ENRICH/outs/possorted_genome_bam.bam -o ko_barcodes.txt --whitelist /home/zjiayu/CROP-seq/new_data_another_method/Count_output/whitelist.txt --search_seq CTTGTGGAAAGGACGAAACACCG
