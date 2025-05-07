#!/usr/bin/env bash


###############################################################
### For 16S
###############################################################


WORKING_DIRECTORY=/scratch_vol0/fungi/valormicro/05_QIIME2
OUTPUT=/scratch_vol0/fungi/valormicro/05_QIIME2/visual

DATABASE=/scratch_vol0/fungi/valormicro/98_database_files
TMPDIR=/scratch_vol0


# Aim: classify reads by taxon using a fitted classifier

# https://docs.qiime2.org/2019.10/tutorials/moving-pictures/
# In this step, you will take the denoised sequences from step 5 (rep-seqs.qza) and assign taxonomy to each sequence (phylum -> class -> …genus -> ). 
# This step requires a trained classifer. You have the choice of either training your own classifier using the q2-feature-classifier or downloading a pretrained classifier.

# https://docs.qiime2.org/2019.10/tutorials/feature-classifier/


# Aim: Import data to create a new QIIME 2 Artifact
# https://gitlab.com/IAC_SolVeg/CNRT_BIOINDIC/-/blob/master/snk/12_qiime2_taxonomy


###############################################################
### For Bacteria
###############################################################

cd $WORKING_DIRECTORY

eval "$(conda shell.bash hook)"
conda activate qiime2-2021.4

# I'm doing this step in order to deal the no space left in cluster :
export TMPDIR='/scratch_vol0/fungi'
echo $TMPDIR

# Make the directory (mkdir) only if not existe already(-p)
mkdir -p taxonomy/16S
mkdir -p export/taxonomy/16S

# I'm doing this step in order to deal the no space left in cluster :
export TMPDIR='/scratch_vol0/fungi'
echo $TMPDIR

###### All this step was for "old" database, now we uysed new ones 
######
######
######
######qiime tools import --type 'FeatureData[Taxonomy]' \
######  --input-format HeaderlessTSVTaxonomyFormat \
######  --input-path /Users/pierre-louisstenger/Documents/PostDoc_02_MetaBarcoding_IAC/02_Data/07_valormicro/valormicro/98_database_files/silva_nr99_v138_wSpecies_train_set.fa \
######  --output-path taxonomy/RefTaxo.qza
######
######qiime tools import --type 'FeatureData[Sequence]' \
######  --input-path /Users/pierre-louisstenger/Documents/PostDoc_02_MetaBarcoding_IAC/02_Data/07_valormicro/valormicro/98_database_files/silva_nr99_v138_wSpecies_train_set.fa \
######  --output-path taxonomy/DataSeq.qza
######
######   
####### Aim: Extract sequencing-like reads from a reference database.
####### Warning: For v4 only !!! Not for its2 !!! 
######
####### The --p-trunc-len parameter should only be used to trim reference sequences,
####### if query sequences are trimmed to this same length or shorter.
###### 
####### Paired sequences that successfully join will typically be variable in length.
####### Single reads not truncated at a specific length may also be variable in length.
######
####### For classification of paired-end reads and untrimmed single-end reads,
####### we recommend training a classifier on sequences that have been extracted
####### at the appropriate primer sites, but are not trimmed !!!
####### -----
####### The primer sequences used for extracting reads should be the actual DNA-binding
####### (i.e., biological) sequence contained within a primer construct.
######
####### It should NOT contain any non-biological, non-binding sequence,
####### e.g., adapter, linker, or barcode sequences.
######
####### If you aren't sure what section of your primer sequences are actual DNA-binding
####### you should consult whoever constructed your sequencing library, your sequencing
####### center, or the original source literature on these primers.
######
####### If your primer sequences are > 30 nt long, they most likely contain some
####### non-biological sequence !
######
######qiime feature-classifier extract-reads --i-sequences taxonomy/DataSeq.qza \
######        --p-f-primer 'GTGCCAGCMGCCGCGGTAA' \
######        --p-r-primer 'TCCTCCGCTTATTGATATGC' \
######        --o-reads taxonomy/RefSeq.qza 
######        
######        #--p-trunc-len {params.length} \
######
####### Aim: Create a scikit-learn naive_bayes classifier for reads
######
######qiime feature-classifier fit-classifier-naive-bayes \
######  --i-reference-reads taxonomy/RefSeq.qza \
######  --i-reference-taxonomy taxonomy/RefTaxo.qza \
######  --o-classifier taxonomy/Classifier.qza

# With new database :

# See here for only V4 : https://www.dropbox.com/sh/nz7c5asn6b3hr1j/AADMAR-YZOBkpUQJLumZ9w3wa/ver_0.02?dl=0&subfolder_nav_tracking=1
# See here for all 16S : https://www.dropbox.com/sh/ibpy9j0clw8dzwm/AAAIVuYnqUzAOxlg2fijePQna/ver_0.02?dl=0&subfolder_nav_tracking=1

# See this thread https://forum.qiime2.org/t/silva-138-classifiers/13131 (found because of this thread : https://forum.qiime2.org/t/silva-138-for-qiime2/12957/4)

#cp $DATABASE/SILVA-138-SSURef-full-length-classifier.qza taxonomy/Classifier.qza
cp /home/fungi/Mayotte_microorganism_colonisation/98_database_files/SILVA-138-SSURef-Full-Seqs.qza taxonomy/16S/DataSeq.qza
cp /home/fungi/Mayotte_microorganism_colonisation/98_database_files/Silva-v138-full-length-seq-taxonomy.qza taxonomy/16S/RefTaxo.qza

# Script Nolwenn
#R1_Primers = c("GTGCCAGCMGCCGCGGTAA","GTGYCAGCMGCCGCGGTAA")
#R2_Primers = c("GGACTACHVGGGTWTCTAAT","GGACTACNVGGGTWTCTAAT")

# Here only for V4 --> forward: 'GTGCCAGCMGCCGCGGTAA'  # 515f & reverse: 'GGACTACHVGGGTWTCTAAT' # 806r
#qiime feature-classifier extract-reads --i-sequences taxonomy/DataSeq.qza \
#        --p-f-primer 'GTGCCAGCMGCCGCGGTAA' \
#        --p-r-primer 'TCCTCCGCTTATTGATATGC' \
#        --o-reads taxonomy/RefSeq.qza 

# Here for V1V2V3V4 --> 27F 'AGAGTTTGATCCTGGCTCAG' & reverse: 'GGACTACHVGGGTWTCTAAT' # 806r
#qiime feature-classifier extract-reads --i-sequences taxonomy/DataSeq.qza \
#        --p-f-primer 'AGAGTTTGATCCTGGCTCAG' \
#        --p-r-primer 'TCCTCCGCTTATTGATATGC' \
#        --o-reads taxonomy/RefSeq.qza         


# If necessary :
# https://forum.qiime2.org/t/available-pre-trained-classifier-of-v3-v4-341f-805r-region-with-gg-99/3275
# Available: Pre-trained classifier of V3-V4 (341F, 805R) region with gg_99

# 16S : V3/V4 : V3V4 (amplified with primers 341F–805R)
###qiime feature-classifier extract-reads --i-sequences taxonomy/DataSeq.qza \
###        --p-f-primer 'CCTACGGGNGGCWGCAG' \
###        --p-r-primer 'GACTACHVGGGTATCTAATCC' \
###        --o-reads taxonomy/RefSeq.qza 
###

# According ADNiD: Caporaso et al. (1), 515f Original and 806r Original


# From Clarisse and from bioinformatics <bioinformatics@microsynth.ch>
# The sequences for the V34 primers were:
# 341F    CCTACGGGNGGCWGCAG       805R    GACTACHVGGGTATCTAATCC
qiime feature-classifier extract-reads --i-sequences taxonomy/16S/DataSeq.qza \
        --p-f-primer 'CCTACGGGNGGCWGCAG' \
        --p-r-primer 'GACTACHVGGGTATCTAATCC' \
        --o-reads taxonomy/16S/RefSeq.qza 


# Aim: Create a scikit-learn naive_bayes classifier for reads

qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads taxonomy/16S/RefSeq.qza \
  --i-reference-taxonomy taxonomy/16S/RefTaxo.qza \
  --o-classifier taxonomy/16S/Classifier.qza
  
# Aim: Create a scikit-learn naive_bayes classifier for reads

qiime feature-classifier classify-sklearn \
   --i-classifier taxonomy/16S/Classifier.qza \
   --i-reads core/ConRepSeq.qza \
   --o-classification taxonomy/16S/taxonomy_reads-per-batch_ConRepSeq.qza
   
qiime feature-classifier classify-sklearn \
  --i-classifier taxonomy/16S/Classifier.qza \
  --i-reads core/RepSeq.qza \
  --o-classification taxonomy/16S/taxonomy_reads-per-batch_RepSeq.qza

qiime feature-classifier classify-sklearn \
  --i-classifier taxonomy/16S/Classifier.qza \
  --i-reads core/RarRepSeq.qza \
  --o-classification taxonomy/16S/taxonomy_reads-per-batch_RarRepSeq.qza

# Switch to https://chmi-sops.github.io/mydoc_qiime2.html#step-9-assign-taxonomy
# --p-reads-per-batch 0 (default)

qiime metadata tabulate \
  --m-input-file taxonomy/16S/taxonomy_reads-per-batch_RarRepSeq.qza \
  --o-visualization taxonomy/16S/taxonomy_reads-per-batch_RarRepSeq.qzv

qiime metadata tabulate \
  --m-input-file taxonomy/16S/taxonomy_reads-per-batch_ConRepSeq.qza \
  --o-visualization taxonomy/16S/taxonomy_reads-per-batch_ConRepSeq.qzv
  
qiime metadata tabulate \
  --m-input-file taxonomy/16S/taxonomy_reads-per-batch_RepSeq.qza \
  --o-visualization taxonomy/16S/taxonomy_reads-per-batch_RepSeq.qzv  

# Now create a visualization of the classified sequences.
  
qiime taxa barplot \
  --i-table core/Table.qza \
  --i-taxonomy taxonomy/16S/taxonomy_reads-per-batch_RepSeq.qza \
  --m-metadata-file $DATABASE/sample-metadata.tsv \
  --o-visualization taxonomy/16S/16Staxa-bar-plots_reads-per-batch_RepSeq.qzv

qiime taxa barplot \
  --i-table core/ConTable.qza \
  --i-taxonomy taxonomy/16S/taxonomy_reads-per-batch_ConRepSeq.qza \
  --m-metadata-file $DATABASE/sample-metadata.tsv \
  --o-visualization taxonomy/16S/16Staxa-bar-plots_reads-per-batch_ConRepSeq.qzv
  
qiime taxa barplot \
  --i-table core/RarTable.qza \
  --i-taxonomy taxonomy/16S/taxonomy_reads-per-batch_RarRepSeq.qza \
  --m-metadata-file $DATABASE/sample-metadata.tsv \
  --o-visualization taxonomy/16S/16Staxa-bar-plots_reads-per-batch_RarRepSeq.qzv  

qiime tools export --input-path taxonomy/16S/Classifier.qza --output-path export/taxonomy/16S/Classifier
qiime tools export --input-path taxonomy/16S/RefSeq.qza --output-path export/taxonomy/16S/RefSeq
qiime tools export --input-path taxonomy/16S/16SDataSeq.qza --output-path export/taxonomy/16S/16SDataSeq
qiime tools export --input-path taxonomy/16S/RefTaxo.qza --output-path export/taxonomy/16S/RefTaxo
  
qiime tools export --input-path taxonomy/16S/16Staxa-bar-plots_reads-per-batch_RarRepSeq.qzv --output-path export/taxonomy/16S/16Staxa-bar-plots_reads-per-batch_RarRepSeq
qiime tools export --input-path taxonomy/16S/16Staxa-bar-plots_reads-per-batch_ConRepSeq.qzv --output-path export/taxonomy/16S/16Staxa-bar-plots_reads-per-batch_ConRepSeq
qiime tools export --input-path taxonomy/16S/16Staxa-bar-plots_reads-per-batch_RepSeq.qzv --output-path export/taxonomy/16S/16Staxa-bar-plots_reads-per-batch_RepSeq

qiime tools export --input-path taxonomy/16S/taxonomy_reads-per-batch_RepSeq.qzv --output-path export/taxonomy/16S/taxonomy_reads-per-batch_RepSeq_visual
qiime tools export --input-path taxonomy/16S/taxonomy_reads-per-batch_ConRepSeq.qzv --output-path export/taxonomy/16S/taxonomy_reads-per-batch_ConRepSeq_visual
qiime tools export --input-path taxonomy/16S/taxonomy_reads-per-batch_RarRepSeq.qzv --output-path export/taxonomy/16S/taxonomy_reads-per-batch_RarRepSeq_visual

qiime tools export --input-path taxonomy/16S/taxonomy_reads-per-batch_RepSeq.qza --output-path export/taxonomy/16S/taxonomy_reads-per-batch_RepSeq
qiime tools export --input-path taxonomy/16S/taxonomy_reads-per-batch_ConRepSeq.qza --output-path export/taxonomy/16S/taxonomy_reads-per-batch_ConRepSeq
qiime tools export --input-path taxonomy/16S/taxonomy_reads-per-batch_RarRepSeq.qza --output-path export/taxonomy/16S/taxonomy_reads-per-batch_RarRepSeq

