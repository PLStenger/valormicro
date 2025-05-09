#!/usr/bin/env bash


###############################################################
### For 16S
###############################################################


WORKING_DIRECTORY=/scratch_vol0/fungi/valormicro/05_QIIME2
OUTPUT=/scratch_vol0/fungi/valormicro/05_QIIME2/visual

DATABASE=/scratch_vol0/fungi/valormicro/98_database_files
TMPDIR=/scratch_vol0

#QIIME="singularity exec --cleanenv -B /scratch_vol0:/scratch_vol0 /scratch_vol0/fungi/qiime2_images/qiime2-2024.5.sif qiime"


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
conda activate /scratch_vol0/fungi/envs/qiime2-amplicon-2024.10
#conda activate qiime2-2021.4

#export PYTHONPATH="${PYTHONPATH}:/scratch_vol0/fungi/.local/lib/python3.9/site-packages/"
#echo $PYTHONPATH

# I'm doing this step in order to deal the no space left in cluster :
export TMPDIR='/scratch_vol0/fungi'
echo $TMPDIR

# Make the directory (mkdir) only if not existe already(-p)
mkdir -p taxonomy/16S
mkdir -p export/taxonomy/16S

#singularity exec --cleanenv --env TMPDIR=/scratch_vol0/fungi

####### OK ! Je # pour gagner du temps ici
###### Code from: https://forum.qiime2.org/t/processing-filtering-and-evaluating-the-silva-database-and-other-reference-sequence-data-with-rescript/15494
######qiime rescript get-silva-data \
######    --p-version '138.2' \
######    --p-target 'SSURef_NR99' \
######    --o-silva-sequences silva-138.2-ssu-nr99-rna-seqs.qza \
######    --o-silva-taxonomy silva-138.2-ssu-nr99-tax.qza

####### If you'd like to be able to jump to steps that only take FeatureData[Sequence] as input you can convert your data to FeatureData[Sequence] like so:
######qiime rescript reverse-transcribe \
######    --i-rna-sequences silva-138.2-ssu-nr99-rna-seqs.qza \
######    --o-dna-sequences silva-138.2-ssu-nr99-seqs.qza
######
####### “Culling” low-quality sequences with cull-seqs
####### Here we’ll remove sequences that contain 5 or more ambiguous bases (IUPAC compliant ambiguity bases) and any homopolymers that are 8 or more bases in length. These are the default parameters.
######qiime rescript cull-seqs \
######    --i-sequences silva-138.2-ssu-nr99-seqs.qza \
######    --o-clean-sequences silva-138.2-ssu-nr99-seqs-cleaned.qza
######
####### Filtering sequences by length and taxonomy
######qiime rescript filter-seqs-length-by-taxon \
######    --i-sequences silva-138.2-ssu-nr99-seqs-cleaned.qza \
######    --i-taxonomy silva-138.2-ssu-nr99-tax.qza \
######    --p-labels Archaea Bacteria Eukaryota \
######    --p-min-lens 900 1200 1400 \
######    --o-filtered-seqs silva-138.2-ssu-nr99-seqs-filt.qza \
######    --o-discarded-seqs silva-138.2-ssu-nr99-seqs-discard.qza 
######
#######Dereplicating in uniq mode
######qiime rescript dereplicate \
######    --i-sequences silva-138.2-ssu-nr99-seqs-filt.qza  \
######    --i-taxa silva-138.2-ssu-nr99-tax.qza \
######    --p-mode 'uniq' \
######    --o-dereplicated-sequences silva-138.2-ssu-nr99-seqs-derep-uniq.qza \
######    --o-dereplicated-taxa silva-138.2-ssu-nr99-tax-derep-uniq.qza
######
######
######qiime feature-classifier fit-classifier-naive-bayes \
######  --i-reference-reads  silva-138.2-ssu-nr99-seqs-derep-uniq.qza \
######  --i-reference-taxonomy silva-138.2-ssu-nr99-tax-derep-uniq.qza \
######  --o-classifier silva-138.2-ssu-nr99-classifier.qza    
######
######
####### Here only for V4 --> forward: 'GTGCCAGCMGCCGCGGTAA'  # 515f & reverse: 'GGACTACHVGGGTWTCTAAT' # 806r
#######qiime feature-classifier extract-reads --i-sequences taxonomy/DataSeq.qza \
#######        --p-f-primer 'GTGCCAGCMGCCGCGGTAA' \
#######        --p-r-primer 'TCCTCCGCTTATTGATATGC' \
#######        --o-reads taxonomy/RefSeq.qza 
######
####### Here for V1V2V3V4 --> 27F 'AGAGTTTGATCCTGGCTCAG' & reverse: 'GGACTACHVGGGTWTCTAAT' # 806r
#######qiime feature-classifier extract-reads --i-sequences taxonomy/DataSeq.qza \
#######        --p-f-primer 'AGAGTTTGATCCTGGCTCAG' \
#######        --p-r-primer 'TCCTCCGCTTATTGATATGC' \
#######        --o-reads taxonomy/RefSeq.qza         
######
######
####### If necessary :
####### https://forum.qiime2.org/t/available-pre-trained-classifier-of-v3-v4-341f-805r-region-with-gg-99/3275
####### Available: Pre-trained classifier of V3-V4 (341F, 805R) region with gg_99
######
####### 16S : V3/V4 : V3V4 (amplified with primers 341F–805R)
#########qiime feature-classifier extract-reads --i-sequences taxonomy/DataSeq.qza \
#########        --p-f-primer 'CCTACGGGNGGCWGCAG' \
#########        --p-r-primer 'GACTACHVGGGTATCTAATCC' \
#########        --o-reads taxonomy/RefSeq.qza 
#########
######
####### According ADNiD: Caporaso et al. (1), 515f Original and 806r Original
######
####### From Clarisse and from bioinformatics <bioinformatics@microsynth.ch>
####### The sequences for the V34 primers were:
####### 341F    CCTACGGGNGGCWGCAG       805R    GACTACHVGGGTATCTAATCC
#######qiime feature-classifier extract-reads --i-sequences taxonomy/16S/DataSeq.qza \
#######        --p-f-primer 'CCTACGGGNGGCWGCAG' \
#######        --p-r-primer 'GACTACHVGGGTATCTAATCC' \
#######        --o-reads taxonomy/16S/RefSeq.qza 
######
######
####### Original from tuto:
#######    --p-f-primer GTGYCAGCMGCCGCGGTAA \
#######    --p-r-primer GGACTACNVGGGTWTCTAAT \
######
####### From Clarisse and from bioinformatics <bioinformatics@microsynth.ch>
####### The sequences for the V34 primers were:
####### 341F    CCTACGGGNGGCWGCAG       805R    GACTACHVGGGTATCTAATCC
######
######qiime feature-classifier extract-reads \
######    --i-sequences silva-138.2-ssu-nr99-seqs-derep-uniq.qza \
######    --p-f-primer CCTACGGGNGGCWGCAG \
######    --p-r-primer GACTACHVGGGTATCTAATCC \
######    --p-n-jobs 2 \
######    --p-read-orientation 'forward' \
######    --o-reads silva-138.2-ssu-nr99-seqs-341f-805r.qza
######
######qiime rescript dereplicate \
######    --i-sequences silva-138.2-ssu-nr99-seqs-341f-805r.qza \
######    --i-taxa silva-138.2-ssu-nr99-tax-derep-uniq.qza \
######    --p-mode 'uniq' \
######    --o-dereplicated-sequences silva-138.2-ssu-nr99-seqs-341f-805r-uniq.qza \
######    --o-dereplicated-taxa  silva-138.2-ssu-nr99-tax-341f-805r-derep-uniq.qza
######
####### Aim: Create a scikit-learn naive_bayes classifier for reads
#######qiime feature-classifier fit-classifier-naive-bayes \
#######  --i-reference-reads taxonomy/16S/RefSeq.qza \
#######  --i-reference-taxonomy taxonomy/16S/RefTaxo.qza \
#######  --o-classifier taxonomy/16S/Classifier.qza
######  
####### Aim: Create a scikit-learn naive_bayes classifier for reads
######qiime feature-classifier fit-classifier-naive-bayes \
######    --i-reference-reads silva-138.2-ssu-nr99-seqs-341f-805r-uniq.qza \
######    --i-reference-taxonomy silva-138.2-ssu-nr99-tax-341f-805r-derep-uniq.qza \
######    --o-classifier silva-138.2-ssu-nr99-341f-805r-classifier.qza
    
scp -r /scratch_vol0/fungi/dugong_microbiome/05_QIIME2/silva-138.2-ssu-nr99-341f-805r-classifier.qza taxonomy/16S/Classifier.qza

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
