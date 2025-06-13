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

# V4V5: 16S rDNA 515F-Y (5′-GTGYCAGCMGCCGCGGTAA-3′) and 926R (5′-CCGYCAATTYMTTTRAGTTT-3′).

qiime feature-classifier extract-reads \
    --i-sequences /scratch_vol0/fungi/dugong_microbiome/05_QIIME2/silva-138.2-ssu-nr99-seqs-derep-uniq.qza \
    --p-f-primer GTGYCAGCMGCCGCGGTAA \
    --p-r-primer CCGYCAATTYMTTTRAGTTT \
    --p-n-jobs 2 \
    --p-read-orientation 'forward' \
    --o-reads silva-138.2-ssu-nr99-seqs-515f-926r.qza

qiime rescript dereplicate \
    --i-sequences silva-138.2-ssu-nr99-seqs-515f-926r.qza \
    --i-taxa silva-138.2-ssu-nr99-tax-derep-uniq.qza \
    --p-mode 'uniq' \
    --o-dereplicated-sequences silva-138.2-ssu-nr99-seqs-515f-926r-uniq.qza \
    --o-dereplicated-taxa  silva-138.2-ssu-nr99-tax-515f-926r-derep-uniq.qza

# Aim: Create a scikit-learn naive_bayes classifier for reads
#qiime feature-classifier fit-classifier-naive-bayes \
#  --i-reference-reads taxonomy/16S/RefSeq.qza \
#  --i-reference-taxonomy taxonomy/16S/RefTaxo.qza \
#  --o-classifier taxonomy/16S/Classifier.qza
  
# Aim: Create a scikit-learn naive_bayes classifier for reads
qiime feature-classifier fit-classifier-naive-bayes \
    --i-reference-reads silva-138.2-ssu-nr99-seqs-515f-926r-uniq.qza \
    --i-reference-taxonomy silva-138.2-ssu-nr99-tax-515f-926r-derep-uniq.qza \
    --o-classifier silva-138.2-ssu-nr99-515f-926r-classifier.qza
    
scp -r silva-138.2-ssu-nr99-515f-926r-classifier.qza taxonomy/16S/Classifier.qza
    
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
