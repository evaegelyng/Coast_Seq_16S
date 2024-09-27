from gwf import Workflow
import os, sys
import math
from glob import glob

project_name = "CoastSeq"

gwf = Workflow(defaults={"account": "edna"}) 


gwf.target("qiime1", inputs = ["/home/mpavila/eDNA/faststorage/Velux/CoastSequence/Autumn/Bac16s/base/results/DADA2_nochim.otus"], outputs = ['tmp/DADA_ASVs.qza'],
           cores=4, memory='64g', walltime="24:00:01") << """
source ~/miniconda3/etc/profile.d/conda.sh
conda activate qiime2-2019.10
qiime tools import \
  --type 'FeatureData[Sequence]' \
  --input-path /home/mpavila/eDNA/faststorage/Velux/CoastSequence/Autumn/Bac16s/base/results/DADA2_nochim.otus \
  --output-path /home/mpavila/eDNA/faststorage/Velux/CoastSequence/Autumn/Bac16s/both_silva/tmp/DADA_ASVs.qza
"""

gwf.target("qiime2", inputs = ['tmp/DADA_ASVs.qza','results/classifier_v4.qza'], outputs = ['tmp/taxonomy_long.qza'],
           cores=4, memory='132g', walltime="48:00:01") << """
source ~/miniconda3/etc/profile.d/conda.sh
conda activate qiime2-2019.10
qiime feature-classifier classify-sklearn \
  --i-classifier /home/mpavila/eDNA/faststorage/Velux/CoastSequence/Autumn/Bac16s/both_silva/results/classifier_v4.qza \
  --i-reads /home/mpavila/eDNA/faststorage/Velux/CoastSequence/Autumn/Bac16s/both_silva/tmp/DADA_ASVs.qza \
  --o-classification /home/mpavila/eDNA/faststorage/Velux/CoastSequence/Autumn/Bac16s/both_silva/tmp/taxonomy_long.qza \
  --p-n-jobs -1 \
  --p-confidence 0.85
"""

gwf.target("qiime3", inputs = ['tmp/taxonomy_long.qza'], outputs = ['results/taxonomy_long/taxonomy.tsv'],
           cores=4, memory='64g', walltime="24:00:01") << """
source ~/miniconda3/etc/profile.d/conda.sh
conda activate qiime2-2019.10
qiime tools export \
  --input-path /home/mpavila/eDNA/faststorage/Velux/CoastSequence/Autumn/Bac16s/both_silva/tmp/taxonomy_long.qza \
  --output-path /home/mpavila/eDNA/faststorage/Velux/CoastSequence/Autumn/Bac16s/both_silva/results/taxonomy_long
"""

gwf.target("make_metadata", inputs = ['/home/mpavila/eDNA/faststorage/Velux/CoastSequence/Autumn/Bac16s/base/results/DADA2_nochim.otus','/home/mpavila/eDNA/faststorage/Velux/CoastSequence/Autumn/Bac16s/both_silva/results/metadata/extraction_refs_both_seasons.txt','/home/mpavila/eDNA/faststorage/Velux/CoastSequence/Autumn/Bac16s/both_silva/results/metadata/PSU_refs_both.txt'], outputs = ['/home/mpavila/eDNA/faststorage/Velux/CoastSequence/Autumn/Bac16s/both_silva/results/metadata/metadata_both.txt','/home/mpavila/eDNA/faststorage/Velux/CoastSequence/Autumn/Bac16s/both_silva/results/otu_silva.txt','/home/mpavila/eDNA/faststorage/Velux/CoastSequence/Autumn/Bac16s/both_silva/results/tax_silva.txt'],
           cores=1, memory='256g', walltime="24:00:01") << """
source ~/miniconda3/etc/profile.d/conda.sh
conda activate aldex
Rscript /home/mpavila/eDNA/faststorage/Velux/CoastSequence/Autumn/Bac16s/both_silva/scripts/make_metadata.R
"""

gwf.target("cleanup", inputs = ['/home/mpavila/eDNA/faststorage/Velux/CoastSequence/Autumn/Bac16s/both_silva/results/otu_silva.txt'], outputs = ['/home/mpavila/eDNA/faststorage/Velux/CoastSequence/Autumn/Bac16s/both_silva/results/metadata/merged_seq_run_cleaned_metadata_ASV_wise.txt','/home/mpavila/eDNA/faststorage/Velux/CoastSequence/Autumn/Bac16s/both_silva/results/cleaned_otu_table_ASV_wise.txt','/home/mpavila/eDNA/faststorage/Velux/CoastSequence/Autumn/Bac16s/both_silva/results/cleanup_ASV_wise/summary/summary_clean_up_both_seasons.txt','/home/mpavila/eDNA/faststorage/Velux/CoastSequence/Autumn/Bac16s/both_silva/results/cleanup_ASV_wise/summary/summary_clean_up_autumn.txt','/home/mpavila/eDNA/faststorage/Velux/CoastSequence/Autumn/Bac16s/both_silva/results/cleanup_ASV_wise/summary/summary_clean_up_spring.txt'],
           cores=1, memory='500g', walltime="10:00:01") << """
source ~/miniconda3/etc/profile.d/conda.sh
conda activate aldex
Rscript /home/mpavila/eDNA/faststorage/Velux/CoastSequence/Autumn/Bac16s/both_silva/scripts/clean_up_ASV_wise.R
"""

gwf.target("nosing", inputs = ['/home/mpavila/eDNA/faststorage/Velux/CoastSequence/Autumn/Bac16s/both_silva/results/metadata/merged_seq_run_cleaned_metadata_ASV_wise.txt','/home/mpavila/eDNA/faststorage/Velux/CoastSequence/Autumn/Bac16s/both_silva/results/cleaned_otu_table_ASV_wise.txt'], outputs = ['/home/mpavila/eDNA/faststorage/Velux/CoastSequence/Autumn/Bac16s/both_silva/results/metadata/no_control_no_sing_samples_cleaned_metadata_ASV_wise.txt','/home/mpavila/eDNA/faststorage/Velux/CoastSequence/Autumn/Bac16s/both_silva/results/no_sing_cleaned_otu_table_ASV_wise.txt'],
           cores=1, memory='256g', walltime="10:00:01") << """
source ~/miniconda3/etc/profile.d/conda.sh
conda activate aldex
Rscript /home/mpavila/eDNA/faststorage/Velux/CoastSequence/Autumn/Bac16s/both_silva/scripts/no_sing_ASV_wise.R
"""

gwf.target("cleantax", inputs = ['/home/mpavila/eDNA/faststorage/Velux/CoastSequence/Autumn/Bac16s/both_silva/results/metadata/no_control_no_sing_samples_cleaned_metadata_ASV_wise.txt','/home/mpavila/eDNA/faststorage/Velux/CoastSequence/Autumn/Bac16s/both_silva/results/no_sing_cleaned_otu_table_ASV_wise.txt','/home/mpavila/eDNA/faststorage/Velux/CoastSequence/Autumn/Bac16s/both_silva/results/tax_silva.txt'], outputs = ['/home/mpavila/eDNA/faststorage/Velux/CoastSequence/Autumn/Bac16s/both_silva/results/metadata/cleaned_silva_metadata.txt','/home/mpavila/eDNA/faststorage/Velux/CoastSequence/Autumn/Bac16s/both_silva/results/cleaned_otu_silva.txt','/home/mpavila/eDNA/faststorage/Velux/CoastSequence/Autumn/Bac16s/both_silva/results/cleaned_tax_silva.txt'],
           cores=1, memory='256g', walltime="4:00:01") << """
source ~/miniconda3/etc/profile.d/conda.sh
conda activate aldex
Rscript /home/mpavila/eDNA/faststorage/Velux/CoastSequence/Autumn/Bac16s/both_silva/scripts/clean_tax_ASV_silva.R
"""

gwf.target("depth_cutoff", inputs = ['/home/mpavila/eDNA/faststorage/Velux/CoastSequence/Autumn/Bac16s/both_silva/results/metadata/cleaned_silva_metadata.txt','/home/mpavila/eDNA/faststorage/Velux/CoastSequence/Autumn/Bac16s/both_silva/results/cleaned_otu_silva.txt','/home/mpavila/eDNA/faststorage/Velux/CoastSequence/Autumn/Bac16s/both_silva/results/cleaned_tax_silva.txt'], outputs = ['/home/mpavila/eDNA/faststorage/Velux/CoastSequence/Autumn/Bac16s/both_silva/results/metadata/f_silva_metadata.txt.txt','/home/mpavila/eDNA/faststorage/Velux/CoastSequence/Autumn/Bac16s/both_silva/results/f_oyu_silva.txt.txt','/home/mpavila/eDNA/faststorage/Velux/CoastSequence/Autumn/Bac16s/both_silva/results/f_tax_silva.txt.txt'],
           cores=1, memory='256g', walltime="4:00:01") << """
source ~/miniconda3/etc/profile.d/conda.sh
conda activate aldex
Rscript /home/mpavila/eDNA/faststorage/Velux/CoastSequence/Autumn/Bac16s/both_silva/scripts/find_depth_cutoff_silva.R
"""

