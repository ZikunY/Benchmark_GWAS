#!/bin/sh
#$ -cwd -S /bin/bash
#$ -l h_vmem=10G
#$ -l h_rt=02:00:00
cd /mnt/mfs/hgrcgrid/homes/zy2412/1000G/data
cd /mnt/mfs/hgrcgrid/shared/GT_ADMIX/Zikun_peru_project/HAPNEST/office
module load PLINK/2.0
module load PLINK/1.9

module load RFMix/2.0.3
module load BCFTOOLS/1.9
 module load Shapeit/2.r837
conda activate my-env

 bcftools convert /mnt/mfs/hgrcgrid/shared/GT_ADMIX/PERU/TOPMED/chr20.dose.vcf.gz --haplegendsample  peru_chr20
 



#################Old file##############
plink2 --vcf /mnt/mfs/hgrcgrid/homes/zy2412/1000G/data/1kGP_high_coverage_Illumina.chr20.filtered.SNV_INDEL_SV_phased_panel.vcf.gz --chr 20 --keep EAS_GBR_iid.txt --recode vcf --snps-only -out /mnt/mfs/hgrcgrid/shared/GT_ADMIX/Zikun_peru_project/HAPNEST/office/EAS_GBR_chr20


plink2 --vcf /mnt/mfs/hgrcgrid/homes/zy2412/1000G/data/1kGP_high_coverage_Illumina.chr20.filtered.SNV_INDEL_SV_phased_panel.vcf.gz --chr 20 --keep /mnt/mfs/hgrcgrid/homes/zy2412/1000G/data/JPT_GBR_YRI_iid.txt --recode vcf --snps-only -out /mnt/mfs/hgrcgrid/shared/GT_ADMIX/Zikun_peru_project/HAPNEST/office/chr_20_JPT_GBR_YRI_plink2

plink2 --vcf /mnt/mfs/hgrcgrid/homes/zy2412/1000G/data/1kGP_high_coverage_Illumina.chr21.filtered.SNV_INDEL_SV_phased_panel.vcf.gz --chr 21 --keep /mnt/mfs/hgrcgrid/homes/zy2412/1000G/data/JPT_GBR_YRI_iid.txt --recode vcf --snps-only -out /mnt/mfs/hgrcgrid/shared/GT_ADMIX/Zikun_peru_project/HAPNEST/office/chr_21_JPT_GBR_YRI_plink2

plink2 --vcf /mnt/mfs/hgrcgrid/shared/GT_ADMIX/Zikun_peru_project/HAPNEST/office/hapnest_20.vcf --chr 20 --indep-pairwise 1500 150 0.1 -out plin2_chr20

plink2 --vcf /mnt/mfs/hgrcgrid/shared/GT_ADMIX/Zikun_peru_project/HAPNEST/office/hapnest_21.vcf --chr 21 --indep-pairwise 1500 150 0.1 -out plin2_chr21

plink2 --vcf /mnt/mfs/hgrcgrid/shared/GT_ADMIX/Zikun_peru_project/HAPNEST/office/hapnest_20.vcf --chr 20 --extract /mnt/mfs/hgrcgrid/shared/GT_ADMIX/Zikun_peru_project/HAPNEST/office/plin2_chr20.prune.in --recode vcf -out hapnest_20_prune

plink2 --vcf /mnt/mfs/hgrcgrid/shared/GT_ADMIX/Zikun_peru_project/HAPNEST/office/hapnest_21.vcf --chr 21 --extract /mnt/mfs/hgrcgrid/shared/GT_ADMIX/Zikun_peru_project/HAPNEST/office/plin2_chr21.prune.in --recode vcf -out hapnest_21_prune

plink2 --vcf hapnest_20_prune.vcf --chr 20 --keep Hapnest_AMR_id_shorten.txt --recode vcf -out hapnest_20_prune_shorten

rfmix -f /mnt/mfs/hgrcgrid/shared/GT_ADMIX/Zikun_peru_project/Simulation/Peru_genotypes_seq_4000_clean.vcf \
        -r /mnt/mfs/hgrcgrid/homes/zy2412/1000G/data/PEL/chr_1_JPT_GBR_plink2.vcf \
        -m /mnt/mfs/hgrcgrid/homes/zy2412/1000G/data/PEL/JPT_GBR_195.tsv \
        -g /mnt/mfs/hgrcgrid/homes/zy2412/1000G/data/PEL/genetic_map_chr1_modified_b38.txt \
        -o /mnt/mfs/hgrcgrid/shared/GT_ADMIX/Zikun_peru_project/Simulation/Peru_genotypes_seq_4000_clean \
        --chromosome=1
        
rfmix -f hapnest_20_prune_shorten.vcf \
        -r chr_20_JPT_GBR_YRI_plink2.vcf \
        -m JPT_GBR_YRI.txt \
        -g /mnt/mfs/hgrcgrid/homes/zy2412/1000G/data/genetic_map_files_38/genetic_map_modified_chr20_b38.txt \
        -o hapnest_three_cohorts_clean_chr20 \
        --chromosome=20

rfmix -f hapnest_21_prune.vcf \
        -r chr_21_JPT_GBR_YRI_plink2.vcf \
        -m JPT_GBR_YRI.tsv \
        -g /mnt/mfs/hgrcgrid/homes/zy2412/1000G/data/genetic_map_files_38/genetic_map_modified_chr21_b38.txt \
        -o hapnest_three_cohorts_clean_chr21 \
        --chromosome=21
        
rfmix -f vcf/hapnest_20_prune_chunk_two_way_1.vcf \
        -r EAS_GBR_chr20.vcf \
        -m EAS_GBR.txt \
        -g /mnt/mfs/hgrcgrid/homes/zy2412/1000G/data/genetic_map_files_38/genetic_map_modified_chr20_b38.txt \
        -o RFMix/hapnest_two_cohorts_prune_chr20_EAS_EUR_${CHRARG} \
        --chromosome=20

        
###########Run null model scenario########

plink2 --vcf vcf/hapnest_20_prune_chunk_two_way_1.vcf --make-bed --chr 20 -out hapnest_20_prune_chunk_two_way_1

plink2 --vcf vcf/hapnest_20_prune_chunk_two_way_1.vcf --pca --make-rel --chr 20 -out vcf/hapnest_20_prune_chunk_two_way_1

########Saige##########
bcftools view -O z -o vcf/hapnest_20_prune_chunk_two_way_1_phased.vcf.gz vcf/hapnest_20_prune_chunk_two_way_1_phased.vcf

 bcftools index -t vcf/hapnest_20_prune_chunk_two_way_1.vcf.gz
 
 ######Tractor###########
 
 bcftools convert ALL.chr20.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz --haplegendsample  1000G_chr20_bcf_v2a
 
 1kGP_high_coverage_Illumina.chr20.filtered.SNV_INDEL_SV_phased_panel.vcf.gz
 
 plink2 --vcf  1kGP_high_coverage_Illumina.chr20.filtered.SNV_INDEL_SV_phased_panel.vcf.gz --export haps --chr 20 -out EAS_GBR_chr20_hap

 plink2 --vcf  1kGP_high_coverage_Illumina.chr20.filtered.SNV_INDEL_SV_phased_panel.vcf.gz --export hapslegend --chr 20 -out EAS_GBR_chr20_hap


 plink2 --vcf  EAS_GBR_chr20.vcf --export hapslegend --chr 20 -out EAS_GBR_chr20_hapslegend


  plink2 --vcf  1kGP_high_coverage_Illumina.chr20.filtered.SNV_INDEL_SV_phased_panel.vcf.gz --export hapslegend --chr 20 -out 1000G_chr20_hap
 
 plink2 --vcf /mnt/mfs/hgrcgrid/shared/GT_ADMIX/Zikun_peru_project/HAPNEST/office/vcf/hapnest_20_prune_chunk_two_way_1.vcf --chr 20 --extract /mnt/mfs/hgrcgrid/shared/GT_ADMIX/Zikun_peru_project/HAPNEST/office/shapeit_prune.in --recode vcf -out hapnest_20_prune_chunk_two_way_1_algined

 
     ######1KG
 shapeit -check \
        --input-vcf vcf/hapnest_20_prune_chunk_two_way_1.vcf \
        --input-map /mnt/mfs/hgrcgrid/homes/zy2412/1000G/data/genetic_map_files_38/genetic_map_combined_chr20_b38.txt \
        --input-ref 1000G_chr20_bcf.hap.gz 1000G_chr20_bcf.legend.gz 1000G_chr20_bcf.samples \
        --output-log alignments_bcf
 

 shapeit --input-vcf vcf/hapnest_20_prune_chunk_two_way_1.vcf \
         --input-map /mnt/mfs/hgrcgrid/homes/zy2412/1000G/data/genetic_map_files_38/genetic_map_combined_chr20_b38.txt \
        --input-ref 1000G_chr20_bcf.hap.gz 1000G_chr20_bcf.legend.gz 1000G_chr20_bcf.samples \
        -O vcf/hapnest_20_prune_chunk_two_way_1_phased \
        --exclude-snp alignments_bcf.snp.strand.exclude
###### V2a
   shapeit -check \
        --input-vcf vcf/hapnest_20_prune_chunk_two_way_1.vcf \
        --input-map /mnt/mfs/hgrcgrid/homes/zy2412/1000G/data/genetic_map_files_38/genetic_map_combined_chr20_b38.txt \
        --input-ref 1000G_chr20_bcf_v2a.hap.gz 1000G_chr20_bcf_v2a.legend.gz 1000G_chr20_bcf_v2a.samples \
        --output-log alignments_bcf_v2a
 

 shapeit --input-vcf vcf/hapnest_20_prune_chunk_two_way_1.vcf \
        --input-map /mnt/mfs/hgrcgrid/homes/zy2412/1000G/data/genetic_map_files_38/genetic_map_combined_chr20_b38.txt \
        --input-ref 1000G_chr20_bcf_v2a.hap.gz 1000G_chr20_bcf_v2a.legend.gz 1000G_chr20_bcf_v2a.samples \
        -O vcf/hapnest_20_prune_chunk_two_way_1_phased_v2a \
        --exclude-snp alignments_bcf_v2a.snp.strand.exclude

  
  #####Shapeit convert



shapeit -convert \
        --input-haps vcf/hapnest_20_prune_chunk_two_way_1_phased_quick \
        --output-vcf vcf/hapnest_20_prune_chunk_two_way_1_phased_quick.vcf

 rfmix -f vcf/hapnest_20_prune_chunk_two_way_1_phased_quick.vcf \
        -r EAS_GBR_chr20.vcf \
        -m EAS_GBR.txt \
        -g /mnt/mfs/hgrcgrid/homes/zy2412/1000G/data/genetic_map_files_38/genetic_map_modified_chr20_b38.txt \
        -o RFMix/hapnest_two_cohorts_prune_chr20_EAS_EUR_phased_quick \
        --chromosome=20
       
######Tracts
python /mnt/mfs/hgrcgrid/homes/zy2412/1000G/data/tractor_example/Tractor/ExtractTracts.py \
      --msp RFMix/hapnest_two_cohorts_prune_chr20_EAS_EUR_phased_quick \
      --vcf vcf/hapnest_20_prune_chunk_two_way_1_phased_quick \
      --num-ancs 2


 shapeit -check \
        --input-vcf vcf/hapnest_20_prune_chunk_two_way_1_algined.vcf \
        --input-map /mnt/mfs/hgrcgrid/homes/zy2412/1000G/data/PEL/genetic_map_chr1_modified_b38.txt \
        --input-ref 1000G_chr20_hap.haps 1000G_chr20_hap.legend 1000G_chr20_hap.sample \
        --output-log alignments
        


        
 shapeit --input-vcf vcf/hapnest_20_prune_chunk_two_way_1.vcf.gz \
         --input-map /mnt/mfs/hgrcgrid/homes/zy2412/1000G/data/PEL/genetic_map_chr1_modified_b38.txt \
         --input-ref EAS_GBR_chr20_hapslegend.haps EAS_GBR_chr20_hapslegend.legend EAS_GBR_chr20_hapslegend.sample \
         -O vcf/hapnest_20_prune_chunk_two_way_1_phased \
         --exclude-snp alignments.snp.strand.exclude
 
python /mnt/mfs/hgrcgrid/homes/zy2412/1000G/data/tractor_example/Tractor/ExtractTracts.py \
      --msp RFMix/hapnest_two_cohorts_prune_chr20_EAS_EUR_1 \
      --vcf vcf/hapnest_20_prune_chunk_two_way_1 \
      --zipped \
      --num-ancs 2

python /mnt/mfs/hgrcgrid/homes/zy2412/1000G/data/tractor_example/Tractor/ExtractTracts.py \
      --msp /mnt/mfs/hgrcgrid/shared/GT_ADMIX/Zikun_peru_project/Simulation/Peru_genotypes_seq_4000_clean \
      --vcf /mnt/mfs/hgrcgrid/shared/GT_ADMIX/Zikun_peru_project/Simulation/Peru_genotypes_seq_4000_clean \
      --num-ancs 2

for asnp in {'ultra.rare.snp','rare.snp','uncommon.snp','common.snp'}; do
for ad in {1.5,2,2.5,3}; do
mkdir -p /ifs/scratch/msph/eigen/zy2412/GT_admix/Simulation/office/power_phenotype/${asnp}/${ad}
done
done
