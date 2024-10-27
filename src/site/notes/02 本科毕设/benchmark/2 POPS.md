---
{"dg-publish":true,"permalink":"/02 本科毕设/benchmark/2 POPS/"}
---


[GitHub - FinucaneLab/pops](https://github.com/FinucaneLab/pops)
包括基因密度、有效基因大小、等位基因计数等多种协变量。
```shell
# pops full feature
/home/fenglijun/data3/project/data/human/pops_feature_results/pops_gene_feature/PoPS.features.txt
/home/fenglijun/data3/project/data/human/pops_feature_results/pops_gene_feature/feature_100.txt

# 特征选取
python /home/fenglijun/data3/project/2_gene_geneset/1_pops/munge_feature_directory.py \
 --gene_annot_path /home/fenglijun/data3/project/2_gene_geneset/1_pops/example/data/utils/pops_15436_tss_annot.txt \
 --feature_dir /home/data3/zhanghao/pathway/data/pig/PigGTEx_v0.Gene.TPM_matrix \
 --save_prefix exampa \
 --max_cols 500
 
 # 程序运行
 python /home/fenglijun/data3/project/2_gene_geneset/1_pops/pops.py \
 --gene_annot_path example/data/utils/gene_annot_jun10.txt \
 --feature_mat_prefix example/data/features_munged/pops_features \
 --num_feature_chunks 2 \
 --magma_prefix example/data/magma_scores/PASS_Schizophrenia \
 --control_features_path example/data/utils/features_jul17_control.txt \
 --out_prefix example/out/PASS_Schizophrenia
 
  python /home/fenglijun/data3/project/2_gene_geneset/1_pops/pops.py \
 --gene_annot_path ../pops_15436_tss_annot.txt \
 --feature_mat_prefix /home/fenglijun/data3/project/2_gene_geneset/1_pops/example/data/utils/q/exampa \
 --num_feature_chunks 2 \
 --magma_prefix ../../../../../benchmark/9/5_magma/t9 \
 --control_features_path ../features_jul17_control.txt \
 --out_prefix ./out/test
```