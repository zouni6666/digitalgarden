---
{"dg-publish":true,"permalink":"/02 本科毕设/benchmark/benchmark/","tags":["gardenEntry"]}
---


# 1 benchmark
Benchmarker 是一种无偏、数据驱动的评估基因优先考虑方法的方法。基于靠近致病基因的 SNP 应该富集于表型遗传性的假设，Benchmarker 使用分层 LD 得分回归（S-LDSC）56来估计靠近优先基因的 SNP 对每个 SNP 遗传性的平均贡献。使用 S-LDSC，Benchmarker 在“基线模型”中联合建模与优先基因对应的 SNP 注释以及其他53个注释，其中包括基因、调控和保守区域。为了评估性能，我们使用回归系数τ及其 P 值来检验假设τ > 0；τ衡量了在控制基线注释后，靠近优先基因的 SNP 对每个 SNP 遗传性的贡献。为了使τ在不同表型之间可比，我们通过每个表型的平均每个 SNP 遗传性对τ进行了归一化处理，并将该量称为归一化τ。对于我们的分析，我们选择了每个表型中 PoP 得分最高的500个基因作为优先基因集，并在每个基因的转录起始位点两侧使用100 kb 的窗口将 SNP 映射到基因。 [https://github.com/RebeccaFine/benchmarker?tab=readme-ov-file](https://github.com/RebeccaFine/benchmarker?tab=readme-ov-file)

先跳过
```shell
sum_file=AG.fastGWA
trait=AG
n_sample=517

# 1 
cat ../$sum_file | awk 'BEGIN{OFS=FS="\t"}{print $2,$1,$3}'  > ${trait}_snp.loc
cat ../$sum_file | awk 'BEGIN{OFS=FS="\t"}{print $2,$10}'  > ${trait}_p_value.loc

magma \
    --annotate \
    --snp-loc ${trait}_snp.loc \
    --gene-loc ../data/depict_noMHC.gene.loc \
    --out $trait

# 2 
magma \
--bfile /home/fenglijun/data3/project/data/human/reference/panel/g1000_eur \
--pval ${trait}_p_value.loc  N=$n_sample --gene-annot bmi_magma.genes.annot --out bmi_magma


# 3 
python 2_withhold_chromosome.py bmi_magma


# 4 
magma \
    --gene-results ${trait}_noChr${chrom}.genes.raw \
    --gene-covar ../data/GPL570-GPL96-GPL1261-GPL1355TermGeneZScores-MGI_MF_CC_RT_IW_BP_KEGG_z_z_GTExGenesOnly.txt onesided=greater \
    --out ${trait}_noChr${chrom}
 
   
trait=bmi_magma 
python ../src/prioritize_genes_from_enrichment_results.py \
--output_directory . \
--trait ${trait} \
--results_file ${trait}_magma_noChr@.gcov.out \
--gene_boundary_file ../data/GPL570ProbeENSGInfo+HGNC_reformatted_noMHC_depictAndGtexGeneIntersection_RF.txt \
--set_definitions_file ../data/Top50GenesPerGeneSet_DEPICTGenes_Ensembl_depictAndGtexGeneIntersection.txt \
--set_definitions_label Top50GenesPerGeneSet \
--output_label magma \
--percentage_cutoff 10 \
--results_file_set_col COVAR \
--results_file_col_to_sort_on P \ 
--sort_direction ascending \
--results_file_separator '\s+'
```
先做 leadsnp 数量前 20 的性状
1. 先对 leadsnp 数量前 20 的性状进行 [[02 本科毕设/benchmark/benchmark\|02 本科毕设/benchmark/benchmark]] :
/home/data 6/zhanghao/pathway/pathway/pigGWAS/lead_snp/unique_study_symbol_counts. csv

2. 选蛋白质编码的基因
取 leadsnp 位点上下游各500 Kb 的范围找基因（在 gtf 文件中带有 ptotein_id 的）
构建金标准数据集使用 leadsnnp 文件，如果 leadsnp 在那个 exon 即与 uniq_gene_exon. txt 文件有交集，基因是有效的属于正基因集

# golden standard
描述了他们如何构建一个新的评估集来比较基于相似性和基因座的基因优先级方法。
为了构建新的评估集，研究人员首先利用他们最近对来自UK Biobank的95个特征进行的统计精细定位结果，
1 将可能的致病基因定义为携带精细定位（PIP > 0.5）的蛋白编码变异体。
2 然后，根据特征进行匹配，他们在这些蛋白编码基因的500 kb内识别出独立的非编码可信区间。这种方法识别出1,348个非编码可信区间，其中包含物理上相邻但独立的编码变异体信号。
3 研究人员从这1,348个基因座中创建了一个评估集，其中包括与基因座定义的非编码可信区间在500 kb内的所有基因（每个基因座中位数为13个基因）。对于具有相同特征的独立、精细定位编码变异体的基因座，将其标记为阳性，而没有这种变异体的基因则标记为阴性。这种分配直接编码了研究人员的假设：（1）携带精细定位编码变异体的基因与特征相关，（2）一个基因座中的多个独立关联很可能通过同一个基因发挥作用。
通过使用这个新的评估集，研究人员比较了基于相似性的 PoPS 方法与其他基于相似性的方法之间的性能差异。他们发现，与 DEPICT、NetWAS 和 MAGMA-sim 相比，PoPS 在使用非编码遗传信号时具有更高的召回率和更高的精确度。此外，研究人员还展示了使用完整的非编码遗传数据时，PoPS 比在 LOCO 框架中运行的 PoPS 具有更高的召回率和更高的精确度。
```shell
# 获取外显子bed文件 即有效的基因列表和位置（只有在外显子上的基因才有效 而且要是蛋白质编码的）这个是只含有外显子的  q是蛋白编码基因的列表
zcat /home/data6/zhanghao/pathway/pathway/pigGWAS/lead_snp/Sus_scrofa.Sscrofa11.1.100.gtf.gz |grep --file /home/data6/zhanghao/pathway/pathway/pigGWAS/lead_snp/q |awk 'BEGIN{OFS=FS="\t"}{if ($3=="exon")print $0}' |cut -f 1,4,5,9 |sort |uniq > /home/data6/zhanghao/pathway/pathway/pigGWAS/lead_snp/uniq_gene_exon.txt


# 1 构建金标准数据集 提取不同性状的含有lead_snp的位点
trait=M_TNUM
cat /home/data6/zhanghao/pathway/pathway/pigGWAS/lead_snp/PigBiobank_lead_SNP_in_300_studies.csv |awk 'BEGIN{OFS=FS=","}{print $4,$5,$5,$1}' |sed 's/,/\t/g' |grep "${trait}" |sed '1d' >fine_mapping_all_snp.bed

#提取在外显子中且能转录的基因
awk 'BEGIN{FS=OFS="\t"} {gsub(/gene_id "/, "", $4); gsub(/".*$/, "", $4); print}' uniq_gene_exon.txt |uniq > uniq_gene_exon2.txt

# 2 获取正例基因及其基因座
# # 构建金标准数据集使用 lead_snp 文件，如果 lead_snp 在那个 exon 即与 uniq_gene_exon.txt 文件有交集，基因是有效的属于正基因集
bedtools intersect -a fine_mapping_all_snp.bed -b /home/data6/zhanghao/pathway/pathway/pigGWAS/lead_snp/uniq_gene_exon2.txt -wa -wb |cut -f8|sort|uniq >positive_gene.list
#有效基因的位点
bedtools intersect -a fine_mapping_all_snp.bed -b /home/data6/zhanghao/pathway/pathway/pigGWAS/lead_snp/uniq_gene_exon2.txt -wa -wb |cut -f 5,6,7,8|sort|uniq|awk 'BEGIN{OFS=FS="\t"}{print $1,$2-500000,$3+500000,$4}' >golden_locus.bed

# 获取全部基因座 -b中的文件应该是所有的基因
#**获取全部蛋白编码基因位置**
zcat /home/data6/zhanghao/pathway/pathway/pigGWAS/lead_snp/Sus_scrofa.Sscrofa11.1.100.gtf.gz |grep --file /home/data6/zhanghao/pathway/pathway/pigGWAS/lead_snp/q |awk 'BEGIN{OFS=FS="\t"} {print $0}' |cut -f 1,4,5,9 |sort |uniq > uniq_gene_all.bed

# 3 获取基因座上包含正例和负例的所有基因
bedtools intersect -a golden_locus.bed -b /home/data6/zhanghao/pathway/pathway/pigGWAS/lead_snp/uniq_gene_all2.bed -wa -wb |cut -f8  |sort|uniq >all_golden_gene_in_locus.txt
wc all_golden_gene_in_locus.txt

# 平均每个位点所含基因数量
  bedtools intersect -a fine_mapping_all_snp.bed -b uniq_gene_all2.bed -wa -wb |sort |uniq > gene_snp.txt

awk 'BEGIN{OFS=FS="\t"}{print $1,$2-500000,$2+500000}' fine_mapping_all_snp.bed | awk 'BEGIN{OFS=FS="\t"}{if($2<0) print $1,0,$3; else print $0}'  >q

bedtools  intersect -a q -b sus11_protein_gene.bed  -wa -wb|wc

```
计算精确率和召回率
真正阳性（TP）是一个被优先考虑的基因，它是条件阳性的；假阳性（FP）是一个被优先考虑的基因，它是条件阴性的；真正阴性（TN）是一个未被优先考虑且条件阴性的基因；假阴性（FN）是一个未被优先考虑且条件阳性的基因。我们的两个问题的答案分别由精确度和召回率给出，

精确度为TP的数量/（TP的数量+FP的数量），召回率为TP的数量/（TP的数量+FN的数量）。

简而言之，研究人员使用他们构建的评估基因集来评估每种方法的精确度和召回率。精确度表示当一个基因被优先考虑时，它确实与条件相关的置信度有多高。召回率表示在所有真正相关的基因中，方法优先考虑了多少个。他们根据TP、FP、TN和FN的数量计算了精确度和召回率。

```shell
# magma

l ../../5_magma/t10.genes.out |sed 's/  */\t/g'|  awk  'BEGIN{OFS=FS="\t"}{if($9<0.05) print $0}'|sort -k9g| cut -f1 >gene_magma.txt

rm q
python 	/home/data3/fenglijun/script/post_gwas/presion_recall.py    all_golden_gene_in_locus.txt  positive_gene.list  ../5_magma/magma_prior_gene.txt  magma >>q 2>w
python 	/home/data3/fenglijun/script/post_gwas/presion_recall.py    all_golden_gene_in_locus.txt  positive_gene.list  ../5_magma/geneset_prior_gene.txt  geneset >>q 2>w
#python 	/home/data3/fenglijun/script/post_gwas/presion_recall.py    all_golden_gene_in_locus.txt  positive_gene.list  ../5_magma/geneset_prior_gene.txt  DEPICT
python 	/home/data3/fenglijun/script/post_gwas/presion_recall.py    all_golden_gene_in_locus.txt  positive_gene.list   ../2_twas/twas_prior_gene.txt  twas >>q 2>w 
python 	/home/data3/fenglijun/script/post_gwas/presion_recall.py    all_golden_gene_in_locus.txt  positive_gene.list  ../3_smr/smr_prior_gene.txt  smr >>q 2>w
python 	/home/data3/fenglijun/script/post_gwas/presion_recall.py    all_golden_gene_in_locus.txt  positive_gene.list  ../4_coloc/coloc_prior_gene.txt  coloc >>q 2>w
python 	/home/data3/fenglijun/script/post_gwas/presion_recall.py    all_golden_gene_in_locus.txt  positive_gene.list  ../7_epi/abc_prior_gene.txt  abc >>q 2>w
#python 	/home/data3/fenglijun/script/post_gwas/presion_recall.py    all_golden_gene_in_locus.txt  positive_gene.list  ../5_magma/geneset_prior_gene.txt  pops

cut -f1 q
echo '======================'
cut -f6 q
echo '======================'
cut -f7 q

```

# 3 Closest gene enrichment 
使用了正态近似来近似一个测试统计量c的零分布。测试统计量c表示在一个位点中，被PoPS优先考虑的基因数以及最接近该位点主导变异体的基因。在零假设下，PoPS随机选择该位点中最近的基因，概率为1/nl，其中nl是该位点中的基因数。在所有L个位点上，c的分布被近似为一组具有不同偏倚的独立伯努利分布的总和。为了计算方便，当L很大时，我们使用了一个匹配矩的正态分布来近似这个分布。
在零假设下，我们对 c > ∑ 1∶L 1/nl 进行了单侧检验，即检验c是否大于所有位点中的基因数之和的期望值。此外，我们计算了被PoPS优先考虑的基因中最近基因的富集程度，通过计算观察到的富集比率c/∑ 1∶L 1/nl，并估计了这个富集的标准误差。为了估计富集的标准误差，我们使用了自助法（bootstrap），对每个表型的L个位点进行了1,024次自助重复采样。

```shell
# 在我们的问题中，进行的是右尾检验，因为我们对 c > ∑(1/n_l) 进行检验，即我们关注的是观察到的c值是否大于预期值

trait=t_10

rm closest_gene.txt
rm loc_gene_num.txt

ca py310
mkdir -p 2_cloest
cd 2_cloest
while FS=OFS='\t' read -r col1 col2 col3; do
  echo "process locus_${col1}_${col2}_${col3}"
  echo -e "${col1}\t${col2}\t${col3}\ta" > sig_loc_bed.txt
	a=${col1}_${col2}_${col3}
	 l ../../1_finemapping/${trait}_locus_cs_snp_${a} |head -n 1 |awk  'BEGIN{OFS=FS="\t"}{print $1,$3,$3,$2}' >snp.bed
	 bedtools closest -a snp.bed -b  /home/fenglijun/data3/project/data/human/reference/genome/g37/sorted_genes.bed |cut -f8 >> closest_gene.txt
	 bedtools intersect -a sig_loc_bed.txt -b /home/fenglijun/data3/project/data/human/reference/genome/g37/uniq_gene_15436.bed | wc |sed 's/  */\t/g'|cut -f2 >>loc_gene_num.txt
done < ../../4_coloc/locus.txt

# 打印出p_value enrichment se
rm q
python  /home/fenglijun/data3/script/post_gwas/cloest_gene_enrichment.py closest_gene.txt loc_gene_num.txt ../../5_magma/magma_prior_gene.txt magma >>q 2>w
python  /home/fenglijun/data3/script/post_gwas/cloest_gene_enrichment.py closest_gene.txt loc_gene_num.txt  ../../5_magma/geneset_prior_gene.txt  geneset >>q 2>w
python  /home/fenglijun/data3/script/post_gwas/cloest_gene_enrichment.py closest_gene.txt loc_gene_num.txt  ../../2_twas/twas_prior_gene.txt twas  >>q 2>w
python  /home/fenglijun/data3/script/post_gwas/cloest_gene_enrichment.py closest_gene.txt loc_gene_num.txt ../../3_smr/smr_prior_gene.txt  smr >>q 2>w
python  /home/fenglijun/data3/script/post_gwas/cloest_gene_enrichment.py closest_gene.txt loc_gene_num.txt  ../../4_coloc/coloc_prior_gene.txt  coloc >>q 2>w
python  /home/fenglijun/data3/script/post_gwas/cloest_gene_enrichment.py closest_gene.txt loc_gene_num.txt ../../7_epi/abc_prior_gene.txt  abc >>q 2>w

cut -f1 q
echo '======================'
cut -f2 q
echo '======================'
cut -f3 q
echo '======================'
cut -f4 q

```
## 新
```shell
cd /home/data6/zhanghao/pathway/pathway/pigGWAS/benchmark/golden_standard/smr/

locus.txt是locus_${trait}.txt：lead_snp上下游500kb

sorted_genes_bed：按照位置排序后的全部蛋白编码的基因位置 sus11_protein_gene.bed

${trait}_locus_cs_snp  /home/data6/zhanghao/pathway/pathway/pigGWAS/benchmark/golden_standard/smr/fine_mapping_all_snp_D_ADG.bed:每个性状的lead_snp位点

uniq_gene_15436.bed  /home/data6/zhanghao/pathway/pathway/pigGWAS/lead_snp/re/sus11_protein_gene.bed：全部蛋白编码的基因位置

能转录为蛋白质的（即gtf文件中有蛋白编号的）一定在外显子区域；在外显子区域的基因不一定能转录为蛋白质


${trait}_locus_cs_snp_${a}在finemapping中能找到  （重新做finemapping） 看less -SN /home/data6/zhanghao/pathway/pathway/pigGWAS/finemapping/AG_FINAL/AG_locus_cs_snp_15_1178484_2178484中含义

```
```shell

trait=

rm closest_gene.txt
rm loc_gene_num.txt

mkdir -p 2_cloest
cd 2_cloest

while FS=OFS='\t' read -r col1 col2 col3; do
  echo "process locus_${col1}_${col2}_${col3}"
  echo -e "${col1}\t${col2}\t${col3}\ta" > sig_loc_bed.txt
	a=${col1}_${col2}_${col3}
	 l ../../1_finemapping/${trait}_locus_cs_snp_${a} |head -n 1 |awk  'BEGIN{OFS=FS="\t"}{print $1,$3,$3,$2}' >snp.bed
	 bedtools closest -a snp.bed -b  /home/fenglijun/data3/project/data/human/reference/genome/g37/sorted_genes.bed |cut -f8 >> closest_gene.txt
	 bedtools intersect -a sig_loc_bed.txt -b /home/fenglijun/data3/project/data/human/reference/genome/g37/uniq_gene_15436.bed | wc |sed 's/  */\t/g'|cut -f2 >>loc_gene_num.txt
done < ../../4_coloc/locus.txt


```


