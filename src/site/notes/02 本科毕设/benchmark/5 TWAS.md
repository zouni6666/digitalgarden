---
{"dg-publish":true,"permalink":"/02 本科毕设/benchmark/5 TWAS/"}
---


[TWAS / FUSION](http://gusevlab.org/projects/fusion/) 
**TWAS（全转录组关联分析）是把转录调控（expression）作为遗传变异（genotype）和表型（phenotype）之间的中介，将单个遗传变异与表型的关联转换成基因/转录本与表型的关联。**
TWAS依赖于两个主要数据来源：
1. **全基因组关联研究（GWAS）数据**：这些研究通常会鉴定出特定的SNP（单核苷酸多态性）与疾病风险之间的统计关联。
2. **基因表达数据**：通常来源于大型生物信息数据库，如GTEx（基因组组织表达项目），这些数据显示不同基因在不同组织中的表达水平。

TWAS利用基因表达数据建立基因表达量与遗传变异之间的预测模型。然后，这些模型用于估计特定遗传变异在特定样本中对基因表达的潜在影响。通过这种方式，TWAS可以推断出特定SNP的表达量变异如何与疾病风险相关。
> [!NOTE]+ 基因分型
>   基因分型是一种生物技术，用于分析个体的基因组中特定位置的遗传变异，包括单核苷酸多态性（SNP）、插入/缺失（indels）、拷贝数变异（CNVs）等。基因分型对于研究遗传病、人类遗传学、种群遗传学以及个体化医疗等领域都是非常重要的。

"FUSION" 是一套用于进行转录组广泛关联研究（TWAS）和调控组广泛关联研究（RWAS）的工具。这套工具的主要功能是<font color="#ffc000">构建预测模型</font>，这些模型能够预测基因对某一功能性或分子表型的影响（使用参考样本（例如来自GTEx项目的数据）中的基因表达数据或其他功能数据（如蛋白质水平、表观遗传标记等）与遗传变异数据（如SNP）来构建模型，预测特定遗传变异如何影响这些功能表型）。构建好的模型利用<font color="#ffc000">GWAS的总结摘要数据</font>来预测并检测这种基因组分的与疾病的关联（发现基因型与大规模疾病研究中发现的表型之间的统计关联）。最终目的是揭示那些可能只在少数样本中被测量到，但与更广泛人群中的疾病风险相关的功能表型。通过这种方式，FUSION帮助研究人员利用已有的遗传和表型数据，以及全基因组关联研究的成果，深入挖掘基因功能如何影响疾病的风险，进一步促进个性化医疗和精准医疗的发展。

```shell
ca ldsc
sum_file=174.11_PheCode.v1.0.fastGWA.gz
trait=t2

#/home/fenglijun/data3/soft/gwas/ldsc/munge_sumstats.py --sumstats ../$sum_file  --out $trait
#/home/fenglijun/data3/soft/gwas/ldsc/munge_sumstats.py --sumstats 26192919-GCST003044-EFO_0000384.h.tsv.gz --snp variant_id  --N 10000  --out t1
#/home/fenglijun/data3/soft/gwas/ldsc/munge_sumstats.py --sumstats cd_build37_40266_20161107.txt.gz  --out t2 --N 10000
#/home/fenglijun/data3/soft/gwas/ldsc/munge_sumstats.py --sumstats 250.1_PheCode.v1.0.fastGWA.gz  --out scz

ca ldsc
zcat ../$sum_file  | head -n 1 >q
for i in $(seq 1 22)
do
  zcat ../$sum_file | awk -v i="$i" 'BEGIN{OFS=FS="\t"} $1==i {print $0}'  > ${i}_snp.txt
  cat q ${i}_snp.txt > ${i}_snp_fi.txt
  /home/fenglijun/data3/soft/gwas/ldsc/munge_sumstats.py --sumstats ${i}_snp_fi.txt --out ${trait}_${i} &
done

ca deepmap2

for i in $(seq 1 22)
do
Rscript /home/fenglijun/data3/soft/gwas/fusion_twas-master/FUSION.assoc_test.R --sumstats ${trait}_${i}.sumstats.gz \
 --weights /home/fenglijun/data3/soft/gwas/fusion_twas-master/examples/weight2/sCCA_weights_v8/sCCA1.pos  --weights_dir /home/fenglijun/data3/soft/gwas/fusion_twas-master/examples/weight2/sCCA_weights_v8/   \
--ref_ld_chr /home/fenglijun/data3/soft/gwas/fusion_twas-master/LDREF/1000G.EUR.  --chr $i  \
--out  ${trait}_${i}.dat 
done 

# 筛选基因
cat *.dat >q
cat  q|cut -f 3,20|awk 'BEGIN{OFS=FS="\t"}{if($2<0.9) print $0}'|sed 's/\..*\t/\t/g' >q2
awk 'BEGIN{OFS=FS="\t"}ARGIND==1{a[$1]=$2}ARGIND==2{print $0,a[$4]}' q2 ~/data3/deep/data/human/reference/genome/g37/uniq_gene_15436.bed| awk 'BEGIN{OFS=FS="\t"} {if ($5 != None && $5 != "NA") print $0}'|grep -v 'NA' >q3

rm q4
while FS=OFS='\t' read -r col1 col2 col3; do
  a=locus_${col1}_${col2}_${col3}
  echo "${col1},${col2},${col3}" > $a
  sed -i 's/,/\t/g'  $a
  bedtools intersect -a $a -b q3 -wa -wb >2${a}
  l  2${a}|sort -k8g |head -n1|cut -f7 >>q4
done < ../4_coloc/locus.txt
sort q4|uniq >twas_prior_gene.txt

```


先将 `/home/data6/zhanghao/pathway/pathway/pigGWAS/twas/Exon_TWAS_SMultiXcan.csv.gz` 文件按照性状名进行划分文件，使用 `/home/data6/zhanghao/pathway/pathway/pigGWAS/data/pig/conbination.R` 进行

```shell
trait=

awk 'BEGIN{FS=OFS="\t"} {t=$8; $8=$3; $3=t; print $0}' /home/data6/zhanghao/pathway/pathway/pigGWAS/data/pig/conbination_twas/${trait}.txt | sort -k3,3 |awk 'BEGIN{FS=OFS="\t"} {sub(/:.*/, "", $2); print}' |cut -f1,2,3,4 |uniq > /home/data6/zhanghao/pathway/pathway/pigGWAS/benchmark/golden_standard/twas/${trait}.txt 

awk 'BEGIN{OFS=FS="\t"} ARGIND==1{a[$2]=$4} ARGIND==2{if($4 in a) print $0, a[$4]; else print $0, "NO MATCH"}' /home/data6/zhanghao/pathway/pathway/pigGWAS/benchmark/golden_standard/twas/${trait}.txt /home/data6/zhanghao/pathway/pathway/pigGWAS/lead_snp/re/sus11_protein_gene.bed |awk 'BEGIN{OFS=FS="\t"} {if ($5 != "NO MATCH" && $5 != "NA") print $0}'|grep -v 'NA' >/home/data6/zhanghao/pathway/pathway/pigGWAS/benchmark/golden_standard/twas/q3

rm q4 
while FS=OFS='\t' read -r col1 col2 col3; do 
	a=locus_${col1}_${col2}_${col3}
	echo "${col1},${col2},${col3}" > $a 
	sed -i 's/,/\t/g' $a 
	bedtools intersect -a $a -b q3 -wa -wb >2${a} 
	cat 2${a}|sort -k8g |head -n1|cut -f7 >>q4 
done < /home/data6/zhanghao/pathway/pathway/pigGWAS/benchmark/golden_standard/smr/locus_${trait}.txt 
sort q4|uniq > /home/data6/zhanghao/pathway/pathway/pigGWAS/benchmark/golden_standard/twas/q/${trait}_twas_prior_gene.txt 
```