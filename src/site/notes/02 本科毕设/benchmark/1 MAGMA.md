---
{"dg-publish":true,"permalink":"/02 本科毕设/benchmark/1 MAGMA/","tags":["gardenEntry"]}
---


[MAGMA 基因及通路分析 – GWASLab – GWAS实验室](https://gwaslab.org/2021/12/13/magma/)
利用MAGMA计算每个基因的基因水平关联统计量及基因间相关性。
```shell
# SNP注释 Annotation
# snp-loc 文件应当包含三列，SNPID，CHR以及POS，这里也可以直接使用plink的bim文件
# gene-loc 第一列为ensg编号 第二列染色体编号 第三列第四列染色体位置 第五列正负链 第六列gene name
# 注释完成后会得到 .genes.annot文件,内容为第一列为geneID 之后为在这个基因内的SNP
sum_file=AG.fastGWA
trait=AG
n_sample=517

cat ../../$sum_file |awk 'BEGIN{OFS=FS="\t"}{print $2,$1,$3}' > ${trait}_snp.loc
cat ../../$sum_file |awk 'BEGIN{OFS=FS="\t"}{print $2,$10}' > ${trait}_p_value.loc
magma \
    --annotate \
    --snp-loc ${trait}_snp.loc \
    --gene-loc /home/data6/zhanghao/pathway/pathway/pigGWAS/panel/reference/uniq_gene.loc \
    --out $trait


# 1 基因关联分析：将GWAS摘要统计数据映射到基因上。这个步骤涉及将SNP的关联信号聚合到基因水平，以评估每个基因与表型的关联程度。
# bfile为原始数据或是参考LD面板，如果数据量不大可以直接使用自己的plink的bed格式原始数据，在原始数据无法获得的时候可以使用magma提供的1000 genome参考数据，b
# iobank级别的数据的情况下，可以随机抽取某个族裔无亲缘关系的一定人数（例如20000人）来构建自己的参考面板。

# pval为SNP的p值文件，包含两列 SNP 以及 P， 这里N为样本量，
# gene-annot为上一步注释后得到的文件。
# 完成后有两个文件输出；1 .genes.out 2 ..genes.raw
n_sample=517
magma \
    --bfile /home/data6/zhanghao/pathway/pathway/pigGWAS/panel/pig --gene-annot ${trait}.genes.annot --pval ${trait}_p_value.loc N=${n_sample} --out $trait
# .genes.out 基因关联分析的结果
#.genes.raw 在第三步通路分析时会使用到这个文件
                 
#2 基因集分析：评估特定基因组（如生物通路或功能群体）是否比随机期望更显著地与表型相关联。
#gmt文件第一列为通路名，第二列为链接，但这里也用通路名与第一列一致，第三列及后面为基因symbol以\t分隔（与uniq_gene.loc中一致） 

magma \
    --gene-results ${trait}.genes.raw \
    --set-annot /home/data6/zhanghao/pathway/pathway/pathway/pig.gmt  \
    --out ${trait}

#3 基因属性分析（可选）：分析基因属性（如特定生物标记物的表达水平）与基因关联强度之间的关系。
#magma --gene-results [基因关联分析结果路径] --gene-covar [基因属性文件路径] --out [输出文件前缀]


# 4 筛选基因
sed 's/  */\t/g'  ${trait}.genes.out |awk  'BEGIN{OFS=FS="\t"}{if($9<0.05) print $0}' |cut -f1 |sort | uniq >magma_prior_gene.txt
wc  magma_prior_gene.txt

less ${trait}.gsa.out |sed 's/  */\t/g'|tail -n+5|awk 'BEGIN{OFS=FS="\t"}{if($7<0.05) print $0}'


less ${trait}.gsa.out |sed 's/  */\t/g'|tail -n+5|awk 'BEGIN{OFS=FS="\t"}{if($7<0.05) print $0}'|sort -k7g|head -n 50 |cut -f1 >q
less /home/data6/zhanghao/pathway/pathway/pathway/pig.gmt |grep --file q -w|cut -f3-|tr '\t'  '\n'|sort |uniq > geneset_prior_gene.txt
wc geneset_prior_gene.txt

```


##  final
```shell
#使用GWAS结果/home/data6/zhanghao/pathway/pathway/pigGWAS/magma/PigBiobank_gene-based.csv
#更改分隔符
sed -i 's/,/\t/g' PigBiobank_gene-based.csv PigBiobank_gene-based.csv 

less PigBiobank_gene-based.csv |cut -f1,2,4,5,6,11 |grep "${trait}" |sed '1d' >${trait}.txt

rm ./q/${trait}_magma_prior_gene.txt 

while FS=OFS='\t' read -r col1 col2 col3; do 
a=locus_${col1}_${col2}_${col3} 
echo "${col1},${col2},${col3}" > ./q/$a 
sed -i 's/,/\t/g' ./q/$a 
bedtools intersect -a ./q/$a -b /home/data6/zhanghao/pathway/pathway/pigGWAS/lead_snp/re/sus11_protein_gene.bed -wa -wb|cut -f7 >./q/${a}_gene.txt 
cat ${trait}.txt |grep --file ./q/${a}_gene.txt -w|cut -f2,6|sort -k2g|cut -f1|head -n1 >> ./q/${trait}_magma_prior_gene.txt 
done < locus_${trait}.txt 

sort ./q/${trait}_magma_prior_gene.txt |uniq >./q/fi_${trait}_magma_prior_gene.txt 

python /home/data6/zhanghao/pathway/pathway/presion_recall.py /home/data6/zhanghao/pathway/pathway/pigGWAS/benchmark/golden_standard/magma/all_golden_gene_in_locus_${trait}.txt /home/data6/zhanghao/pathway/pathway/pigGWAS/benchmark/golden_standard/magma/positive_gene_${trait}.list /home/data6/zhanghao/pathway/pathway/pigGWAS/benchmark/golden_standard/magma/q/fi_${trait}_magma_prior_gene.txt ${trait} >>q1 2>w



```