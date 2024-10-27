---
{"dg-publish":true,"permalink":"/02 本科毕设/benchmark/4 coloc/","tags":["gardenEntry"]}
---


共定位方法汇总 [https://cloud.tencent.com/developer/article/2375207](https://cloud.tencent.com/developer/article/2375207)
e [caviar](https://github.com/fhormoz/caviar)
[https://github.com/fhormoz/caviar](https://github.com/fhormoz/caviar)
根据对来自英国生物库（UK Biobank）的95种复杂性状和GTEx v.8（第35参考文献）中的49种组织的eQTL（表达数量遗传性状）进行的精细定位结果，我们计算了类似于eCAVIAR软件报告的CLPP（共定位概率估计）。对于每个被精细定位的变异i，针对复杂性状g和eQTL性状e，CLPP被计算为P(Cig,Cie) = P(Cig) * P(Cie)，其中P(Cig)是变异i在复杂性状g中的PIP（后验概率）值，P(Cie)是变异i在eQTL性状e中的PIP值。这个量是变异对于复杂性状和基因表达性状都是原因的概率的估计。在每个精细定位的信号区域（CS）内和每个基因中，我们选取所有变异和GTEx组织中的最大CLPP值。
简而言之，CLPP 是一种用于估计变异在复杂性状和基因表达性状中是否都具有因果作用的概率的量。它基于精细定位结果和每个变异的后验概率值。通过计算 CLPP，我们可以评估变异对于复杂性状和基因表达性状的共定位概率，并确定是否存在潜在的因果关系。
coloc包
[https://zhuanlan.zhihu.com/p/392589375](https://zhuanlan.zhihu.com/p/392589375)
[https://blog.csdn.net/weixin_46587777/article/details/132368982](https://blog.csdn.net/weixin_46587777/article/details/132368982)


![02 本科毕设/benchmark/attachment/Untitled 1.png](/img/user/02%20%E6%9C%AC%E7%A7%91%E6%AF%95%E8%AE%BE/benchmark/attachment/Untitled%201.png)

左侧的图代表SNP在GWAS和QTL中-log10(p)的分布情况，p值越小越在Y轴上方 右侧两个各自分表代表QTL和GWAS自己的分布情况（横坐标是snp的pos位点位置） 纵坐标代表着SNP在该GWAS/QTL数据中的-log10(p)值，越高代表p值越小，lead SNP是在最顶。 R2为该数据在相应人群中某SNP与lead SNP之间的连锁程度。 主要展示数据中SNP的连锁不平衡情况 标出来的SNP，则是两个数据中PPH4最大的数值，具体每个SNP都有对应H1~H4的数据，查看结果中的共定位数据


先做血液的
```shell
trait=AG
sum_file=AG.fastGWA

ll ../finemapping/${trait}/${trait}_locus_snp_*|sed 's/_/\t/g'|cut -f 5,6,7 >locus.txt
qtl_file=/home/fenglijun/data3/project/data/human/gtex/hg19/v7/eQTL/GTEx_V7_cis_eqtl_summary/Whole_Blood.txt 
# 每个位点的gwas数据提取 
zcat ../${sum_file} |head -n 1 >q 
while IFS=
![02 本科毕设/benchmark/attachment/微信图片_20240320154953.png](/img/user/02%20%E6%9C%AC%E7%A7%91%E6%AF%95%E8%AE%BE/benchmark/attachment/%E5%BE%AE%E4%BF%A1%E5%9B%BE%E7%89%87_20240320154953.png)


```shell
sed 's/ \+/\t/g' Blood.D_ADG.enloc.snp.out > newfile.txt

sort -t



```shell

#将分隔符转为\t
#sed -e 's/ \+/\t/g' /home/data6/zhanghao/pathway/pathway/pigGWAS/data/pig/conbination_enloc/${trait}.txt |awk 'BEGIN{FS=OFS="\t"} {sub(/:.*/, "", $3); print}'>${trait}.txt

less -SN  /home/data6/zhanghao/pathway/pathway/pigGWAS/data/pig/conbination_enloc/${trait}.txt |awk 'BEGIN{FS=OFS="\t"} {sub(/:.*/, "", $3); print}'> /home/data6/zhanghao/pathway/pathway/pigGWAS/benchmark/golden_standard/encol/${trait}.txt
rm ./q/${trait}_enloc_prior_gene.txt 

while FS=OFS='\t' read -r col1 col2 col3; do 
a=locus_${col1}_${col2}_${col3} 
echo "${col1},${col2},${col3}" > ./q/$a 
sed -i 's/,/\t/g' ./q/$a 
bedtools intersect -a ./q/$a -b /home/data6/zhanghao/pathway/pathway/pigGWAS/lead_snp/re/sus11_protein_gene.bed -wa -wb|cut -f7 >./q/${a}_gene.txt  #按照位置排序后的全部蛋白编码的基因位置 sus11_protein_gene.bed
cat ${trait}.txt |grep --file ./q/${a}_gene.txt -w|cut -f3,8|sort -k2g|cut -f1|tail -n1|cut -f1 >> ./q/${trait}_enloc_prior_gene.txt 
done < locus_${trait}.txt 

sort ./q/${trait}_enloc_prior_gene.txt |uniq >./q/fi_${trait}_enloc_prior_gene.txt 

python /home/data6/zhanghao/pathway/pathway/presion_recall.py /home/data6/zhanghao/pathway/pathway/pigGWAS/benchmark/golden_standard/encol/all_golden_gene_in_locus_${trait}.txt /home/data6/zhanghao/pathway/pathway/pigGWAS/benchmark/golden_standard/encol/positive_gene_${trait}.list /home/data6/zhanghao/pathway/pathway/pigGWAS/benchmark/golden_standard/encol/q/fi_${trait}_enloc_prior_gene.txt ${trait} >>q1 2>w


rm closest_gene.txt
rm loc_gene_num.txt

mkdir -p 2_cloest
cd 2_cloest
while FS=OFS='\t' read -r col1 col2 col3; do
  echo "process locus_${col1}_${col2}_${col3}"
  echo -e "${col1}\t${col2}\t${col3}\ta" > sig_loc_bed.txt
	a=${col1}_${col2}_${col3}
	 l /home/data6/zhanghao/pathway/pathway/pigGWAS/benchmark/golden_standard/smr/fine_mapping_all_snp_${trait}.bed_${a} |head -n 1 |awk  'BEGIN{OFS=FS="\t"}{print $1,$3,$3,$2}' >snp.bed
	 bedtools closest -a snp.bed -b  /home/data6/zhanghao/pathway/pathway/pigGWAS/lead_snp/re/sus11_protein_gene.bed |cut -f8 >> closest_gene.txt
	 bedtools intersect -a sig_loc_bed.txt -b /home/data6/zhanghao/pathway/pathway/pigGWAS/lead_snp/re/sus11_protein_gene.bed | wc |sed 's/  */\t/g'|cut -f2 >>loc_gene_num.txt
done < /home/data6/zhanghao/pathway/pathway/pigGWAS/benchmark/golden_standard/smr/locus_${trait}.txt

# 打印出p_value enrichment se
rm q

python  /home/data6/zhanghao/pathway/pathway/cloest_gene_enrichment.py closest_gene.txt loc_gene_num.txt  /home/data6/zhanghao/pathway/pathway/pigGWAS/benchmark/golden_standard/encol/q/fi_${trait}_enloc_prior_gene.txt  coloc >>q 2>w



```
{ #041d5d}

\t' read -r col1 col2 col3; do 
a="locus_${col1}_${col2}_${col3}" 
zcat ../${sum_file} |awk -v col1="$col1" -v col2="$col2" -v col3="$col3" ' 
BEGIN { OFS = FS = "\t" } $1 == col1 && $3 >= col2 && $3 <= col3 { print $0 } 
' > "${a}.txt" & 
done < locus.txt 

# 每个位点的qtl数据提取 
cat ${qtl_file} |head -n 1 >q_qtl 
while IFS=
![02 本科毕设/benchmark/attachment/微信图片_20240320154953.png](/img/user/02%20%E6%9C%AC%E7%A7%91%E6%AF%95%E8%AE%BE/benchmark/attachment/%E5%BE%AE%E4%BF%A1%E5%9B%BE%E7%89%87_20240320154953.png)


{{CODE_BLOCK_1}}



{{CODE_BLOCK_2}}
{ #041d5d}

\t' read -r col1 col2 col3; do 
a="locus_${col1}_${col2}_${col3}" 
cat "${qtl_file}" | awk -v col1="$col1" -v col2="$col2" -v col3="$col3" ' 
BEGIN { OFS = FS = "\t" } 
$2 == col1 && $3 >= col2 && $3 <= col3 { print $0 } 
' > "q_${a}.txt" & 
done < locus.txt 

while IFS=
![02 本科毕设/benchmark/attachment/微信图片_20240320154953.png](/img/user/02%20%E6%9C%AC%E7%A7%91%E6%AF%95%E8%AE%BE/benchmark/attachment/%E5%BE%AE%E4%BF%A1%E5%9B%BE%E7%89%87_20240320154953.png)


{{CODE_BLOCK_1}}



{{CODE_BLOCK_2}}
{ #041d5d}

\t' read -r col1 col2 col3; do 
a="locus_${col1}_${col2}_${col3}" 
cat q_qtl q_${a}.txt >"qtl_${a}.txt" 
rm q_${a}.txt 
done < locus.txt 

while IFS=
![02 本科毕设/benchmark/attachment/微信图片_20240320154953.png](/img/user/02%20%E6%9C%AC%E7%A7%91%E6%AF%95%E8%AE%BE/benchmark/attachment/%E5%BE%AE%E4%BF%A1%E5%9B%BE%E7%89%87_20240320154953.png)


{{CODE_BLOCK_1}}



{{CODE_BLOCK_2}}
{ #041d5d}

\t' read -r col1 col2 col3; do 
a="locus_${col1}_${col2}_${col3}" 
cat "q" "${a}.txt" > "gwas_${a}.txt" 
rm "${a}.txt" 
done < locus.txt 

## 运行共定位R脚本 
while IFS=
![02 本科毕设/benchmark/attachment/微信图片_20240320154953.png](/img/user/02%20%E6%9C%AC%E7%A7%91%E6%AF%95%E8%AE%BE/benchmark/attachment/%E5%BE%AE%E4%BF%A1%E5%9B%BE%E7%89%87_20240320154953.png)


{{CODE_BLOCK_1}}



{{CODE_BLOCK_2}}
{ #041d5d}

\t' read -r col1 col2 col3; do 
a="locus_${col1}_${col2}_${col3}" 
/home/fenglijun/data3/soft/miniconda3/envs/deepmap2/bin/Rscript /home/fenglijun/data3/script/post_gwas/coloc_fi.R \ 
${a} \ 
gwas_${a}.txt \ 
qtl_${a}.txt \ 
./ & done < locus.txt 

## 筛选每个位点的共定位基因 
rm coloc_prior_gene.txt 
while IFS=
![02 本科毕设/benchmark/attachment/微信图片_20240320154953.png](/img/user/02%20%E6%9C%AC%E7%A7%91%E6%AF%95%E8%AE%BE/benchmark/attachment/%E5%BE%AE%E4%BF%A1%E5%9B%BE%E7%89%87_20240320154953.png)


{{CODE_BLOCK_1}}



{{CODE_BLOCK_2}}
{ #041d5d}

\t' read -r col1 col2 col3; do 
a="coloc_locus_${col1}_${col2}_${col3}.txt" 
cat ${a}|head -n 2|tail -n1 |cut -f1 |sed 's/\..*//g' >>coloc_prior_gene.txt 
done < locus.txt 

sort coloc_prior_gene.txt|uniq >q 
mv q coloc_prior_gene.txt
```
![02 本科毕设/benchmark/attachment/微信图片_20240320154953.png](/img/user/02%20%E6%9C%AC%E7%A7%91%E6%AF%95%E8%AE%BE/benchmark/attachment/%E5%BE%AE%E4%BF%A1%E5%9B%BE%E7%89%87_20240320154953.png)


{{CODE_BLOCK_1}}



{{CODE_BLOCK_2}}
{ #041d5d}

\t' -k6,6g newfile.txt > result.txt
```



{{CODE_BLOCK_2}}
{ #041d5d}

\t' read -r col1 col2 col3; do 
a="locus_${col1}_${col2}_${col3}" 
zcat ../${sum_file} |awk -v col1="$col1" -v col2="$col2" -v col3="$col3" ' 
BEGIN { OFS = FS = "\t" } $1 == col1 && $3 >= col2 && $3 <= col3 { print $0 } 
' > "${a}.txt" & 
done < locus.txt 

# 每个位点的qtl数据提取 
cat ${qtl_file} |head -n 1 >q_qtl 
while IFS=
![02 本科毕设/benchmark/attachment/微信图片_20240320154953.png](/img/user/02%20%E6%9C%AC%E7%A7%91%E6%AF%95%E8%AE%BE/benchmark/attachment/%E5%BE%AE%E4%BF%A1%E5%9B%BE%E7%89%87_20240320154953.png)


{{CODE_BLOCK_1}}



{{CODE_BLOCK_2}}
{ #041d5d}

\t' read -r col1 col2 col3; do 
a="locus_${col1}_${col2}_${col3}" 
cat "${qtl_file}" | awk -v col1="$col1" -v col2="$col2" -v col3="$col3" ' 
BEGIN { OFS = FS = "\t" } 
$2 == col1 && $3 >= col2 && $3 <= col3 { print $0 } 
' > "q_${a}.txt" & 
done < locus.txt 

while IFS=
![02 本科毕设/benchmark/attachment/微信图片_20240320154953.png](/img/user/02%20%E6%9C%AC%E7%A7%91%E6%AF%95%E8%AE%BE/benchmark/attachment/%E5%BE%AE%E4%BF%A1%E5%9B%BE%E7%89%87_20240320154953.png)


{{CODE_BLOCK_1}}



{{CODE_BLOCK_2}}
{ #041d5d}

\t' read -r col1 col2 col3; do 
a="locus_${col1}_${col2}_${col3}" 
cat q_qtl q_${a}.txt >"qtl_${a}.txt" 
rm q_${a}.txt 
done < locus.txt 

while IFS=
![02 本科毕设/benchmark/attachment/微信图片_20240320154953.png](/img/user/02%20%E6%9C%AC%E7%A7%91%E6%AF%95%E8%AE%BE/benchmark/attachment/%E5%BE%AE%E4%BF%A1%E5%9B%BE%E7%89%87_20240320154953.png)


{{CODE_BLOCK_1}}



{{CODE_BLOCK_2}}
{ #041d5d}

\t' read -r col1 col2 col3; do 
a="locus_${col1}_${col2}_${col3}" 
cat "q" "${a}.txt" > "gwas_${a}.txt" 
rm "${a}.txt" 
done < locus.txt 

## 运行共定位R脚本 
while IFS=
![02 本科毕设/benchmark/attachment/微信图片_20240320154953.png](/img/user/02%20%E6%9C%AC%E7%A7%91%E6%AF%95%E8%AE%BE/benchmark/attachment/%E5%BE%AE%E4%BF%A1%E5%9B%BE%E7%89%87_20240320154953.png)


{{CODE_BLOCK_1}}



{{CODE_BLOCK_2}}
{ #041d5d}

\t' read -r col1 col2 col3; do 
a="locus_${col1}_${col2}_${col3}" 
/home/fenglijun/data3/soft/miniconda3/envs/deepmap2/bin/Rscript /home/fenglijun/data3/script/post_gwas/coloc_fi.R \ 
${a} \ 
gwas_${a}.txt \ 
qtl_${a}.txt \ 
./ & done < locus.txt 

## 筛选每个位点的共定位基因 
rm coloc_prior_gene.txt 
while IFS=
![02 本科毕设/benchmark/attachment/微信图片_20240320154953.png](/img/user/02%20%E6%9C%AC%E7%A7%91%E6%AF%95%E8%AE%BE/benchmark/attachment/%E5%BE%AE%E4%BF%A1%E5%9B%BE%E7%89%87_20240320154953.png)


{{CODE_BLOCK_1}}



{{CODE_BLOCK_2}}
{ #041d5d}

\t' read -r col1 col2 col3; do 
a="coloc_locus_${col1}_${col2}_${col3}.txt" 
cat ${a}|head -n 2|tail -n1 |cut -f1 |sed 's/\..*//g' >>coloc_prior_gene.txt 
done < locus.txt 

sort coloc_prior_gene.txt|uniq >q 
mv q coloc_prior_gene.txt
```
![02 本科毕设/benchmark/attachment/微信图片_20240320154953.png](/img/user/02%20%E6%9C%AC%E7%A7%91%E6%AF%95%E8%AE%BE/benchmark/attachment/%E5%BE%AE%E4%BF%A1%E5%9B%BE%E7%89%87_20240320154953.png)


{{CODE_BLOCK_1}}



{{CODE_BLOCK_2}}
{ #041d5d}

