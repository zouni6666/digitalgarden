---
{"dg-publish":true,"permalink":"/02 本科毕设/benchmark/3 SMR/","tags":["gardenEntry"]}
---


软件官网链接 [https://yanglab.westlake.edu.cn/software/smr/#DataManagement](https://yanglab.westlake.edu.cn/software/smr/#DataManagement)
[https://zhuanlan.zhihu.com/p/337969763](https://zhuanlan.zhihu.com/p/337969763)
看这个 smr+coloc [https://zhuanlan.zhihu.com/p/659654054](https://zhuanlan.zhihu.com/p/659654054)
我们使用SMR软件工具，并使用来自GTEx v.7（参考文献33）的48个人类组织的预计算cis-eQTL摘要数据，应用了SMR 34。对于每个基因，cis-eQTL摘要数据已经预先筛选出位于距离转录起始位点1Mb范围内的SNP。对于精确度-召回率分析，SMR运行时排除了编码信号。
细心的朋友会发现，刚才我们并没有下载 eqtl 的 txt 列表格式的文件，理论上和 gwas 一样，eqtl 也应该有一个非常大数据量的列表文件，那么如何得到呢，依然是靠 smr 软件

先做血液的
```shell
# 获取txt的数据eqtl文件（不用运行）
smr_linux_x86_64  --beqtl-summary Whole_Blood --query 0.05 --out Whole_Blood

ca gwas(conda activate R)
sum_file=AG.fastGWA
trait=t10

cat ../$sum_file |awk 'BEGIN{OFS=FS="\t"}{print $2,$4,$5,$7,$11,$12,$13,$6}' >smr_${sum_file}.input
# 软件运行
#mygwas.ma为GWAS的summary文件，包含SNP、A1、A2、freq、b、se、p、N，N指样本数（N没有的话，可以直接赋予NA），内容如下所示
#Note that the column "n" will not be used in either SMR or HEIDI analysis and thus can be replaced by "NA" if not available
# For a case-control study, the effect size should be log(odds ratio) with its corresponding standard error. 
nohup smr_linux_x86_64 \
--bfile /home/data6/zhanghao/pathway/pathway/pigGWAS/panel/pig \
--gwas-summary smr_$sum_file.input \
--beqtl-summary /home/data6/zhanghao/pathway/pathway/pigGWAS/eqtl/PigGTEx_v0.finemapped_eQTL/Blood.susieinf.gz \
--out ${trait}_smr \
--thread-num 60 &

## 筛选基因
zcat ./SMR/Blood.pcg.smr.gz | awk 'BEGIN{OFS=FS="\t"} {print $3,$4,$6,$21}' | sort -k4,4g | awk 'BEGIN{OFS=FS="\t"} {print $2,$3,$3,$1,$4}' > smr_prior_gene_min.txt


while FS=OFS='\t' read -r col1 col2 col3; do 
a=locus_${col1}_${col2}_${col3} 
echo "${col1},${col2},${col3}" > $a 
sed -i 's/,/\t/g' $a 
cat /home/data6/zhanghao/pathway/pathway/pigGWAS/smr/smr_prior_gene_min.txt |bedtools intersect -a $a -b stdin -wa -wb >2${a} 
less 2${a}|sort -k8g |head -n1|cut -f7 >>q4 
done < /home/data6/zhanghao/pathway/pathway/pigGWAS/finemapping/AG_FINAL/locus.txt
sort q4|uniq >twas_prior_gene.txt

```
##### 
######
使用 [PigGTEx-Portal](https://piggtex.farmgtex.org/) 中下载的 SMR 结果
先看血液的  `/home/data6/zhanghao/pathway/pathway/pigGWAS/smr/SMR/Blood. pcg. smr. gz`

```shell
zcat Blood.pcg.smr.gz |cut -f18 |head

zcat ./SMR/Blood.pcg.smr.gz | awk 'BEGIN{OFS=FS="\t"; min=) {min=$19; gene=$3} print $0} END{print gene}' | tail -n 1 > smr_prior_gene_min.txt

```

```shell
zcat ./SMR/Blood.pcg.smr.gz | awk 'BEGIN{OFS=FS="\t"} {print $3,$4,$6,$21}' | sort -k4,4g | awk 'BEGIN{OFS=FS="\t"} {print $2,$3,$3,$1,$4}' > smr_prior_gene_min.txt

bedtools intersect -a smr_prior_gene_min2.txt -b /home/data6/zhanghao/pathway/pathway/pigGWAS/finemapping/AG_FINAL/locus.txt -wa -wb > results.txt

```
结果：/home/data 6/zhanghao/pathway/pathway/pigGWAS/smr/results. txt



#### zuixin
/home/data 6/zhanghao/pathway/pathway/pigGWAS/data/pig/conbination_smr   
```shell
less -SN M_TNB_smr.txt | awk 'BEGIN{OFS=FS="\t"} {print $3,$4,$6,$19}' |sort -k4,4g | awk 'BEGIN{OFS=FS="\t"} {print $2,$3,$3,$1,$4}' >/home/data6/zhanghao/pathway/pathway/pigGWAS/benchmark/golden_standard/smr/smr_prior_gene_min_M_TNB.txt


#每个位点下选择一个最显著p值最小的基因 有多少位点多少基因 有相同的合并 finemap
trait=
cat /home/data6/zhanghao/pathway/pathway/pigGWAS/lead_snp/PigBiobank_lead_SNP_in_300_studies.csv |awk 'BEGIN{OFS=FS=","}{print $4,$5,$5,$1}' |sed 's/,/\t/g' |grep "${trait}" |sed '1d' >fine_mapping_all_snp_${trait}.bed
awk 'BEGIN{OFS=FS="\t"}{print $1,$2-500000,$2+500000}' fine_mapping_all_snp_${trait}.bed | awk 'BEGIN{OFS=FS="\t"}{if($2<0) print $1,0,$3; else print $0}'> /home/data6/zhanghao/pathway/pathway/pigGWAS/benchmark/golden_standard/smr/locus_${trait}.txt



#bedtools intersect -a /home/data6/zhanghao/pathway/pathway/pigGWAS/benchmark/golden_standard/smr/smr_prior_gene_min_M_TNB.txt -b /home/data6/zhanghao/pathway/pathway/pigGWAS/benchmark/golden_standard/smr/locus_M_TNUM.txt -wa -wb >/home/data6/zhanghao/pathway/pathway/pigGWAS/benchmark/golden_standard/smr/results_M_TNUM.txt





```

```shell
trait=
cat /home/data6/zhanghao/pathway/pathway/pigGWAS/lead_snp/PigBiobank_lead_SNP_in_300_studies.csv |awk 'BEGIN{OFS=FS=","}{print $4,$5,$5,$1}' |sed 's/,/\t/g' |grep "${trait}" |sed '1d' >fine_mapping_all_snp_${trait}.bed   #finemap结果
awk 'BEGIN{OFS=FS="\t"}{print $1,$2-500000,$2+500000}' fine_mapping_all_snp_${trait}.bed | awk 'BEGIN{OFS=FS="\t"}{if($2<0) print $1,0,$3; else print $0}'> /home/data6/zhanghao/pathway/pathway/pigGWAS/benchmark/golden_standard/smr/locus_${trait}.txt

rm ./q/${trait}_smr_prior_gene.txt 

while FS=OFS='\t' read -r col1 col2 col3; do 
	a=locus_${col1}_${col2}_${col3} 
	echo "${col1},${col2},${col3}" > ./q/$a 
	sed -i 's/,/\t/g' ./q/$a 
	bedtools intersect -a ./q/$a -b /home/data6/zhanghao/pathway/pathway/pigGWAS/lead_snp/re/sus11_protein_gene.bed -wa -wb|cut -f7 >./q/${a}_gene.txt 
	cat /home/data6/zhanghao/pathway/pathway/pigGWAS/data/pig/conbination_smr/${trait}_smr.txt |grep --file ./q/${a}_gene.txt -w|cut -f 3,21|sort -k2g|cut -f1|head -n 1 >> ./q/${trait}_smr_prior_gene.txt 
done < locus_${trait}.txt 
cat ./q/${trait}_smr_prior_gene.txt |sort|uniq >./q/fi_${trait}_smr_prior_gene.txt

trait=
bedtools intersect -a /home/data6/zhanghao/pathway/pathway/pigGWAS/benchmark/golden_standard/smr/fine_mapping_all_snp_${trait}.bed -b /home/data6/zhanghao/pathway/pathway/pigGWAS/lead_snp/uniq_gene_exon2.txt -wa -wb |cut -f8|sort|uniq >/home/data6/zhanghao/pathway/pathway/pigGWAS/benchmark/golden_standard/smr/positive_gene_${trait}.list

bedtools intersect -a /home/data6/zhanghao/pathway/pathway/pigGWAS/benchmark/golden_standard/smr/fine_mapping_all_snp_${trait}.bed -b /home/data6/zhanghao/pathway/pathway/pigGWAS/lead_snp/uniq_gene_exon2.txt -wa -wb |cut -f 5,6,7,8|sort|uniq|awk 'BEGIN{OFS=FS="\t"}{print $1,$2-500000,$3+500000,$4}' >/home/data6/zhanghao/pathway/pathway/pigGWAS/benchmark/golden_standard/smr/golden_locus_${trait}.bed

bedtools intersect -a /home/data6/zhanghao/pathway/pathway/pigGWAS/benchmark/golden_standard/smr/golden_locus_${trait}.bed -b /home/data6/zhanghao/pathway/pathway/pigGWAS/lead_snp/uniq_gene_all2.bed -wa -wb |cut -f8  |sort|uniq >/home/data6/zhanghao/pathway/pathway/pigGWAS/benchmark/golden_standard/smr/all_golden_gene_in_locus_${trait}.txt

python /home/data6/zhanghao/pathway/pathway/presion_recall.py /home/data6/zhanghao/pathway/pathway/pigGWAS/benchmark/golden_standard/smr/all_golden_gene_in_locus_${trait}.txt /home/data6/zhanghao/pathway/pathway/pigGWAS/benchmark/golden_standard/smr/positive_gene_${trait}.list /home/data6/zhanghao/pathway/pathway/pigGWAS/benchmark/golden_standard/smr/q/fi_${trait}_smr_prior_gene.txt ${trait} >>q1 2>w
```



