


Given a 6-column TSV (chr, start, end, plus, minus, asym where asym=(plus−minus)/(plus+minus)),
write R code that:
- filters total>0
- identify regions with strong/moderate plus/minus strand bias & color in plot
- plots asymmetry vs genomic position with ggplot2, faceted by chromosome
e.g. input: NC_085930.1 53800000 53900000 8 0 -1
Return only R code.

If I have bed file, can I find local regions of high strand asymmetry within chromosome. Count Strand asymmetry (SA) is calculated by equation = (counts in +ve - counts in -ve strand)/(counts in +ve strand - counts in -ve strand).
SA value close to -1 means -ve strand bias & close to +1 means +ve strand bias. A higher asymmetry is observed in +ve strand of borangutan, but I want to to narrow down to the exact regions contributing to this asymmetry.
e.g. 
species  count  strand  count  strand  count_difference  asymmetry
borangutan  113486  +  51260  -  1  62226  0.377709

e.g.: input bed file to identify local regions of asymmetry
NC_085930.1	121463	121853	3	1	+
NC_072373.2	100110165	100110309	-1	1	-

#🔹 Step 1: Make fixed genomic windows
#bedtools makewindows -g /media/aswin/SCFR/SCFR-main/genome_sizes/borangutan.genome -w 100000 > borangutan.windows.bed
bedtools makewindows -g /media/aswin/SCFR/SCFR-main/genome_sizes/borangutan.genome -w 100000 > borangutan_windows.bed

#🔹 Step 2: Intersect BED with windows
awk 'NR>1{print$1,$2,$3,$4,$5,$6}' OFS="\t" ../borangutan_single_exon.tsv > borangutan_single_exon_scfr.bed
bedtools intersect -a borangutan_windows.bed -b borangutan_single_exon_scfr.bed -wa -wb > window_hits.bed

#🔹 Step 3: Count + and − per window
awk '{key = $1 FS $2 FS $3
  if ($NF == "+") plus[key]++
  else if ($NF == "-") minus[key]++}
END { for (k in plus) {
    p = plus[k]
    m = minus[k] + 0
    if (p + m > 0)
      print k, p, m, (p-m)/(p+m)}
  for (k in minus) if (!(k in plus)) {
    p = 0; m = minus[k]
    print k, p, m, (p-m)/(p+m)}}' window_hits.bed > window_SA.tsv

awk '
{
  key = $1 FS $2 FS $3
  if ($9 == "+") plus[key]++
  else if ($9 == "-") minus[key]++
}
END {
  for (k in plus) {
    p = plus[k] + 0
    m = minus[k] + 0
    if (p + m > 0)
      print k, p, m, (m - p) / (m + p)
  }
}
' window_hits.bed > window_asymmetry.tsv

#🔹 Step 4: Extract highly asymmetric windows
#awk '$6 > 0.5 {print}' window_asymmetry.tsv > neg_strand_hotspots.bed
awk '$6>=0.6 && ($4+$5)>=10' window_SA.tsv > pos_hotspots.bed
awk '$6<=-0.6 && ($4+$5)>=10' window_SA.tsv > neg_hotspots.bed


#######################

while read i
do
t=$(echo $i | awk '{print$1}')
c=$(echo $i | awk '{print$2,$3}')
cs=$(awk -v t="$t" '{if($0~t) print$4,$5}' /media/aswin/SCFR/SCFR-main/genes/human/GCF_009914755.1_T2T-CHM13v2.0_genomic.gtf)
m=$(echo "$cs" | grep "$c")
if [[ -z $m ]]
then
f="multi"
else
f="single"
fi
echo $t $c $f
unset t c cs m f
done < <(awk '$14<2 {print$11,$8,$9}' human_multi_exon.tsv) > check_single_exons_in_multi_scfrs

##########################

bedtools intersect -a scfr.bed -b cds.bed -wo -s > scfr_cds_overlaps.bed

while read scfr 
do
f=$(echo $scfr | awk '{print$5}'
grep "$scfr" scfr_cds_overlaps.bed | awk '$5==$11 && $2<=$8 && $3>=$9 {print$0,$8-$2,$3-$9}'
done < scfr.bed


#Refomrat 
#1m58.476s
time awk '{if($4~"-") print$0,"1","-"; else print$0,"1","+"}' OFS="\t" /media/aswin/SCFR/SCFR-main/SCFR_all/human_SCFR_all.out > human_scfr_all.bed
awk '{if($6~"-") print$1,$2,$3,$4,$6,"-",$5,$7,$8,$9,$10,$11,$12; else print$1,$2,$3,$4,$6,"+",$5,$7,$8,$9,$10,$11,$12}' OFS="\t" ../human_coding_exons.bed > human_coding_exons_reformatted.bed
#Get overlaps ()
time bedtools intersect -a human_scfr_all.bed -b human_coding_exons_reformatted.bed -wo -s > human_scfr_cds_all_overlaps.bed

#Convert gtf to cds bed (headers: chrom start end gene tx strand frame exon_number gene_tx_count exon_tx_count exon_order sharing splicing)
time python3 gtf_to_cds_with_transcript_exon_metadata_2.py -i /media/aswin/SCFR/SCFR-main/genes/human/GCF_009914755.1_T2T-CHM13v2.0_genomic.gtf -o human_coding_exons.bed
#Get overlaps (2m14.940s) (headers: chrom start end frame filler strand chrom start end gene tx strand frame exon_number gene_tx_count exon_tx_count exon_order sharing splicing)
time bedtools intersect -a human_scfr_all.bed -b human_coding_exons.bed -wo -s > human_scfr_cds_all_overlaps.bed

while read scfr 
do
f=$(echo $scfr | awk '{print$4}' | tr -d "-")
#restrict exon overlaps to within scfr or at scfr boundaries & have same frame
grep "$scfr" human_scfr_cds_all_overlaps.bed | awk '{gsub(/-/,"",$4)}1' OFS="\t" | awk '$4==$13 && $2<=$8 && $3>=$9 {print$0}' > temp.bed
#Total number of exons
tnoe=$(awk '{print$7,$8,$9,$10,$11,$12}' OFS="\t" temp.bed | sort -u | wc -l)
#Overlapping exons
mc=$(bedtools sort -i temp.bed | bedtools merge -i - | wc -l)
if [[ $tnoe == "1" ]] || [[ $mc == "1" ]]
then
cat temp.bed | awk '{if($6~"-") print$0,$3-$9,$8-$2; else print$0,$8-$2,$3-$9}' OFS="\t" | awk '!seen[$7,$8,$9,$10,$12,$13,$17,$19,$20,$21,$22]++' >> single_exon.tsv

c=$(bedtools sort -i temp.bed | wc -l)

grep "$scfr" human_scfr_cds_all_overlaps.bed | awk '$4==$11 && $2<=$8 && $3>=$9 {print$0,$8-$2,$3-$9}'
done < <(awk '{print$1,$2,$3,$4,$5,$6}' human_scfr_cds_all_overlaps.bed | sort -u | tr " " "\t")

scfr=`grep KLHL17 human_scfr_cds_all_overlaps.bed | grep  392806 | grep  393139 | grep XM_054336255.1 | grep last`
scfr=`grep BRCA1 human_scfr_cds_all_overlaps.bed | grep  43926117 | grep 43926252`
scfr=`grep TTC23L human_scfr_cds_all_overlaps.bed | grep  35087777 | grep 35088032 | awk '{print$1,$2,$3,$4,$5,$6}' OFS="\t" | sort -u`
scfr=$(grep DBI human_scfr_cds_all_overlaps.bed | grep 119801067 | awk '{print$1,$2,$3,$4,$5,$6}' OFS="\t" | sort -u)

scfr=`grep KLHL17 human_scfr_cds_all_overlaps.bed | grep 393139 | grep 392806 | grep XM_054336255.1 | grep 392897`

awk '!seen[$7,$8,$9,$10,$12,$13,$17,$19,$20,$21,$22]++' 

###############################

#get zero based scfr bed file (2m1.660s)
time awk '{if($4~"-") print$0,"1","-"; else print$0,"1","+"}' OFS="\t" /media/aswin/SCFR/SCFR-main/SCFR_all/human_SCFR_all.out > human_scfr_all.bed
#Convert gtf to zero based cds bed (0m19.108s) (headers: chrom start end gene tx strand frame exon_number gene_tx_count exon_tx_count exon_order sharing splicing)
time python3 gtf_to_cds_with_transcript_exon_metadata_2.py -i /media/aswin/SCFR/SCFR-main/genes/human/GCF_009914755.1_T2T-CHM13v2.0_genomic.gtf -o human_coding_exons.bed
#Get overlaps (2m14.940s) (headers: chrom start end frame filler strand chrom start end gene tx strand frame exon_number gene_tx_count exon_tx_count exon_order sharing splicing)
time bedtools intersect -a human_scfr_all.bed -b human_coding_exons.bed -wo -s > human_scfr_cds_all_overlaps.bed
#Get unique SCFRs containing xoplete exons within or at boundaries
awk '$2<=$8 && $3>=$9 {print$0}' human_scfr_cds_all_overlaps.bed > human_scfr_containing_cds.bed
awk '$2<=$8 && $3>=$9 {print$0}' human_scfr_cds_all_overlaps.bed | awk '{print$1,$2,$3,$4,$5,$6}' | sort -u | tr " " "\t" > human_scfr_containing_cds.bed

while read scfr 
do
#strand
s=$(echo $scfr | awk '{print$6}')

f=$(echo $scfr | awk '{print$4}' | tr -d "-")
#restrict exon overlaps to within scfr or at scfr boundaries & have same frame
grep "$scfr" human_scfr_containing_cds.bed > temp.bed
chr=$(echo $scfr | awk '{print$1}')
gr=$(readlink -f /media/aswin/SCFR/SCFR-main/genome_reports/* | grep "$species")
cs=$(awk -F "\t" -v a="$chr" '$9==a {print$11}' $gr)

/media/aswin/SCFR/SCFR-main/chrs/human/"$chr".fasta
awk '$3=="CDS"' /media/aswin/SCFR/SCFR-main/genes/human/GCF_009914755.1_T2T-CHM13v2.0_genomic.gtf | grep KDM5D | grep 20642300 | awk -v a="$cs" '{if($7=="+") print$0,$4-($4-($4%3))+1; else if($7=="-") print$0,(a-$5+1)-(a-$5+1)-((a-$5+1)%3)-$8}'


#Total number of exons
tnoe=$(awk '{print$7,$8,$9,$10,$11,$12}' OFS="\t" temp.bed | sort -u | wc -l)
#Overlapping exons
mc=$(awk '{print$7,$8,$9,$10,$11,$12}' OFS="\t" temp.bed | bedtools sort -i - | bedtools merge -i - -s -c 5 -o count | wc -l)
mc=$(bedtools sort -i temp.bed | bedtools merge -i - | wc -l)
if [[ $tnoe == "1" ]] || [[ $mc == "1" ]]
then
cat temp.bed | awk '{if($6~"-") print$0,$3-$9,$8-$2; else print$0,$8-$2,$3-$9}' OFS="\t" | awk '!seen[$7,$8,$9,$10,$12,$13,$17,$19,$20,$21,$22]++' >> single_exon.tsv



#####

head -n 1000 human_scfr_containing_cds.bed | awk 'BEGIN{FS="\t"; OFS="\t"} {print$1,$2,$3,$4,$5,$6}' | awk -v a="$cs" '{if($6=="+") print$0,$2-($2-($2%3))+1; else if($6=="-") print$0,"-"(((a-$3+1)-((a-$3+1)-((a-$3+1)%3))+2)%3)+1}' | awk '{if($4==$7) print$0,"same"; else print$0,"diff"}' | less

#correct
awk '$3=="CDS"' /media/aswin/SCFR/SCFR-main/genes/human/GCF_009914755.1_T2T-CHM13v2.0_genomic.gtf | grep APOBEC1 | grep 7667091 | awk -v a="$cs" '{if($7=="+") print$0,$4-($4-($4%3))+1; else if($7=="-") print$0,"-"(((a-$5+1)-((a-$5+1)-((a-$5+1)%3))+2)%3)+1+$8}'
#correct
awk '$3=="CDS"' /media/aswin/SCFR/SCFR-main/genes/human/GCF_009914755.1_T2T-CHM13v2.0_genomic.gtf | grep APOBEC1 | grep 7667091 | awk -v a="$cs" '{if($7=="+") print$0,$4-($4-($4%3))+1; else if($7=="-") print$0,(a-$5+1)-(a-$5+1)-((a-$5+1)%3)-$8}'
#incorrect
#awk '$3=="CDS"' /media/aswin/SCFR/SCFR-main/genes/human/GCF_009914755.1_T2T-CHM13v2.0_genomic.gtf | grep -w CHD4 | grep 6612859 | awk -v a="$cs" '{if($7=="+") print$0,$4-($4-($4%3))+1; else if($7=="-") print$0,"-"(((a-$5+1)-((a-$5+1)-((a-$5+1)%3))+2)%3)+1+$8}'
#correct
awk '$3=="CDS"' /media/aswin/SCFR/SCFR-main/genes/human/GCF_009914755.1_T2T-CHM13v2.0_genomic.gtf | grep -w CHD4 | grep 6612859 | awk -v a="$cs" '{if($7=="+") print$0,$4-($4-($4%3))+1; else if($7=="-") print$0,(a-$5+1)-(a-$5+1)-((a-$5+1)%3)-$8}'
#correct
awk '$3=="CDS"' /media/aswin/SCFR/SCFR-main/genes/human/GCF_009914755.1_T2T-CHM13v2.0_genomic.gtf | grep -w CHD4 | grep 6612859 | awk -v a="$cs" '{if($7=="+") print$0,$4-($4-($4%3))+1; else if($7=="-") print$0,"-"(a-$5+1)-(a-$5+1)-((a-$5+1)%3)+$8}'
#correct
awk '$3=="CDS"' /media/aswin/SCFR/SCFR-main/genes/human/GCF_009914755.1_T2T-CHM13v2.0_genomic.gtf | grep -w TNFRSF1A | grep 6339687 | awk -v a="$cs" '{if($7=="+") print$0,$4-($4-($4%3))+1; else if($7=="-") print$0,"-"(((a-$5+1)-((a-$5+1)-((a-$5+1)%3))+2)%3)+1+$8}'
#incorrect
#awk '$3=="CDS"' /media/aswin/SCFR/SCFR-main/genes/human/GCF_009914755.1_T2T-CHM13v2.0_genomic.gtf | grep -w TNFRSF1A | grep 6339687 | awk -v a="$cs" '{if($7=="+") print$0,$4-($4-($4%3))+1; else if($7=="-") print$0,"-"(a-$5+1)-(a-$5+1)-((a-$5+1)%3)+$8}'
#correct
awk '$3=="CDS"' /media/aswin/SCFR/SCFR-main/genes/human/GCF_009914755.1_T2T-CHM13v2.0_genomic.gtf | grep -w TNFRSF1A | grep 6339687 | awk -v a="$cs" '{if($7=="+") print$0,$4-($4-($4%3))+1; else if($7=="-") print$0,(a-$5+1)-(a-$5+1)-((a-$5+1)%3)-$8}'
#correct
awk '$3=="CDS"' /media/aswin/SCFR/SCFR-main/genes/human/GCF_009914755.1_T2T-CHM13v2.0_genomic.gtf | grep -w SCNN1A | grep 6373769 | awk -v a="$cs" '{if($7=="+") print$0,$4-($4-($4%3))+1; else if($7=="-") print$0,(a-$5+1)-(a-$5+1)-((a-$5+1)%3)-$8}'
#correct
awk '$3=="CDS"' /media/aswin/SCFR/SCFR-main/genes/human/GCF_009914755.1_T2T-CHM13v2.0_genomic.gtf | grep -w SCNN1A | grep 6371663 | awk -v a="$cs" '{if($7=="+") print$0,$4-($4-($4%3))+1; else if($7=="-") print$0,(a-$5+1)-(a-$5+1)-((a-$5+1)%3)-$8}'
#correct
awk '$3=="CDS"' /media/aswin/SCFR/SCFR-main/genes/human/GCF_009914755.1_T2T-CHM13v2.0_genomic.gtf | grep -w MYL7 | grep 44298123 | awk -v a="$cs" '{if($7=="+") print$0,$4-($4-($4%3))+1; else if($7=="-") print$0,(a-$5+1)-(a-$5+1)-((a-$5+1)%3)-$8}'
#correct
awk '$3=="CDS"' /media/aswin/SCFR/SCFR-main/genes/human/GCF_009914755.1_T2T-CHM13v2.0_genomic.gtf | grep -w MYL7 | grep 44299629 | awk -v a="$cs" '{if($7=="+") print$0,$4-($4-($4%3))+1; else if($7=="-") print$0,(a-$5+1)-(a-$5+1)-((a-$5+1)%3)-$8}'
#correct
awk '$3=="CDS"' /media/aswin/SCFR/SCFR-main/genes/human/GCF_009914755.1_T2T-CHM13v2.0_genomic.gtf | grep KDM5D | grep 20642300 | awk -v a="$cs" '{if($7=="+") print$0,$4-($4-($4%3))+1; else if($7=="-") print$0,(a-$5+1)-(a-$5+1)-((a-$5+1)%3)-$8}'

#Plus strand
#correct m1
awk '$3=="CDS"' /media/aswin/SCFR/SCFR-main/genes/human/GCF_009914755.1_T2T-CHM13v2.0_genomic.gtf | grep SAMD11 | grep 369889 | awk -v a="$cs" '{if($7=="+") print$0,$4-($4-($4%3))+$8; else if($7=="-") print$0,(a-$5+1)-(a-$5+1)-((a-$5+1)%3)-$8}'
#incorrect m2
#awk '$3=="CDS"' /media/aswin/SCFR/SCFR-main/genes/human/GCF_009914755.1_T2T-CHM13v2.0_genomic.gtf | grep SAMD11 | grep 369889 | awk -v a="$cs" '{if($7=="+") print$0,$4+$8-(($4+$8)-(($4+$8)%3)); else if($7=="-") print$0,(a-$5+1)-(a-$5+1)-((a-$5+1)%3)-$8}'
#incorrect m3
#awk '$3=="CDS"' /media/aswin/SCFR/SCFR-main/genes/human/GCF_009914755.1_T2T-CHM13v2.0_genomic.gtf | grep SAMD11 | grep 369889 | awk -v a="$cs" '{if($7=="+") print$0,$4-($4-($4%3))+1-$8; else if($7=="-") print$0,(a-$5+1)-(a-$5+1)-((a-$5+1)%3)-$8}'
#incorrect m1
#awk '$3=="CDS"' /media/aswin/SCFR/SCFR-main/genes/human/GCF_009914755.1_T2T-CHM13v2.0_genomic.gtf | tail -n 1000 | grep 22396565 | awk -v a="$cs" '{if($7=="+") print$0,$4-($4-($4%3))+$8; else if($7=="-") print$0,(a-$5+1)-(a-$5+1)-((a-$5+1)%3)-$8}'
#correct m2
awk '$3=="CDS"' /media/aswin/SCFR/SCFR-main/genes/human/GCF_009914755.1_T2T-CHM13v2.0_genomic.gtf | tail -n 1000 | grep 22396565 | awk -v a="$cs" '{if($7=="+") print$0,$4+$8-(($4+$8)-(($4+$8)%3)); else if($7=="-") print$0,(a-$5+1)-(a-$5+1)-((a-$5+1)%3)-$8}'
#correct m3
awk '$3=="CDS"' /media/aswin/SCFR/SCFR-main/genes/human/GCF_009914755.1_T2T-CHM13v2.0_genomic.gtf | tail -n 1000 | grep 22396565 | awk -v a="$cs" '{if($7=="+") print$0,$4-($4-($4%3))+1-$8; else if($7=="-") print$0,(a-$5+1)-(a-$5+1)-((a-$5+1)%3)-$8}'
awk '$3=="CDS"' /media/aswin/SCFR/SCFR-main/genes/human/GCF_009914755.1_T2T-CHM13v2.0_genomic.gtf | grep USP9Y | grep 13644275 | awk -v a="$cs" '{if($7=="+") print$0,$4-($4-($4%3))+$8}'
awk '$3=="CDS"' /media/aswin/SCFR/SCFR-main/genes/human/GCF_009914755.1_T2T-CHM13v2.0_genomic.gtf | grep USP9Y | grep 13644688 | awk -v a="$cs" '{if($7=="+") print$0,$4-($4-($4%3))+$8}'
awk '$3=="CDS"' /media/aswin/SCFR/SCFR-main/genes/human/GCF_009914755.1_T2T-CHM13v2.0_genomic.gtf | grep CD99 | awk -v a="$cs" '{if($7=="+") print$0,$4-($4-($4%3))+$8}' | grep 2391385 
awk '$3=="CDS"' /media/aswin/SCFR/SCFR-main/genes/human/GCF_009914755.1_T2T-CHM13v2.0_genomic.gtf | grep CD99 | awk -v a="$cs" '{if($7=="+") print$0,$4-($4-($4%3))+$8}' | grep 2394572
awk '$3=="CDS"' /media/aswin/SCFR/SCFR-main/genes/human/GCF_009914755.1_T2T-CHM13v2.0_genomic.gtf | grep CD99 | awk -v a="$cs" '{if($7=="+") print$0,$4-($4-($4%3))+$8}' | egrep "2415161|2415158"
awk '$3=="CDS"' /media/aswin/SCFR/SCFR-main/genes/human/GCF_009914755.1_T2T-CHM13v2.0_genomic.gtf | grep KIF2C | awk -v a="$cs" '{if($7=="+") print$0,$4-($4-($4%3))+$8}' | grep 30203
awk '$3=="CDS"' /media/aswin/SCFR/SCFR-main/genes/human/GCF_009914755.1_T2T-CHM13v2.0_genomic.gtf | grep -w CYP4X1 | awk -v a="$cs" '{if($7=="+") print$0,$4-($4-($4%3))+$8}' | grep 10644
awk '$3=="CDS"' /media/aswin/SCFR/SCFR-main/genes/human/GCF_009914755.1_T2T-CHM13v2.0_genomic.gtf | grep -w PRKAA2 | awk -v a="$cs" '{if($7=="+") print$0,$4-($4-($4%3))+$8}' | grep 70249
awk '$3=="CDS"' /media/aswin/SCFR/SCFR-main/genes/human/GCF_009914755.1_T2T-CHM13v2.0_genomic.gtf | grep -w CACHD1 | awk -v a="$cs" '{if($7=="+") print$0,$4-($4-($4%3))}' | grep 10710
awk '$3=="CDS"' /media/aswin/SCFR/SCFR-main/genes/human/GCF_009914755.1_T2T-CHM13v2.0_genomic.gtf | grep -w DAP3 | awk -v a="$cs" '{if($7=="+") print$0,$4-($4-($4%3))+$8}' | grep 64255
awk '$3=="CDS"' /media/aswin/SCFR/SCFR-main/genes/human/GCF_009914755.1_T2T-CHM13v2.0_genomic.gtf | grep -w ATP1A4 | awk -v a="$cs" '{if($7=="+") print$0,$4-($4-($4%3))+$8}' | grep 08333




#All bed entries
awk '$3=="CDS"' /media/aswin/SCFR/SCFR-main/genes/human/GCF_009914755.1_T2T-CHM13v2.0_genomic.gtf | awk -v a="$cs" '{if($7=="+") print$0,$4-($4-($4%3))}' | awk -v a="$cs" '{if($7=="+") print$0,$4-($4-($4%3))+$8}' | awk '{print$1,$4,$5,$10,$12,$7,$8,$(NF-1),$NF}' | less -S
awk '$3=="CDS"' /media/aswin/SCFR/SCFR-main/genes/human/GCF_009914755.1_T2T-CHM13v2.0_genomic.gtf | awk -v a="$cs" '{if($7=="+") print$0,$4-($4-($4%3))}' | awk -v a="$cs" '{if($7=="+") print$0,$4-($4-($4%3))+$8}' | awk '{print$1,$4,$5,$10,$12,$7,$8,$(NF-1),$NF}' | sort -k7,7n -k8,8n -k9,9n | less
