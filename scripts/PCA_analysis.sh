#######################################################################################################################################################################################################################################################################################################
#PCA
#######################################################################################################################################################################################################################################################################################################


#######################################################################################################################################################################################################################################################################################################
#Length bin:

#Run PCA & k-mean clustering for different length bins & see if clustering gets better with SCFR length bins (189m8.519s)
process_species() {
  species="$1"

  echo ">$species"

  SCFR="SCFR_all/${species}_SCFR_all.out"
  GENOME_FASTA="genomes/${species}/GC*.fna"
  BASEOUT="/media/aswin/SCFR/SCFR-main/Length_bin_PCA_kmeans/${species}"

  mkdir -p "$BASEOUT"

  awk '
    BEGIN {
      OFS = "\t"

      lo["1000_2500"]=1000;   hi["1000_2500"]=2500
      lo["2500_5000"]=2500;   hi["2500_5000"]=5000
      lo["5000_7500"]=5000;   hi["5000_7500"]=7500
      lo["7500_10000"]=7500;  hi["7500_10000"]=10000

      # NEW BIN (>10000)
      lo["gt10000"]=10001;    hi["gt10000"]=1e12
    }
    {
      len = $3 - $2
      strand = ($4 ~ "-") ? "-" : "+"

      for (b in lo) {
        if (len >= lo[b] && len <= hi[b]) {
          print $1, $2, $3, $4, 1, strand > "'"$BASEOUT"'/" b ".bed"
          break
        }
      }
    }
  ' "$SCFR"

  for lenbin in 1000_2500 2500_5000 5000_7500 7500_10000 gt10000
  do
    OUTDIR="$BASEOUT/$lenbin"
    mkdir -p "$OUTDIR"

    BED="$BASEOUT/$lenbin.bed"
    FASTA="$OUTDIR/${species}_${lenbin}.fasta"

    [[ ! -s "$BED" ]] && continue

    sort -k1,1 -k2,2n "$BED" \
    | bedtools getfasta \
        -fi $GENOME_FASTA \
        -bed stdin \
        -name+ -s \
    > "$FASTA"
python3 /media/aswin/SCFR/SCFR-main/my_scripts/PCA/corrected_rscu_calc.py "$OUTDIR" "$OUTDIR"
Rscript /media/aswin/SCFR/SCFR-main/my_scripts/PCA/pca_script.R "$OUTDIR" "$OUTDIR"
Rscript /media/aswin/SCFR/SCFR-main/my_scripts/PCA/kmeans_script.R "$OUTDIR" "$OUTDIR"
Rscript /media/aswin/SCFR/SCFR-main/my_scripts/PCA/plot_by_cluster.R "$OUTDIR" "$OUTDIR"
Rscript /media/aswin/SCFR/SCFR-main/my_scripts/PCA/plot_color_by_clade.R "$OUTDIR" "$OUTDIR"
  done
}
export -f process_species

#Run for 7 species
#time parallel -j 7 process_species ::: human bonobo chimpanzee gorilla borangutan sorangutan gibbon
time for species in human bonobo chimpanzee gorilla borangutan sorangutan gibbon; do
    process_species "$species"
done

#Check run
cd /media/aswin/SCFR/SCFR-main/Length_bin_PCA_kmeans
for i in $(ls -d */ | tr -d "/")
do
for f in 1000_2500 2500_5000 5000_7500 7500_10000 gt10000
do
pdf=$(find $i/$f -mindepth 1 -maxdepth 1 -name "pca_cluster_PC1_PC2.pdf" -type f)
if [[ -z $pdf ]]; then p="not-run"; else p="run"; fi
ns=$(grep ">" -c $i/$f/$i"_"$f".fasta")
ts=$(grep -v ">" $i/$f/$i"_"$f".fasta" | wc | awk '{print$3-$1}')
echo $i $f $p $ns $ts
unset pdf p ns ts
done
unset f
done | column -t

#######################################################################################################################################################################################################################################################################################################
#Length threshold: (243m58.707s)

process_species_threshold() {
  species="$1"

  echo ">$species"

  SCFR="SCFR_all/${species}_SCFR_all.out"
  GENOME_FASTA="genomes/${species}/GC*.fna"
  BASEOUT="/media/aswin/SCFR/SCFR-main/Length_threshold_PCA_kmeans/${species}"

  mkdir -p "$BASEOUT"

  awk '
    BEGIN {
      OFS = "\t"
      thr["gt1000"]=1000
      thr["gt2500"]=2500
      thr["gt5000"]=5000
      thr["gt7500"]=7500
      thr["gt10000"]=10000
    }
    {
      len = $3 - $2
      strand = ($4 ~ "-") ? "-" : "+"

      for (t in thr) {
        if (len > thr[t]) {
          print $1, $2, $3, $4, 1, strand > "'"$BASEOUT"'/" t ".bed"
        }
      }
    }
  ' "$SCFR"

  for threshold in gt1000 gt2500 gt5000 gt7500 gt10000
  do
    OUTDIR="$BASEOUT/$threshold"
    mkdir -p "$OUTDIR"

    BED="$BASEOUT/$threshold.bed"
    FASTA="$OUTDIR/${species}_${threshold}.fasta"

    [[ ! -s "$BED" ]] && continue

    sort -k1,1 -k2,2n "$BED" \
    | bedtools getfasta \
        -fi $GENOME_FASTA \
        -bed stdin \
        -name+ -s \
    > "$FASTA"
python3 /media/aswin/SCFR/SCFR-main/my_scripts/PCA/corrected_rscu_calc.py "$OUTDIR" "$OUTDIR"
Rscript /media/aswin/SCFR/SCFR-main/my_scripts/PCA/pca_script.R "$OUTDIR" "$OUTDIR"
Rscript /media/aswin/SCFR/SCFR-main/my_scripts/PCA/kmeans_script.R "$OUTDIR" "$OUTDIR"
Rscript /media/aswin/SCFR/SCFR-main/my_scripts/PCA/plot_by_cluster.R "$OUTDIR" "$OUTDIR"
Rscript /media/aswin/SCFR/SCFR-main/my_scripts/PCA/plot_color_by_clade.R "$OUTDIR" "$OUTDIR"
  done
}
export -f process_species_threshold

#Run for 7 species
#time parallel -j 7 process_species_threshold ::: human bonobo chimpanzee gorilla borangutan sorangutan gibbon
time for species in human bonobo chimpanzee gorilla borangutan sorangutan gibbon; do
    process_species_threshold "$species"
done

#Check run
cd /media/aswin/SCFR/SCFR-main/Length_threshold_PCA_kmeans
for i in $(ls -d */ | tr -d "/")
do
for f in gt1000 gt2500 gt5000 gt7500 gt10000
do
pdf=$(find $i/$f -mindepth 1 -maxdepth 1 -name "pca_cluster_PC1_PC2.pdf" -type f)
if [[ -z $pdf ]]; then p="not-run"; else p="run"; fi
ns=$(grep ">" -c $i/$f/$i"_"$f".fasta")
ts=$(grep -v ">" $i/$f/$i"_"$f".fasta" | wc | awk '{print$3-$1}')
echo $i $f $p $ns $ts
unset pdf p ns ts
done
unset f
done | column -t

#######################################################################################################################################################################################################################################################################################################
#label coding vs non-coding SCFRs in PCA

cd /media/aswin/SCFR/SCFR-main
for bed in $(find Length_threshold_PCA_kmeans/ -mindepth 2 -maxdepth 2 -name "gt1000.bed" -type f | xargs readlink -f)
do
species=$(echo $bed | cut -f7 -d "/")
sp=$(grep -i "$species" /media/aswin/SCFR/SCFR-main/genome_reports/species_names | awk '{print$2}' | tr "_" " ")
path=$(echo $bed | sed 's!/gt1000.bed!!g')
cd $path
echo ">"$species
#scfr overlapping exons in same frame
bedtools intersect -a $bed -b /media/aswin/SCFR/SCFR-main/exon_shadow/$species/"$species"_coding_exons.bed -wo | awk '$12$13==$4' | awk '{$NF = $NF + 1; print}' | awk '{print$0,$9-$8+1,$3-$2+1,($20/($9-$8+1))*100}' \
| awk '{print$4"::"$1":"$2"-"$3"("$6")","in_frame_coding",$20,$21,$22,$23}' | sort -k1,1 -k6,6rn | awk '!seen[$1]++' | sed '1i SCFR coding_status overlap_len cds_len scfr_len percent_coding_in_scfr' | sed 's/[ ]\+/\t/g' > scfr_inframe_coding_gt1000.tsv
#scfr overlapping exons in diff frame
bedtools intersect -a $bed -b /media/aswin/SCFR/SCFR-main/exon_shadow/$species/"$species"_coding_exons.bed -wo | awk '$12$13!=$4' | awk '{$NF = $NF + 1; print}' | awk '{print$0,$9-$8+1,$3-$2+1,($20/($9-$8+1))*100}' \
| awk '{print$4"::"$1":"$2"-"$3"("$6")","out_frame_coding",$20,$21,$22,$23}' | sort -k1,1 -k6,6rn | awk '!seen[$1]++' | sed '1i SCFR coding_status overlap_len cds_len scfr_len percent_coding_in_scfr' | sed 's/[ ]\+/\t/g' > scfr_outframe_coding_gt1000.tsv
#scfr_inframe_coding_gt1000.tsv & scfr_outframe_coding_gt1000.tsv have some common scfrs, remove these common ones from scfr_outframe_coding_gt1000.tsv
awk -F'\t' 'NR==FNR {seen[$1]; next} !($1 in seen)' scfr_inframe_coding_gt1000.tsv scfr_outframe_coding_gt1000.tsv | sed '1i SCFR coding_status overlap_len cds_len scfr_len percent_coding_in_scfr' | sed 's/[ ]\+/\t/g'  > scfr_outframe_coding_gt1000_filtered.tsv
#scfr not overlapping exons
awk '{print $1}' scfr_inframe_coding_gt1000.tsv scfr_outframe_coding_gt1000_filtered.tsv > exclude_list.txt
seqkit grep -v -f exclude_list.txt gt1000/"$species"_gt1000.fasta | grep ">" | sed 's/$/ non_coding/g' | sed '1i SCFR coding_status' | sed 's/[ ]\+/\t/g' > scfr_non_coding_gt1000.tsv
#final list of scfrs & coding status
awk '{print$1,$2}' scfr_inframe_coding_gt1000.tsv scfr_outframe_coding_gt1000_filtered.tsv scfr_non_coding_gt1000.tsv | grep -v "coding_status" | sed '1i SCFR coding_status' | sed 's/[ ]\+/\t/g' > scfr_coding_status.tsv
unset species sp path 
cd /media/aswin/SCFR/SCFR-main
done

#######################################################################################################################################################################################################################################################################################################
#Summary

#Summary of Length bin
cd /media/aswin/SCFR/SCFR-main/
time for species in human bonobo borangutan sorangutan chimpanzee gorilla gibbon 
do
cd Length_bin_PCA_kmeans/"$species"/
for lenbin in 2500_5000 5000_7500 7500_10000 gt10000
do
cd $lenbin
i1=$(awk '$1 == "PC1" || $1 == "PC2" {sum += $3} END {print sum}' explained_variance.tsv)
i2=$(awk 'NR>1 && $5 != "NA" {if($5 > max) {max=$5; k=$0}} END {print k}' k_optimization_scores.tsv)
echo $species $lenbin $i1 $i2
unset i1 i2
cd ../
done
unset lenbin
cd /media/aswin/SCFR/SCFR-main/
done | sed '1i Species Length_bin PC1_PC2 k Silhouette DBI WCSS Curvature' | sed 's/[ ]\+/\t/g' > /media/aswin/SCFR/SCFR-main/Length_bin_PCA_kmeans/all_species_pca_clustering_summary.tsv

#Summary of Length threshold
cd /media/aswin/SCFR/SCFR-main/
time for species in human bonobo borangutan sorangutan chimpanzee gorilla gibbon 
do
cd Length_threshold_PCA_kmeans/"$species"/
for lenthr in gt2500 gt5000 gt7500 gt10000
do
cd $lenthr
i1=$(awk '$1 == "PC1" || $1 == "PC2" {sum += $3} END {print sum}' explained_variance.tsv)
i2=$(awk 'NR>1 && $5 != "NA" {if($5 > max) {max=$5; k=$0}} END {print k}' k_optimization_scores.tsv)
echo $species $lenthr $i1 $i2
unset i1 i2
cd ../
done
unset lenthr
cd /media/aswin/SCFR/SCFR-main/
done | sed '1i Species Length_threshold PC1_PC2 k Silhouette DBI WCSS Curvature' | sed 's/[ ]\+/\t/g' > /media/aswin/SCFR/SCFR-main/Length_threshold_PCA_kmeans/all_species_pca_clustering_summary.tsv

/media/aswin/SCFR/SCFR-main/Length_threshold_PCA_kmeans$ Rscript ../Length_bin_PCA_kmeans/plot_pca_clustering_summary.R all_species_pca_clustering_summary.tsv all_species_pca_clustering_summary.png





/media/aswin/SCFR/SCFR-main/exon_shadow/gibbon/gibbon_coding_exons.bed
bedtools intersect -a /media/aswin/SCFR/SCFR-main/exon_shadow/gibbon/gibbon_coding_exons.bed -b 1000_2500.bed -wo > test



cd /media/aswin/SCFR/SCFR-main/
time for species in human bonobo borangutan sorangutan chimpanzee gorilla gibbon 
do
for f in 1000_2500 2500_5000 5000_7500 7500_10000 gt10000
do
cd Length_bin_PCA_kmeans/"$species"/"$f"
bedtools intersect -a ../"$f".bed -b /media/aswin/SCFR/SCFR-main/exon_shadow/$species/"$species"_coding_exons.bed -wo \
 | awk '$6$7==$17' | awk '{print$0,$3-$2,$16-$15,($20/($16-$15))*100}' | awk '{print$17"::"$14":"$15"-"$16"("$19")","in_frame_coding",$20,$21,$22,$23}' | awk 'NR==1; NR>1 {print $0}' | sort -k1,1 -k6,6rn | awk '!seen[$1]++' | sed '1i SCFR coding_status overlap_len cds_len scfr_len percent_coding_in_scfr' | sed 's/[ ]\+/\t/g' > scfr_coding_"$f".tsv
bedtools intersect -a ../"$f".bed -b /media/aswin/SCFR/SCFR-main/exon_shadow/$species/"$species"_coding_exons.bed -wo \
 | awk '$6$7!=$17' | awk '{print$0,$3-$2,$16-$15,($20/($16-$15))*100}' | awk '{print$17"::"$14":"$15"-"$16"("$19")","out_frame_coding",$20,$21,$22,$23}' | awk 'NR==1; NR>1 {print $0}' test | sort -k1,1 -k6,6rn | awk '!seen[$1]++' | sed '1i SCFR coding_status overlap_len cds_len scfr_len percent_coding_in_scfr' | sed 's/[ ]\+/\t/g' > scfr_non_coding_"$f".tsv




find . -maxdepth 3 -mindepth 3  -name "*.fasta" -type f | xargs -n1 bash -c 'paste <(echo $0 | cut -f4 -d "/") <(grep -v ">" $0 | wc | awk "{print\$3-\$1}")' > total_fasta_length

