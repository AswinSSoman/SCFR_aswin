

#Run PCA & k-mean clustering for different length bins & see if clustering gets better with SCFR length bins
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
      lo["1000_2500"]=1000; hi["1000_2500"]=2500
      lo["2500_5000"]=2500; hi["2500_5000"]=5000
      lo["5000_7500"]=5000; hi["5000_7500"]=7500
      lo["7500_10000"]=7500; hi["7500_10000"]=10000
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

  for lenbin in 1000_2500 2500_5000 5000_7500 7500_10000
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

    python3 scripts/rscu_metrics.py "$OUTDIR" "$OUTDIR"
    Rscript my_scripts/plotPCA_kmeans_loadings.R "$OUTDIR" "$OUTDIR"
  done
}
export -f process_species

#Run for 7 species
time parallel -j 7 process_species ::: \
  human bonobo chimpanzee gorilla borangutan sorangutan gibbon


