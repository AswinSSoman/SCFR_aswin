

cd /media/aswin/SCFR/SCFR-main
time for seq in $(find Length_threshold_PCA_kmeans/ -mindepth 3 -name "*_gt1000.fasta" -type f)
do
mkdir -p repeat_masker/"$species"
species=$(echo $seq | cut -f2 -d "/")
echo ">"$species
#sp=$(grep -i "$species" /media/aswin/SCFR/SCFR-main/genome_reports/species_names | awk '{print$1}' | tr "_" " ")
cd repeat_masker
time /media/aswin/programs/RepeatMasker-4.1.7/RepeatMasker/RepeatMasker -pa 8 -a -s -u -gff -html -species "$species" -dir $species -xsmall $seq
unset species
cd /media/aswin/SCFR/SCFR-main
done
