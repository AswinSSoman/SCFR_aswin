

cd /media/aswin/SCFR/SCFR-main
time for seq in $(find Length_threshold_PCA_kmeans/ -mindepth 3 -name "*_gt1000.fasta" -type f | xargs readlink -f)
do
mkdir -p repeat_masker/"$species"
species=$(echo $seq | cut -f7 -d "/")
sp=$(grep -i "$species" /media/aswin/SCFR/SCFR-main/genome_reports/species_names | awk '{print$2}' | tr "_" " ")
echo ">"$species
cd repeat_masker
time /media/aswin/programs/RepeatMasker-4.1.7/RepeatMasker/RepeatMasker -pa 8 -a -s -u -gff -html -species "$sp" -dir $species -xsmall $seq
unset species sp
cd /media/aswin/SCFR/SCFR-main
done
