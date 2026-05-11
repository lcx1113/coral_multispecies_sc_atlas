input_dir="/home/users/Ruijie/tools/samap_directory/extended_longest_species_fasta/cleaned_noSpeciesPrefix"

files=(${input_dir}/*_cleaned_modified.longest.fasta)

for ((i=0; i<${#files[@]}; i++)); do
  for ((j=i+1; j<${#files[@]}; j++)); do

    file1="${files[i]}"
    file2="${files[j]}"

    sp1=$(basename "$file1" | cut -d'_' -f1)
    sp2=$(basename "$file2" | cut -d'_' -f1)

    id1=$(echo "$sp1" | cut -c1-2 | tr 'A-Z' 'a-z')
    id2=$(echo "$sp2" | cut -c1-2 | tr 'A-Z' 'a-z')

    echo "Running $sp1 vs $sp2"

    bash extended_longest_map_genes.sh \
      --tr1 "$file1" --t1 prot --n1 "$id1" \
      --tr2 "$file2" --t2 prot --n2 "$id2"

  done
done
