#!/bin/bash

source /home/sheldon/Applications/miniconda3/etc/profile.d/conda.sh
conda activate qiime1

my_list=""
#base_used=toSilva
base_used=toRrn
to_metadat_file=metadat_files/$base_used'_tax_metadat.tsv'
OTUs_outDir=OTUs_maps_n_tables

for file_path in ./OTUs_maps_n_tables/*$base_used*map
do
    base_filename=$(basename $file_path .map)
    echo Creating OTUs table 'for:' $base_filename.. with $to_metadat_file
    out_file=$OTUs_outDir/$base_filename.biom
    make_otu_table.py -i $file_path -t $to_metadat_file -o $out_file
    summarize_taxa.py -i $out_file -o $OTUs_outDir --level 7 --absolute_abundance --suppress_biom_table_output
    my_list+="$out_file,"
done
my_list=$(echo $my_list | sed 's/,$//')

exit 1

echo ""
echo Merging OTUs tables..
merge_otu_tables.py -i $my_list -o $OTUs_outDir/$base_used.biom
#exit 1

merge_otu_tables.py -i $OTUs_outDir/toRrn.biom,$OTUs_outDir/toSilva.biom \
                                                    -o $OTUs_outDir/all.biom

echo ""
echo Calculating alpha-diversity..
alpha_diversity.py -i $OTUs_outDir/$base_used.biom -o alpha-div.txt \
                                                        -m observed_otus,chao1
cat alpha-div.txt

echo ""
echo Calculating alpha-rarefaction..
rm -r test

#alpha_rarefaction.py -i $OTUs_outDir/$base_used.biom \
#echo "alpha_rarefaction.py -i $OTUs_outDir/$base_used.biom \
alpha_rarefaction.py -i $OTUs_outDir/all.biom \
                        -m sampl_mappping.map -o test \
                        --parameter_fp param_for_plot_with_single_sampl.txt \
                        --parallel --jobs_to_start 15
                        #--retain_intermediate_files


#conda deactivate
