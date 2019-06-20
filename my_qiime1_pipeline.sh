#!/bin/bash

source /home/sheldon/Applications/miniconda3/etc/profile.d/conda.sh
conda activate qiime1

my_list=""
#base_used=toSilva
base_used=toZymo
#base_used=toNewlot_zymo
#base_used=toRrn
to_metadat_file=metadat_files/$base_used'_tax_metadat.tsv'
OTUs_outDir=OTUs_maps_n_tables

#for file_path in ./OTUs_maps_n_tables/*qFilt*$base_used*map2
for file_path in ./OTUs_maps_n_tables/*$base_used*map2
do
    base_filename=$(basename $file_path .map2)
    echo Creating OTUs table 'for:' $base_filename.. with $to_metadat_file
    out_file=$OTUs_outDir/$base_filename.biom
    make_otu_table.py -i $file_path -t $to_metadat_file -o $out_file
    summarize_taxa.py -i $out_file -o $OTUs_outDir --level 7 --absolute_abundance --suppress_biom_table_output
    
    # Change a bit the header:
    sed -i 's/\#OTU\sID/superkingdom;phylum;class;order;family;genus;species/' \
                                        $OTUs_outDir/$base_filename'_L7.txt'
    #    > $OTUs_outDir/$base_filename'_L7.txt2'
    #rm $OTUs_outDir/$base_filename'_L7.txt'
    #mv $OTUs_outDir/$base_filename'_L7.txt2' $OTUs_outDir/$base_filename'_L7.txt'
    
    # For Krona charts:
    #cat $OTUs_outDir/$base_filename'_L7.txt' | \
    #        awk -F "\t" 'NR>1 {first=$1; $1=""; print $0"\t"first;}' | \
    #        sed 's/\#OTU\sID/domain;phylum;class;order;family;genus;species/' \
    #                > $OTUs_outDir/$base_filename'_L7.krona'
                    
    my_list+="$out_file,"
done
my_list=$(echo $my_list | sed 's/,$//')

#exit 1

echo ""
echo Merging OTUs tables..
to_merged_biom=$OTUs_outDir/$base_used.biom
echo WROTE": $to_merged_biom"
#merge_otu_tables.py -i $my_list",$OTUs_outDir/ref_oldLot.biom" -o $OTUs_outDir/'qFilt_'$base_used.biom
merge_otu_tables.py -i $my_list",$OTUs_outDir/ref_oldLot.biom" -o $to_merged_biom
echo $(echo $my_list | sed 's/,/ /')

echo ""
echo Converting biom table into spf format..
summarize_taxa.py -i $to_merged_biom -o $OTUs_outDir --level 7 -a --suppress_biom_table_output
../solve_SILVA.py $OTUs_outDir/$base_used"_L7.txt"
#mv $base_used'_L7.biom' 'old_'$base_used'_L7.biom' 
#cat 'old_'$base_used'_L7.biom' | sed "s/\b\OTU ID\b/taxonomy/" | sed "s/#//" > $base_used'_L7.biom'
exit 1

#python ../biom_to_stamp.py -m taxonomy $to_merged_biom > $OTUs_outDir/$base_used'.spf0'
#python ../checkHierarchy.py $OTUs_outDir/$base_used'.spf0'
#rm $OTUs_outDir/$base_used'.spf0'

biom convert -i $OTUs_outDir/$base_used'.biom' -o felix.tsv --to-tsv
../solve_SILVA.py "$to_merged_biom"
rm felix.tsv

exit 1


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
