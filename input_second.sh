#!/bin/bash


# Usage: ./input_second.sh /path/to/your_fav_SAM.sam /path/to/out/dir/
# Creates all files in './' by default, or in '$2'

if [ -z "$2" ]; then
    outDir=./
else
    outDir=$2
fi


base_used=toZymo # To replace if using a different  
to_metadat_file=metadat_files/$base_used'_tax_metadat.tsv'
if [[ -f $to_metadat_file ]]; then 
    echo;
else
    echo "ERROR corresponding taxonomic metadata file NOT FOUND"
    exit 1
fi

./2-prim_analysis.py -i $1 -l species -w true -o $outDir

base_name=$(basename $1 .sam)
tmp_filename=$base_name'_MINOR_RM_LCA'

echo Creating OTUs table 'for:' $tmp_filename.. 
echo With $to_metadat_file (to produce OTU biom table)
cat $tmp_filename'.map' | sed '$ d' > $outDir/$tmp_filename'.map2'
make_otu_table.py -i $outDir/$tmp_filename'.map2' -t $to_metadat_file \
                                                -o $outDir/$tmp_filename'.biom'
echo Summarizing taxonomy'...'
summarize_taxa.py -i $outDir/$tmp_filename'.biom' -o $outDir \
                    --level 7 --absolute_abundance --suppress_biom_table_output

# Change a bit the header:
sed -i 's/\#OTU\sID/superkingdom;phylum;class;order;family;genus;species/' \
                                                        $outDir/$tmp_filename'_L7.txt'
echo "FINISHED !" ; echo ""
