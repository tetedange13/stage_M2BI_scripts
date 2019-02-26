#!/bin/bash

# This script aims to catch the sequences of all 16S copies of a given 
# organism from the Zymo mock community.
# It uses the 16S seq (FASTA file provided by Zymo) and BLAST it (locally) 
# against the genome of the given organism (FASTA file also provided by Zymo)   


to_genomes=Zymo_genomes-ZR160406
to_16S=Zymo_16S-ZR160406
to_BLASTdbs=BLASTdbs

for genome_file in $to_genomes/*fasta
do
    # CREATE BLAST DB:
    echo ""
    full_sp_name=$(basename $genome_file .fasta)
    echo BUILDING BLAST DB.. "for" $full_sp_name
    genus=$(echo $full_sp_name | cut -d '_' -f 1)
    current_BLASTdb=$to_BLASTdbs/$genus"_BLASTdb"
    
    $apps/blast-2.8.1+-src/ReleaseMT/bin/makeblastdb -in $genome_file \
    -dbtype nucl -parse_seqids -out $current_BLASTdb
    
    
    # RUN BLAST:
    echo ""
    rRNA_filename=$(basename $full_sp_name .fasta)"_16S.fasta"
    echo RUNNING BLAST.. "ON:" $rRNA_filename
    echo "TAILLE 16S:" $(cat $to_16S/$rRNA_filename | grep -v ">" | wc -c)
    out_BLAST=$(cat $to_16S/$rRNA_filename | \
              $apps/blast-2.8.1+-src/ReleaseMT/bin/blastn -task=blastn \
              -num_threads=10 -outfmt 6 -db $current_BLASTdb \
              -max_target_seqs 100 -word_size=15)
              
    
    # EXTRACT RANGES FROM BLAST OUTPUT:
    #ranges=$(echo "$out_BLAST" | awk '{ if ($4 > 1000) printf $9 "-" $10 " " }' | sort -k1,1)
    ranges=$(echo "$out_BLAST" | awk '{ if ($4 > 1000) print $2,$4,$9"-"$10 }' | sort -k1,1)
    #echo "$ranges"
    
    
    # EXTRACT SEQUENCES OF FOUND 16S:
    compt=0
    OLDIFS=$IFS
    IFS=$'\n'
    
    for line in $out_BLAST
    #echo "$out_BLAST" | while read line # Pour parcourir un fichier
    do
        echo $line | awk '{print $2,$4,$9"-"$10 }'
        
        length=$(echo $line | awk '{printf $4}')
        if [ "$length" -gt "1000" ]; then
            #echo $line | awk '{print $2,$4,$9"-"$10}'
            compt=$((compt+1))
        fi
        #seq_16S=$(echo LM1 $range | \
        #          $apps/blast-2.8.1+-src/ReleaseMT/bin/blastdbcmd \
        #          -db Listeria_BLASTdb -entry_batch - -out $genus"_"$range.fa)
    done
    echo "NB 16S FOUND:" $compt
    IFS=$OLDIFS
done

