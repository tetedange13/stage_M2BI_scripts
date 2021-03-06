#!/usr/bin/env python3


import sys, os, gzip
import os.path as osp
import time as t
import multiprocessing as mp
import itertools as ittls
from functools import partial
import src.parallelized as pll
from Bio import SeqIO

eval = pll.eval
taxfoo = pll.eval.taxfoo
pd = pll.eval.pd


global problematic_taxids # taxid that need direct remote on taxo DB
# More often entries that have been renamed and have an other taxid
problematic_taxids = {"65071", "82939", "162153", "307514", "1940813",
                      "1940819"}


def seqid2taxid_from_taxmap(to_taxmap_file, to_headers_list_file):
    """
    """
    dict_acc2taxid = {}
    with open(to_taxmap_file, 'r') as taxmap:
        for line in taxmap:
            pacc, _, _, _, _, taxid = line.rstrip('\n').split('\t')
            taxid = taxid.strip()
            
            if pacc not in dict_acc2taxid.keys():
                dict_acc2taxid[pacc] = taxid
            else:
                assert(dict_acc2taxid[pacc] == taxid)

    
    done_headers = set()
    with open(to_headers_list_file, 'r') as headers_file, \
         open("seqid2taxid", 'w') as seqid2taxid_file:
        for header in headers_file:
            pacc = header.rstrip('\n').split('.')[0] # We get only the 1st part
            if header not in done_headers:
                done_headers.add(header)
                seqid2taxid_file.write(header.rstrip('\n') + '\t' + 
                                       dict_acc2taxid[pacc] + '\n')
            else:
                print("PROBLEM")
                print(header)
                sys.exit(2)


def detect_problems(to_seqid2taxid, to_seqid2acc, to_target_list, dirOut="./"):
    """
    to_target_list = list of seqid that have been mapped
    Detect taxids that are 'unknown' within our current NCBI taxo
    """
    import src.ncbi_taxdump_utils as taxo_utils
    print("Loading taxonomic module...")
    taxfoo = taxo_utils.NCBI_TaxonomyFoo()
    nodes_path = to_dbs + "nt_db/taxo_18feb19/nodes.dmp"
    names_path = to_dbs + "nt_db/taxo_18feb19/names.dmp"
    taxfoo.load_nodes_dmp(nodes_path)
    taxfoo.load_names_dmp(names_path)
    print("Module loaded !")

    dict_seqid2taxid = {}
    with open(to_seqid2taxid, 'r') as seqid2taxid_file:#, \
        for line in seqid2taxid_file:
            seqid, taxid = line.rstrip('\n').split('\t')
            assert(seqid not in dict_seqid2taxid.keys())
            dict_seqid2taxid[seqid] = taxid

    dict_seqid2acc = {}
    with open(to_seqid2acc) as seqid2acc_file:
        for line in seqid2acc_file:
            seqid, pacc = line.rstrip('\n').split()
            dict_seqid2acc[seqid] = pacc

    set_problems = set()
    done_taxids = set()
    set_uniden = set()
    with open(to_target_list, 'r') as target_list:
        for line in target_list:
            ref_name = line.rstrip('\n')
            current_taxid = dict_seqid2taxid[ref_name]

            name_taxid = taxfoo.get_taxid_name(int(current_taxid))
            if name_taxid == "unidentified":
                set_uniden.add(current_taxid)


            if (not name_taxid and current_taxid or 
                current_taxid in problematic_taxids):
                if current_taxid not in done_taxids:
                    set_problems.add(current_taxid + '\t' + 
                                     dict_seqid2acc[ref_name] + '\n')
                    done_taxids.add(current_taxid)

    print()
    if set_problems:
        print("PROBLEMS FOUND..")
        with open(osp.join(dirOut, "problems.txt"), 'w') as pbs_file:
            [pbs_file.write(problem) for problem in set_problems]
    else:
        print("NO PROBLEMS FOUND !")
    print("Taxids associated to 'unidentified':", set_uniden)
    print()


def remote_taxid_search(to_problems_file, dirOut="./"):
    """
    Search taxid using remmote Biopython, used in cases when local search failed 
    """
    import src.remote as remote
    print("REMOTE TAXID SEARCH...")

    dict_wrong2good = {} # Conversion between wrong and good taxid
    to_wrong2good_file = osp.join(dirOut, "wrong2good_taxids")

    if osp.isfile(to_wrong2good_file):
        with open(to_wrong2good_file, 'r') as wrong2good_file:
            for line in wrong2good_file:
                wrong_taxid, good_taxid = line.rstrip('\n').split('\t')
                assert(wrong_taxid != good_taxid)
                dict_wrong2good[wrong_taxid] = good_taxid 

    possible_issues = ("HTTP400", "TAXO_NOT_FOUND", "SEVERAL_TAXIDS", 
                       "GI_NOT_FOUND")
    with open(to_problems_file, 'r') as pbs_file:
        for line in pbs_file:
            wrong_taxid, pacc = line.rstrip('\n').split('\t')

            if wrong_taxid not in dict_wrong2good.keys():
                good_taxid, found_organism = remote.taxid_from_gb_entry(pacc)

                if (wrong_taxid == good_taxid or 
                    wrong_taxid in problematic_taxids):
                    print("PB WITH TAXID FOUND IN GB ENTRY !")
                    print("   --> Requesting it on the NCBI Taxonomy db")
                    good_taxid = remote.query_taxid(found_organism)

                if good_taxid in possible_issues:
                    with open(osp.join(dirOut, "remote_fails"), 'a') as fails:
                        fails.write(pacc + '\t' + good_taxid + '\n')
                else:
                    dict_wrong2good[wrong_taxid] = good_taxid
                    with open(to_wrong2good_file, 'a') as wrong2good_file:
                        wrong2good_file.write(wrong_taxid + '\t' + 
                                      dict_wrong2good[wrong_taxid] + '\n')
    print("FINISHED !\n")


def local_taxid_search(gz_db_name, to_problems_file, dirOut="./", 
                       nb_threads=10):
    """
    Parallel local taxid mapping:
    """
    to_dbs = ("/mnt/72fc12ed-f59b-4e3a-8bc4-8dcd474ba56f/" +
              "metage_ONT_2019/nt_db/taxo_18feb19/")
    to_gz_acc2taxid = to_dbs + "nucl_" + gz_db_name + ".accession2taxid.gz"
    TAXID_START_TIME = t.time()
    dict_acc2wrong = {}

    with open(to_problems_file, 'r') as pbs_file:
        for line in pbs_file:
            wrong_taxid, pacc = line.rstrip('\n').split('\t')
            if wrong_taxid not in dict_acc2wrong.values():
                dict_acc2wrong[pacc] = wrong_taxid

    size_chunks = 1000000
    results = []
    my_pool = mp.Pool(nb_threads)
    with gzip.open(to_gz_acc2taxid, 'rt') as acc2taxid_gz:
        acc2taxid_gz.readline() # Skip header
        worker = partial(pll.taxid_mapping, 
                         set_accession=set(dict_acc2wrong.keys()))
        print("\nMapping taxids against accession numbers...")
        print("DB USED:", gz_db_name)

        compt_cut_done = 1
        while True:
            print("cut", compt_cut_done, ":", nb_threads*size_chunks)
            compt_cut_done += 1
            groups = [list(ittls.islice(acc2taxid_gz, size_chunks)) for i in range(nb_threads)]

            # As soon as 1st element not an empty list
            # (if the 1st one IS empty, the following will be too)
            if groups[0]: 
                result = my_pool.imap_unordered(worker, groups)
                results.extend(result)
            else:
                print("TOUT CONSOMME")
                break
    my_pool.close()

    dict_acc2good = {} # To convert accession number to GOOD taxid
    to_wrong2good_file = osp.join(dirOut, "wrong2good_taxids")
    with open(to_wrong2good_file, 'a') as wrong2good_file:
        for res in results:
            if res: # If  NOT 'None'
                for found_pacc, found_taxid in res:
                    if found_taxid not in dict_acc2good.values():
                        old_taxid = dict_acc2wrong[found_pacc]
                        dict_acc2good[found_pacc] = found_taxid
                        wrong2good_file.write(old_taxid + '\t' + 
                                              found_taxid + '\n')

    nb_found = len(dict_acc2good.keys())
    nb_to_find = len(dict_acc2wrong.keys())
    if nb_found != 0:
        to_need_remote = osp.join(dirOut, gz_db_name +"_need_remote")
        with open(to_need_remote, 'w') as need_remote_file:
            found_acc_numbers = dict_acc2good.keys()
            for acc_to_find in dict_acc2wrong:
                if acc_to_find not in found_acc_numbers:
                    need_remote_file.write(dict_acc2wrong[acc_to_find] + '\t' + 
                                           acc_to_find + '\n')

    print("PARALLEL SEARCH TIME:", t.time() - TAXID_START_TIME)
    print("TO_FIND:", nb_to_find)
    print("--> FOUND:", nb_found)
    print(nb_to_find - nb_found, "pacc couldn't be found locally!")
    print()


def correct_seqid2taxid(to_old_seqid2taxid, to_wrong2good_taxids, dirOut='./'):
    """
    Correct seqid2taxid file, by replacing taxids that have been found to be
    wrong (taxid not found in the nt NCBI taxonomy dump files), by correct 
    taxids (possibly found remotely or locally)
    """
    print("CORRECTING TAXIDS...")
    dict_wrong2good = {}
    with open(to_wrong2good_taxids, 'r') as wrong2good_file:
        for line in wrong2good_file:
            wrong_taxid, good_taxid = line.rstrip('\n').split('\t')
            assert(wrong_taxid not in dict_wrong2good.keys())
            dict_wrong2good[wrong_taxid] = good_taxid

    all_wrong_taxids = dict_wrong2good.keys()
    compt_corr = []
    with open(to_old_seqid2taxid, 'r') as old_seqid2taxid, \
         open(osp.join(dirOut, "seqid2taxid"), 'w') as new_seqid2taxid:
        for line in old_seqid2taxid:
            seqid, old_taxid = line.rstrip('\n').split('\t')
            if old_taxid in all_wrong_taxids:
                compt_corr.append(old_taxid)
                new_seqid2taxid.write(seqid + '\t' + 
                                       dict_wrong2good[old_taxid] + '\n')
            else:
                new_seqid2taxid.write(line)

    print("FINISHED!", len(compt_corr), "CORRECTIONS DONE")
    print("(representing", len(set(compt_corr)), "unique taxids)")
    print()


def write_complete_lineage(to_seqid2taxid):
    """
    Write the file containing the list of all the taxids (complete lineage) of
    a given seqid2taxid file
    """
    with open(to_seqid2taxid, 'r') as seqid2taxid_file:
        set_taxids_db = {line.rstrip('\n').split('\t')[1] 
                         for line in seqid2taxid_file}

    set_complete_lineage = set()
    for taxid_from_db in set_taxids_db:
        lineage = taxfoo.get_lineage_as_taxids(taxid_from_db)
        for taxid in lineage:
            set_complete_lineage.add(str(taxid))
        del taxid
    del taxid_from_db

    # We need to add the 'root' (taxid=1):
    set_complete_lineage.add('1')

    print("LONG COMPLETE LINEAGE:", len(set_complete_lineage))

    with open('taxids_complete_lineage', 'w') as complete_lineage_file:
        [complete_lineage_file.write(taxid + '\n') 
         for taxid in set_complete_lineage]


def parse_and_rewrite_nodes(filename, to_complete_lineage):
    """"
    Parse the NCBI nodes_dmp file. and rewrite it, to keep only species that 
    within a certain database 
    """
    with open(to_complete_lineage, 'r') as complete_lineage_file:
        set_taxids_to_keep = {line.rstrip('\n') 
                              for line in complete_lineage_file}

    with open(filename, 'r') as fp, \
         open("nodes.dmp", 'w') as new_nodes_file:
        for line in fp:

            node_id = line.split('\t|\t')[0]
            if node_id in set_taxids_to_keep:
                new_nodes_file.write(line)
        del line


def parse_and_rewrite_names(filename, to_complete_lineage):
    """
    Parse the NCBI names_dmp file. and rewrite it, to keep only species that 
    within a certain database 
    """
    with open(to_complete_lineage, 'r') as complete_lineage_file:
        set_taxids_to_keep = {line.rstrip('\n') for line in complete_lineage_file}

    with open(filename, 'r') as fp, \
         open("names.dmp", 'w') as new_names_file:

        # Write header of the 'names.dmp' file:
        new_names_file.write(fp.readline())
         
        for line in fp:
            line = line.rstrip('\t|\n')
            taxid, name, uniqname, name_class = line.split('\t|\t')

            if (name_class == 'scientific name' and 
                taxid in set_taxids_to_keep):
                new_names_file.write(line + '\t|\n')
        del line


def write_metadat_file(to_seqid2taxid, base_used):
    """
    Given a complete lineage as the set of taxids, writes (temporary) the 
    complete lineage for this given set of taxids and write associating each 
    taxid from this 'complete_lineage' with its SILVA-formatted taxonomy
    """
    write_complete_lineage(to_seqid2taxid)

    with open('taxids_complete_lineage', 'r') as complete_lineage, \
         open(base_used + '_tax_metadat.tsv', 'w') as metadat_file:
        metadat_file.write('unmapped' + '\t' + ';'.join(['no'] * 7) + 
                           '\n')
        
        for line in complete_lineage:
            taxid_grp = line.rstrip('\n')
            metadat_file.write(taxid_grp + '\t' + 
                               eval.taxo_from_taxid(taxid_grp) + '\n')
        del line

    os.remove('taxids_complete_lineage')


def stats_base(db_to_stat):
    """
    Produce statistics of each databases, concerning the number of different
    taxonomic ranks, present sequences of our 8 expected bacteria etc.
    """
    print("Proceeding...")
    
    to_seqid2taxid = ("/mnt/72fc12ed-f59b-4e3a-8bc4-8dcd474ba56f/" +
                      "metage_ONT_2019/Centri_idxes/" + db_to_stat.lower() +
                      '/seqid2taxid')

    dict_name2taxid = {'Listeria monocytogenes':1639, 
                       'Bacillus subtilis':1423, 
                       'Staphylococcus aureus':1280, 
                       'Escherichia coli':562, 
                       'Lactobacillus fermentum':1613, 
                       'Enterococcus faecalis':1351,
                       'Pseudomonas aeruginosa':287, 
                       'Salmonella enterica':28901,
    # dict_euk = { 
                'Saccharomyces cerevisiae':4932,
                'Cryptococcus neoformans':5207}

    taxids_zymo_sp = list(dict_name2taxid.values())
    taxids_zymo_gen = [taxfoo.get_taxid_parent(taxid) for taxid in taxids_zymo_sp]
    taxids_zymo_gen += [286, 5206, 1386] # Add 'Pseudomonas', 'Bacillus' and 'Cryptococcus'
    
    with open(to_seqid2taxid, 'r') as seqid2taxid_file:
        list_taxids = [line.rstrip('\n').split('\t')[1] 
                       for line in seqid2taxid_file]

    # print(sum(map(pll.is_trash, list_taxids)));sys.exit()
    
    dict_tmp = {}
    dict_bacil = {'Bacillus intestinalis':0, 'Bacillus subtilis':0}
    dict_counts_sp = {a_key:0 for a_key in dict_name2taxid}
    dict_counts_gen = {a_key.split()[0]:0 for a_key in dict_name2taxid}

    for taxid in list_taxids:
        dict_tmp[taxid] = taxfoo.get_taxid_rank(int(taxid))
        if taxfoo.is_strain(int(taxid)): # Strain => rk='no rank' && rk_parent='species'
            dict_tmp[taxid] = 'strain'
        lineage = taxfoo.get_dict_lineage_as_taxids(taxid)
        possible_tax_lvls = lineage.keys()

        if 'species' in possible_tax_lvls:
            sp_taxid = lineage['species']
            sp_name = taxfoo.get_taxid_name(sp_taxid)
            if sp_taxid in taxids_zymo_sp:
                dict_counts_sp[sp_name] += 1
            if sp_name in dict_bacil.keys():
                dict_bacil[sp_name] += 1

        if 'genus' in possible_tax_lvls:
            gen_taxid = lineage['genus']
            if gen_taxid in taxids_zymo_gen:
                dict_counts_gen[taxfoo.get_taxid_name(gen_taxid)] += 1
    del taxid

    print(dict_bacil)
        
    print("STATS OF {} DB:".format(db_to_stat.upper()))
    print(pd.Series(dict_tmp).value_counts())
    print()
    print(pd.Series(dict_counts_gen, 
          name='GENUS lvl').sort_values(ascending=False))
    print()
    print(pd.Series(dict_counts_sp,
          name='SPECIES lvl').sort_values(ascending=False))
    print()


def transform_ONT_CSV(to_initial_csv):
    """
    Transform the CSV produced by Epi2me to make it usable by my scripts
    """
    print("Transforming EPI2ME CSV..")
    base_filename = osp.basename(to_initial_csv).lower()

    if "epi2me" in base_filename:
        print("Detected EPI2ME tool !")
        name_col_readID = 'read_id'
        name_col_taxid = 'taxid'
    elif "wimp" in base_filename:
        print("Detected WIMP tool !")
        name_col_readID = 'readid'
        name_col_taxid = 'taxID'
    else:
        print("ERROR ! Filename must contain either 'EPI2ME' or 'WIMP'")
        sys.exit()

    initial_csv = pd.read_csv(to_initial_csv, header=0, sep=',', 
                             # usecols=['exit_status', 'taxid', 'accuracy', 
                             #          'lca'])
                             usecols=[name_col_readID, 'exit_status', 
                                      name_col_taxid])

    if all(initial_csv[name_col_readID].isnull()): # Actually NO readID...
        print("NO readID actually found --> adding a pseudo-readID..")
        nb_rows = len(initial_csv.index)
        initial_csv[name_col_readID] = [foo[0]+str(foo[1]) 
                                    for foo 
                                    in zip(['read_']*nb_rows, range(nb_rows))]
    else:
        # If 'False' here, will need more code to group entries...
        assert(len(initial_csv[name_col_readID]) == len(initial_csv[name_col_readID].unique()))
        
    initial_csv['lineage'] = initial_csv[name_col_taxid].apply(
                                            lambda val: ';' + str(val) + ';')


    def exitStatus_to_typeAlign(val):
        """
        /!\ This func assert that every status different from 
        'Classif successful' can be considered as 'unmapped' 
        """
        if val in ('Classification successful', 'Classified'):
        # 'Classified' for WIMP and 'Classification successful' for EPI2ME
            return 'normal'
        return 'unmapped'

    # Before applying this func, be sure not another category
    # possible_exit_status = ['Classification below QC threshold', 
    #                         'Classification successful', 'No classification']
    print("All different 'exit_status' in this file:")
    print(initial_csv.exit_status.value_counts())
    # assert(sorted(initial_csv.exit_status.unique()) == possible_exit_status)

    initial_csv['type_align'] = initial_csv.exit_status.apply(
                                                    exitStatus_to_typeAlign)
    not_unmapped = initial_csv.type_align != 'unmapped'
    initial_csv['nb_trashes'] = initial_csv[not_unmapped][name_col_taxid].apply(
                                        lambda val: int(pll.is_trash(val)))

    initial_csv.drop(columns=[name_col_taxid, 'exit_status'], inplace=True)
    print("Done ! Head of the written new CSV:")
    print(initial_csv.head())
    inbase = osp.splitext(osp.basename(to_initial_csv))[0]
    initial_csv.set_index([name_col_readID]).to_csv(inbase+'.csv', 
                                                    index_label=False)
    print("FINISHED ! Wrote: {}.csv".format(inbase))
    print()


def extract_reference_seq(to_extractable, to_zymo_SEGO):
    """
    Given a list of readID with their assignation and the positions of the 
    alignment within the target sequence (genome), produce a (pseudo-)fasta 
    file, containing the sequence of the aligned region within the reference
    genome (so no errors)
    These sequences will be then given to our custom version NanoSim as a 
    reference fasta, in order to simulate reads  
    """

    with open(to_extractable, 'r') as extractable_file:
        header = extractable_file.readline()

        dict_felix = {}
        for line in extractable_file:
            readID, ref_name, start, end = line.rstrip('\n').split(',')

            if ref_name not in dict_felix.keys():
                dict_felix[ref_name] = [(readID, int(start), int(end))]
            else:
                dict_felix[ref_name].append((readID, int(start), int(end)))
        del line

    records = SeqIO.parse(to_zymo_SEGO, "fasta")
    infile_base = osp.splitext(osp.basename(to_extractable))[0]
    with open(infile_base + '.fa', 'w') as out_fasta:
        for entry in records:
            a_ref_name = entry.id

            if a_ref_name in dict_felix.keys():
                for read_info in dict_felix[a_ref_name]:
                    readID, start, end = read_info
                    str_to_write = '>' + readID + '_' + a_ref_name # header fa
                    out_fasta.write(str_to_write + '\n' + 
                                    str(entry.seq[start:end]) + '\n')
                del read_info
        del entry

    print("Wrote:", infile_base + '.fa')



def detect_RRN_litige(to_acc2taxid, to_species_annot):
    """
    Detect in the RRN db, cases where the annotated species name is not 
    coherent with the given GI (and its associated taxid)
    """
    with open(to_acc2taxid, 'r') as acc2taxid_file:
        dict_acc2taxid = {line.rstrip('\n').split('\t')[0]:line.rstrip('\n').split('\t')[1] 
                          for line in acc2taxid_file}

    with open(to_species_annot, 'r') as sp_annot:
        dict_litiges = {}
        for line in sp_annot:
            annot_name, op_id, gi, _ = line.rstrip('\n').split('\t')
            corresponding_taxid = dict_acc2taxid[gi]
            lineage = taxfoo.get_lineage_as_dict(corresponding_taxid)

            if 'species' in lineage.keys():
                sp_name = lineage['species']
                splitted_sp_name = sp_name.split()
                if len(splitted_sp_name) == 2: 
                    if '_'.join(splitted_sp_name) != annot_name:
                        print(gi, ' | ', sp_name, 'VS', annot_name)
                        dict_litiges[gi] = [sp_name, annot_name]
                else:
                    # if '_'.join(splitted_sp_name) != annot_name:
                    print(gi, ' | ', sp_name, 'VS', annot_name)

            else:
                print("NO SPECIES")

    print("TOT NB LITIGES (len==2):", len(dict_litiges.keys()))


def tiny_func_for_STAMP(to_biom):
    """
    Format the result biom convert -i a_biom.biom -o to_tsv_biom --to-tsv, to
    transform it into a spf profile file 
    """
    to_tsv_biom = "felix.tsv"
    tmp_csv = pd.read_csv(to_tsv_biom, sep='\t', header=1)
    tmp_csv.columns = ['taxid'] + list(tmp_csv.columns)[1:]
    not_no_majo_found = tmp_csv.taxid != 'no_majo_found'
    tmp_csv['sp_name'] = tmp_csv[not_no_majo_found].taxid.astype('int').apply(
                                                        taxfoo.get_taxid_name)
    tmp_csv.loc[-not_no_majo_found, 'sp_name'] = 'no_majo_found'
    tmp_csv.set_index('sp_name', inplace=True)
    tmp_csv.drop('taxid', axis='columns', inplace=True)

    out_spf = osp.splitext(to_biom)[0] + '.spf'
    tmp_csv.to_csv(out_spf, sep='\t')


def tiny_func_for_STAMP_2(to_biom_txt):
    """
    Taking a tsv format OTU table output from 'summarize_taxa.py' script 
    (can be a merged OTUs, with mmore than 2 different samples), merge all
    entries different from the eight expected into a 'misassigned' category
    and produce a proper spf file 
    """ 
    df_prok_ref = eval.generate_df_zymo()
    set_taxids = df_prok_ref.index
    my_taxo_csv = pd.read_csv(to_biom_txt, sep='\t', skiprows=[0])
    assert(len(my_taxo_csv.index) >= 8)

    list_taxo = [eval.taxo_from_taxid(taxid) for taxid in set_taxids]

    samples_names = my_taxo_csv.columns[1:]
    nb_samples = len(samples_names)
    print('Detected {} different samples'. format(nb_samples))

    misassigned = pd.Series([0]*nb_samples, index=samples_names)
    future_spf = pd.DataFrame(columns=samples_names)

    i = 0
    for idx, row in my_taxo_csv.iterrows():
        taxo = row["#OTU ID"]

        if taxo not in list_taxo:
            misassigned = misassigned.add(row[1:])
        else:
            splitted_taxo = taxo.split(';')
            sp_name = splitted_taxo[-1]
            if len(splitted_taxo) == 7:
                sp_name = splitted_taxo[-2] + ' ' + splitted_taxo[-1]

            future_spf.loc[sp_name, :] = row
            i += 1
    del idx, row

    future_spf.loc['misassigned', :] = misassigned
    assert(len(future_spf.index) == 9)
    to_out_spf = osp.splitext(to_biom_txt)[0] + '.spf'
    future_spf.index.name = 'sp_name'
    print("Wrote:", to_out_spf) ; print()
    future_spf.to_csv(to_out_spf, sep='\t', index=True)


def correct_spf(to_spf):
    """
    Another func trying to produce a properly formatted input file for STAMP...
    """
    print("Input spf file:", to_spf)
    test = pd.read_csv(to_spf, sep='\t')

    for idx, row in test.iterrows():
        are_other = row.apply(lambda val_from_list: val_from_list=='Other')

        if test.loc[idx, 'Level_1'] == 'None':
            print("Removing Unclassified entry !")
            test.drop(idx, inplace=True)

        if sum(are_other) > 0:
            print(test.loc[idx, :])
            test.loc[idx, are_other] = 'Unclassified'
            print(idx+2, sum(are_other), list(test.loc[idx, test.columns[-5:-1]]))

    test.loc[:, 'Level_7'] = test.Level_6 + ' ' + test.Level_7

    print('New spf file:', to_spf.rstrip('0'))
    test.to_csv(to_spf.rstrip('0'), sep='\t', index=False)


def parse_gbff(to_gbff):
    """
    Parse gbff file (NCBI "16S" db), in order to get all couples (acc_nb;taxid)
    """
    print("Parsing gbff file to get taxid associated with each acc_nb..")
    import re

    # We get all couples using 2 regex:
    reg_acc_nb = re.compile("(?:VERSION\s+)(.*)")
    reg_taxid = re.compile("(?:.*db_xref.*taxon:)([0-9]+)")

    with open(to_gbff, 'r') as gbff_file:
        to_search_in = ''
        list_couples = []

        for line in gbff_file:
            if line != "//\n": # Delimiter between each 
                to_search_in += line
            else:
                matches_acc_nb = reg_acc_nb.findall(to_search_in)
                matches_taxid = reg_taxid.findall(to_search_in)
                assert(len(matches_acc_nb) == 1)
                assert(len(matches_taxid) == 1)

                list_couples.append((matches_acc_nb[0], matches_taxid[0]))

                # We reset the str to search into:
                to_search_in = ''
    
    # We write the corresponding seqid2taxid as output: 
    with open("seqid2taxid", 'w') as seqid2taxid_file:
        for acc_nb, taxid in list_couples:
            seqid2taxid_file.write(acc_nb + '\t' + taxid + '\n')
        del acc_nb, taxid

    print("FINISHED! Found {} couples".format(len(list_couples)))
    print()



# MAIN:
if __name__ == "__main__":
    to_dbs = ("/mnt/72fc12ed-f59b-4e3a-8bc4-8dcd474ba56f/" +
              "metage_ONT_2019/")
    to_dbs_SILVA = to_dbs + "SILVA_refNR132_28feb19/"
    to_dbs_RRN = to_dbs + "rrn_8feb19/"
    to_dbs_nt = to_dbs + "nt_db/taxo_18feb19/"

    # seqid2taxid_from_taxmap(to_dbs_SILVA + "taxmap_embl_ssu_ref_nr99_132.txt", 
    #                        to_dbs_SILVA + "headers.txt")

    # test(to_dbs_SILVA + "wrong2good_taxids", "seqid2taxid")

    # detect_problems("seqid2taxid", to_dbs_SILVA + "seqid2acc", 
    #                 to_dbs_SILVA + "headers.txt")
    #detect_problems("seqid2taxid", "seqid2acc", "headers.txt")

    # local_taxid_search("wgs", "problems.txt")
    # local_taxid_search("gb", "wgs_need_remote")
    # local_taxid_search("gss", "gb_need_remote")

    # remote_taxid_search("problems.txt")
    
    # correct_seqid2taxid("old_seqid2taxid", "wrong2good_taxids")

    # write_complete_lineage("seqid2taxid")
    write_metadat_file(to_dbs + 'Centri_idxes/ncbi16S/seqid2taxid', 'toZymo')
    # parse_and_rewrite_names(to_dbs_nt + "names.dmp", 
    #                         "taxids_complete_lineage")
    # parse_and_rewrite_nodes(to_dbs_nt + "nodes.dmp", 
    #                         "taxids_complete_lineage")

    # stats_base('RRN')
    # detect_RRN_litige(to_dbs_RRN+'acc2taxid', to_dbs_RRN+'species_annotation')
    # transform_ONT_CSV("../2_WIMP_cDNA_run1_trimmed_toNCBIbact.sam")
    # transform_ONT_CSV("../EPI2ME_16Skit_run1_500K_toNCBIbact.sam")
    # extract_reference_seq("extractable_16SkitRun1Zymo.csv", "../../zymo_SEGO.fa")

    # tiny_func_for_STAMP_2("OTUs_maps_n_tables/toZymo_L7.txt")
    # correct_spf("OTUs_maps_n_tables/toRrn.spf0")

    if False:
    # if len(sys.argv) > 1:
        # tiny_func_for_STAMP(sys.argv[1])
        tiny_func_for_STAMP_2(sys.argv[1])

    # parse_gbff(to_dbs + "NCBI_16S_18-06-2019/bacteria.16SrRNA.gbff")
