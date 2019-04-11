#!/usr/bin/env python3

"""
Mapping statistics computation

Usage:
  second_analysis.py (-i <inOTUs>) (-l <taxoCut>)
  
Options:
  -h --help                  help
  --version                  version of the script
  -i --inOTUs=input_otu      input OTU table
  -l --taxoCut=taxo_cutoff   cutoff for the taxonomic level

Remarks: 
1) The OTUs table given as input must be properly formatted (as it is
explained on qiime.org) and be at the SPECIES level
2) Ideally the file resulting from this cmd: 
`summarize_taxa.py -i my_biom_w_metadat.biom -o . --level 7 -a --suppress_biom_table_output`
"""

import os, sys
from docopt import docopt
import matplotlib.pyplot as plt
import src.check_args as check


def gether_counts(a_taxo_csv, taxonomic_cutoff, name_series):
    """
    """
    tmp_dict = {}
    for idx, row in a_taxo_csv.iterrows():
        taxa_name = row[taxonomic_cutoff]
        if taxa_name not in tmp_dict.keys():
            tmp_dict[taxa_name] = int(row['abs_count'])
        else:
            tmp_dict[taxa_name] += int(row['abs_count'])
    del idx, row

    return pd.Series(tmp_dict, name=name_series).sort_values(ascending=False)


def plot_pie_chart(pdSeries_to_plot, print_percents, arg_title=""):
    fig, axis = plt.subplots(subplot_kw=dict(aspect="equal"))
    if print_percents:
        patches, _ , _= axis.pie(pdSeries_to_plot, startangle=90, 
                                 counterclock=False, 
                                 autopct=lambda pct: "{:.1f}%".format(pct))
    else:
        patches, _ = axis.pie(pdSeries_to_plot, startangle=90, 
                                 counterclock=False)
    
    lim_nb_items_leg, nb_items, y_val = 31, len(pdSeries_to_plot.index), 0
    # print(nb_items)
    if nb_items > lim_nb_items_leg:
        y_val = -0.022 * (nb_items - lim_nb_items_leg)
    axis.legend(patches, pdSeries_to_plot.index, loc="center left", 
                bbox_to_anchor=(0.95, y_val, 0.5, 1), fontsize='x-small')
    axis.set_ylabel('') # To remove auto legends sided to the chart
    fig.subplots_adjust(left=-0.4)
    axis.set_title(arg_title)


def calc_L1dist_logModulus(pdSeries_obs_abund, pdSeries_expect_abund):
    """
    Compute both the L1-distance and the vector of log-modulus values between 
    expected and observed abundances

    From: https://genomebiology.biomedcentral.com/track/pdf/10.1186/s13059-017-1299-7
    "While the log-modulus examines a fold-change, the L1 distance shows the 
    distance between relative abundance vectors by dataset"
    """
    log_modulus = lambda difference: (pd.np.sign(difference) * 
                                      pd.np.log10(1 + pd.np.abs(difference)))

    dict_L1dist,dict_log_modulus = {}, {}
    for sp_name in pdSeries_obs_abund.index:
        diff = pdSeries_obs_abund[sp_name] - pdSeries_expect_abund[sp_name]
        dict_L1dist[sp_name] = pd.np.abs(diff)
        dict_log_modulus[sp_name] = log_modulus(diff)

    # print(dict_L1dist)
    return (sum(dict_L1dist.values()), 
            pd.Series(dict_log_modulus).sort_values(ascending = False))



# MAIN
if __name__ == '__main__':
    ARGS = docopt(__doc__, version='0.1')
    to_infile, infile_base, _, _, ext_infile = check.infile(ARGS["--inOTUs"], 
                                                        ['csv', 'tsv', 'txt'])
    print("Loading taxonomic Python module...")
    import src.taxo_eval as eval
    pd = eval.pd
    print("Taxonomic Python module loaded !\n")
    taxo_cutoff = check.acceptable_str(ARGS["--taxoCut"], 
                                       eval.want_taxo + ['kingdom'])

    # With the 'python engine', it is possible to use regex as separator
    # ( here, it says sep = ';' OR '\t' ):
    my_taxo_csv = pd.read_csv(to_infile, sep=';|\t', skiprows=[0, 1], 
                         engine='python', names=eval.want_taxo+['abs_count'])
    counts_sp_obs = gether_counts(my_taxo_csv, taxo_cutoff, 'tot_sp_obs')

    # Draw pie chart of abundances (simple counts):
    plot_pie_chart(counts_sp_obs, False, "MON BEAU CAMEMBERT")
    # plt.show()

    df_prok_ref = eval.generate_df_zymo()
    df_prok_ref = df_prok_ref.assign(abs_count=df_prok_ref['rel_count']*sum(counts_sp_obs)/100)


    # Gether all FP into a 'misassigned' category:
    # /!\ 'no_majo_found' are ignored /!\
    # dict_species2count = {'misassigned':0}
    # for counted_item, count in counts_species.items():
    #     if counted_item != 'no_majo_found':
    #         if dict_species2res[counted_item] == 'FP':
    #             dict_species2count['misassigned'] += count
    #         else:
    #             dict_species2count[counted_item] = count
    # del counted_item

    # counts_only_zymo = pd.Series(dict_species2count).sort_values(ascending=False)
    # print(counts_only_zymo)


    df_counts = pd.DataFrame(gether_counts(df_prok_ref, taxo_cutoff, 
                             'expect_counts')) 
    
    set_expected = df_counts.index
    df_counts.loc['misassigned', 'expect_counts'] = 0 # Add misassigned
    df_counts.expect_counts = pd.to_numeric(df_counts.expect_counts,
                                              downcast='integer')

    # Gether all 'misassigned' together:
    dict_tmp = {'misassigned':0}
    for taxa_name, count in counts_sp_obs.items():
        if taxa_name not in set_expected:
            dict_tmp['misassigned'] += count
        else:
            dict_tmp[taxa_name] = count
    print(dict_tmp)

    df_counts = df_counts.assign(obs_counts=dict_tmp.values())
    print(df_counts)
    sys.exit()

    print()
    plot_pie_chart(counts_only_zymo, True, "MON BEAU CAMEMBERT")
    plt.show()

    if taxo_cutoff == 'species':
        
        # print(counts_only_zymo.drop(['misassigned']))

        not_misassigned = counts_only_zymo.drop(['misassigned'])
        tot_reads_zymo = sum(not_misassigned)
        foo_dict = {sp_name:expected_percents[sp_name]*tot_reads_zymo/100 
                    for sp_name in expected_percents.keys()}
        # print(pd.Series(test))
        a_tupl = calc_L1dist_logModulus(not_misassigned, 
                                        pd.Series(foo_dict))
        L1dist, vect_logModulus = a_tupl
        # print("L1-distance =", L1dist)
        print()
        # print(vect_logModulus)
        # print(counts_species[map(lambda val: val.split()[0] == 'Bacillus', counts_species.index)])
        # plt.show()
        # sys.exit()



