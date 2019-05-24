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
`summarize_taxa.py -i my_biom_w_metadat.biom -o . --level 7 --absolute_abundance --suppress_biom_table_output`
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
        if taxonomic_cutoff == 'species':
            taxa_name = row['genus'] + " " + row[taxonomic_cutoff]
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
    


# MAIN
if __name__ == '__main__':
    ARGS = docopt(__doc__, version='0.1')
    to_infile, infile_base, _, _, ext_infile = check.infile(ARGS["--inOTUs"], 
                                                        ['csv', 'tsv', 'txt'])
    
    print()
    print("Loading taxonomic Python module...")
    import src.taxo_eval as eval
    pd = eval.pd
    print("Taxonomic Python module loaded !\n")
    taxo_cutoff = check.acceptable_str(ARGS["--taxoCut"], 
                                       eval.want_taxo + ['kingdom'])

    print("Processing:", infile_base)

    # With the 'python engine', it is possible to use regex as separator
    # ( here, it says sep = ';' OR '\t' ):
    my_taxo_csv = pd.read_csv(to_infile, sep=';|\t', skiprows=[0, 1], 
                         engine='python', names=eval.want_taxo+['abs_count'])
    print(my_taxo_csv.nlargest(10, 'abs_count')[['genus', 'species', 'abs_count']])

    counts_sp_obs = gether_counts(my_taxo_csv, taxo_cutoff, 'tot_sp_obs')
    # print(counts_sp_obs)

    tot_nb_mapped = sum(counts_sp_obs)
    print()
    print("TOT NB OF MAPPED READS:", tot_nb_mapped)

    # Draw pie chart of abundances (simple counts):
    # plot_pie_chart(counts_sp_obs, False, "MON BEAU CAMEMBERT")
    # plt.show()

    df_prok_ref = eval.generate_df_zymo()
    df_prok_ref = df_prok_ref.assign(abs_count=df_prok_ref.rel_count * 
                                               tot_nb_mapped / 100)
    # print(df_prok_ref.rel_count);sys.exit()
    df_counts = pd.DataFrame(gether_counts(df_prok_ref, taxo_cutoff, 
                             'expect_counts')) 
    set_expected = df_counts.index
    df_counts.loc['0-misassigned', 'expect_counts'] = 0 # Add misassigned
    df_counts.expect_counts = pd.to_numeric(df_counts.expect_counts,
                                            downcast='integer')
    
    # Need to add some (virtual) reads to the "expected counts":
    tot_expect = sum(df_counts.expect_counts)
    if tot_expect != tot_nb_mapped:
        gap = tot_nb_mapped - tot_expect 
        assert(gap < 9 and pd.np.sign(gap) == 1.)
        for i in range(gap): df_counts.expect_counts[i] += 1

    # Gether all 'misassigned' together:
    dict_tmp = {'0-misassigned':0}
    for taxa_name, count in counts_sp_obs.items():
        if taxa_name not in set_expected:
            dict_tmp['0-misassigned'] += count
        else:
            dict_tmp[taxa_name] = count
    del taxa_name, count

    df_counts = df_counts.assign(obs_counts=[dict_tmp[tax_name] 
                                             for tax_name in df_counts.index])
    df_counts = df_counts.assign(difference=df_counts.obs_counts - 
                                            df_counts.expect_counts)

    log_modulus = lambda diff: (pd.np.sign(diff) * 
                                pd.np.log10(1 + pd.np.abs(diff)))
    df_counts = df_counts.assign(log_modulus=lambda x: log_modulus(x.difference),
                                 L1dist=lambda x: pd.np.abs(x.difference))
    df_counts.sort_index(inplace=True)
    
    print(df_counts)
    print()
    print("L1-distance =", sum(df_counts.L1dist))
    print()

    log_modulus_series = df_counts.log_modulus
    print("Values of log-modulus:")
    if taxo_cutoff == 'species':
        print(';'.join(map(lambda a_str: a_str.split()[-1], 
                           log_modulus_series.index)))
    else:
        print(';'.join(log_modulus_series.index))
    print(';'.join(str(round(val, 3)) for val in log_modulus_series.values))
    print()
    
    # Relative abundances:
    print("Relative abundances (freq from raw counts):")
    print(';'.join(str(round(val, 3)) for val in df_counts.obs_counts/sum(df_counts.obs_counts)))
    # print(';'.join(str(round(val, 3)) for val in df_counts.expect_counts/sum(df_counts.expect_counts)))


    plot_stacked_bar = False
    if plot_stacked_bar:
        tmp_series = df_counts.obs_counts.drop('0-misassigned').sort_index(ascending=False)/sum(df_counts.obs_counts)
        # print(tmp_series);sys.exit()
        tmp_df = pd.DataFrame([tmp_series.to_dict(), 
                               {felix:0 for i, felix in enumerate(tmp_series.index)}], 
                              index=[0,1])#.sort_values(by=1, axis='columns')
        # Invert to have the proper order on the stacked barplot:
        tmp_df = tmp_df.reindex(sorted(tmp_df.columns, reverse=True), axis=1)
        print(tmp_df)
        tmp_df.plot(kind='bar', stacked=True, width=0.1, legend=True, 
                    title='Cusco2018 16S_run2', yticks=[0, 0.25, 0.5, 0.75, 1])
        plt.show()

    sys.exit()

    print()
    plot_pie_chart(df_counts.obs_counts.sort_values(ascending=False), True, 
                   "Absolute abundances\n at the {} level".format(taxo_cutoff))
    plt.show()
