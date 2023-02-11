from pathlib import Path
import argparse

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.patches as mpatches

# matplotlib.use('AGG')#or PDF, SVG and PS


def parse_varscan():
    path = Path(__file__).parent.parent/'test'

    filename = 'S26.2499399.varscan.readcounts'
    with open(path/'readcounts'/filename, 'r') as rf:
        with open(path/'readcounts_csv'/filename, 'w') as wf:
            for line in rf:
                line_csv = line.strip().split('\t')
                print(line_csv)
                line_c1_c5 = line_csv[0:5]
                line_c1_c5 = ','.join(line_c1_c5)
                line_c6 = line_csv[5:]
                print(line_c6)
                line_c6 = ';'.join(line_c6)
                line_csv = line_c1_c5 + ',' + line_c6 + '\n'
                wf.write(line_csv)

if __name__ == '__main__':

    Description = "Detect cross contamination through allele frequency of samples "
    Epilog = """Example usage: python3 detect_contamination.py <file_in1> <file_in2> <file_out1> <file_out2>"""
    # python3 ./bin/detect_contamination.py \\
    # /media/data/fuhaoyi/gitDOWN/sars2_surveillance/test/output/5nextclade/test_nextclade.csv \\
    # /media/data/fuhaoyi/gitDOWN/sars2_surveillance/test/output/7readcount \\
    # /media/data/fuhaoyi/gitDOWN/sars2_surveillance/test/output/8contamination/vaf_table.csv \\
    # /media/data/fuhaoyi/gitDOWN/sars2_surveillance/test/output/8contamination/vaf_heatmap.png

    parser = argparse.ArgumentParser(description=Description, epilog=Epilog)
    parser.add_argument("file_in1", help="Path of xxx_nexclade.csv file")
    parser.add_argument("file_in2", help="Path of sample readcount folder")
    parser.add_argument("file_out1", help="Path of vaf_tabel.csv")
    parser.add_argument("file_out2", help="Path of vaf_heatmap.png")


    args = parser.parse_args()
    path1 = Path(args.file_in1)
    path2 = Path(args.file_in2)
    path_vaf_table = Path(args.file_out1)
    path_vaf_heatmap = Path(args.file_out2)

    # path1 = Path('/media/data/fuhaoyi/sequencing_data/20230130/output/5nextclade/20230130_nextclade.csv')
    # path2 = Path('/media/data/fuhaoyi/sequencing_data/20230130/output/7readcount')
    # path_vaf_table = Path('/media/data/fuhaoyi/sequencing_data/20230130/output/8contamination/vaf_table.csv')
    # path_vaf_heatmap = Path('/media/data/fuhaoyi/sequencing_data/20230130/output/8contamination/detect_comtam.png')




    nextclade = pd.read_csv(path1, sep=';', dtype={'seqName': 'string'})
    nextclade = nextclade.set_index('seqName')

    sample_list = nextclade.index.values.tolist()

    colors = mpl.colormaps['tab10'].resampled(10)        # matplotlib 3.6
    # colors = plt.cm.get_cmap('Set1', 10)              # matplotlib 3.5
    clade_list = nextclade['Nextclade_pango'].values.tolist()
    clade_color_map = set(clade_list)
    clade_color_map = {clade:colors(i % 10) for i, clade in enumerate(clade_color_map)}
    clade_color_list = [clade_color_map[i] for i in clade_list]

    substitutions_list = nextclade['substitutions'].dropna()
    substitutions_list = substitutions_list.values.tolist()
    substitutions_list = ','.join(substitutions_list)
    substitutions_list = substitutions_list.split(',')
    substitutions_list = set(substitutions_list)
    substitutions_list = sorted(substitutions_list, key=lambda x: int(''.join([i for i in x if i.isdigit()])))

    # print(sample_list)
    # print(substitutions_list)

    vaf_table = pd.DataFrame(data=0, index=sample_list, columns=substitutions_list)
    vaf_table.index.name = 'seqName'
    # print(vaf_tabel)

    for ix, row in vaf_table.iterrows():
        readcount = pd.read_csv(path2/f'{ix}'/f'{ix}.readcount.csv')
        readcount = readcount.set_index('position')
        # print(readcount)
        for mutation in vaf_table.columns:
            ref_base = mutation[0]
            position = int(mutation[1:-1])
            mut_base = mutation[-1]
            if position in readcount.index:
                vaf_dic = readcount.loc[position, 'vaf']
                vaf_dic = vaf_dic.split()
                vaf_dic = {i.split(':')[0]: float(i.split(':')[1]) for i in vaf_dic}
                vaf_table.loc[ix, mutation] = vaf_dic[mut_base]
            else:
                vaf_table.loc[ix, mutation] = 0

    vaf_table.to_csv(path_vaf_table)

    # plot heatmap
    # plt.style.use('fivethirtyeight')
    # g = sns.clustermap(vaf_table,
    #                    row_colors=clade_list_clor,
    #                    cmap="plasma")
    # patch_list = []
    # for clade, color in clade_color_map.items():
    #     patch = mpatches.Patch(color=color, label=clade)
    #     patch_list.append(patch)
    # plt.gcf().legend(handles=patch_list)
    #
    # plt.savefig(path_vaf_heatmap, dpi=600)


    # plot heatmap
    plt.style.use('fivethirtyeight')
    g = sns.clustermap(vaf_table, cmap="viridis", figsize=(30, 15),
                       xticklabels=True, yticklabels=True)
    ax_heatmap = g.ax_heatmap
    # ax_heatmap.set_xticklabels(ax_heatmap.get_xticklabels())
    # ax_heatmap.set_yticklabels(ax_heatmap.get_yticklabels())

    patch_list = []
    for clade, color in clade_color_map.items():
        patch = mpatches.Patch(color=color, label=clade)
        patch_list.append(patch)
    plt.gcf().legend(handles=patch_list)

    for xtick_label in ax_heatmap.get_yticklabels():
        sample_name = xtick_label.get_text()
        clade = nextclade.loc[sample_name, 'Nextclade_pango']
        label_color = clade_color_map[clade]
        xtick_label.set_color(label_color)

    plt.savefig(path_vaf_heatmap, dpi=600)








