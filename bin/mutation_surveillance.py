import pandas as pd

from utils import *

if __name__ == '__main__':

    data = pd.read_csv('/media/data/fuhaoyi/sequencing_data/20230108/output/6tree/20230108_nextclade.csv',
                       sep=';')
    data = data.fillna(value='')
    #print(data)


    def different_mutaion(data, ref):

        ref_lineage = ref[ref['proportion'] > 0.5]
        ref_lineage = ref_lineage.reset_index(drop=True)

        diff_mut = []
        for i, row in data.iterrows():
            aa_sub = row['aaSubstitutions'].split(',')
            aa_del = row['aaDeletions'].split(',')
            aa_mut = aa_sub+aa_del
            ref_mut = ref_lineage['mutation'].to_list()
            diff = set(aa_mut) ^ set(ref_mut)
            # print(diff)
            diff_mut.append(','.join(diff))

        return diff_mut

    path = Path('../src/ref_lineage_mutation')
    ref_lineage = pd.read_csv(path / 'BA.5.2.csv')
    ba_5_2 = different_mutaion(data, ref_lineage)

    path = Path('../src/ref_lineage_mutation')
    ref_lineage = pd.read_csv(path / 'BF.7.csv')
    bf7 = different_mutaion(data, ref_lineage)

    output_data = pd.DataFrame({
        'seqName': data['seqName'],
        'clade': data['clade'],
        'pango': data['Nextclade_pango'],
        'diff_BA.5.2': ba_5_2,
        'diff_BF.7': bf7,
    })
    print(output_data)
    output_data.to_csv('../results/diff_mutations.csv', index=False)








