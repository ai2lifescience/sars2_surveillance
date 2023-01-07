import collections
from pathlib import Path
import os
from collections import defaultdict
import shutil
import errno
import gzip
import re



def merge_lane(input_path):

    filenames = os.listdir(input_path)
    print(filenames)
    sample_dic = defaultdict(list)
    for i in filenames:
        sample_name = i.split('.')[0]
        sample_name = sample_name.split('_')
        sample_name = [sample_name[0], sample_name[1], sample_name[3]]
        sample_name = '_'.join(sample_name)
        sample_dic[sample_name].append(i)
    #print(sample_dic)

    output_path = input_path.parent / 'combine'
    os.makedirs(output_path, exist_ok=True)

    for k, v in sample_dic.items():

        scr = [input_path / i for i in v]
        dst = output_path / f'{k}.fastq.gz'

        with open(dst, 'wb') as wfp:
            for fn in scr:
                with open(fn, 'rb') as rfp:
                    shutil.copyfileobj(rfp, wfp)

        print(f'-----{k}-----')
        print('from:\n', scr)
        print('to:\n', dst)


if __name__ == '__main__':

    # create combine folder than merge sequencing file in Raw_data from each lane into it
    path = Path('/media/data/fuhaoyi/sequencing_data/20230103/Raw_Data')
    merge_lane(path)









