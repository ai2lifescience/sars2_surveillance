import collections
from pathlib import Path
import os
from collections import defaultdict
import shutil
import errno
import gzip
import re


def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError as exc:  # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise

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
    mkdir_p(output_path)
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

    path = Path('/media/data/fuhaoyi/sequencing_data/20230103/Raw_Data')
    merge_lane(path)









