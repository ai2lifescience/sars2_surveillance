import collections
from pathlib import Path
import os
from collections import defaultdict
import shutil
import errno
import gzip
import re


def check_gzip(path):
    filenames = os.listdir(path)
    filenames = sorted(filenames)
    corrupted = {}
    for name in filenames:
        with gzip.open(path/name, 'rb') as f:
            try:
                line = f.readline()
            except:
                corrupted[name] = True
                print('erro', name)
            else:
                corrupted[name] = False

    counter = collections.Counter(corrupted.values())
    # print(counter)


    return corrupted


def get_sample_list(path):

    corrupted = check_gzip(path)

    filenames = os.listdir(path)
    filenames_filter = [i for i in filenames if corrupted[i]==False]

    sample_list = [i.split('.')[0] for i in filenames_filter]
    sample_list = ['_'.join(i.split('_')[0:-1]) for i in sample_list]
    # paired samples
    counter = collections.Counter(sample_list)
    #print(counter)
    sample_list = [k for k, v in counter.items() if v>1]
    print(sample_list)
    print(len(sample_list))

    return sample_list

def renamefiles(path):
    filenames = os.listdir(path)

    for i in filenames:
        oldname = path/i
        newname = i.replace('_', '___')
        newname = path/newname
        print(oldname, newname)
        os.rename(oldname, newname)

    # filenames = os.listdir(path)
    # print(filenames)



if __name__ == '__main__':

    path = Path('/media/data/fuhaoyi/sequencing_data/20221207/combine')
    sample_list = get_sample_list(path)

    path = Path('/media/data/fuhaoyi/sequencing_data/20221211/combine')
    sample_list = get_sample_list(path)

    path = Path('/media/data/fuhaoyi/sequencing_data/20230103/combine')
    sample_list = get_sample_list(path)








