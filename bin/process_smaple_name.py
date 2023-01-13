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
        newname = i.replace('fq', 'fastq')
        newname = newname.replace('_1', '_R1')
        newname = newname.replace('_2', '_R2')
        newname = newname.replace('_L5', '')



        newname = path/newname
        print(oldname, newname)
        os.rename(oldname, newname)

    # filenames = os.listdir(path)
    # print(filenames)

def move_files(input_path, output_path):
    filenames, all_files, all_dirs = [], [], []
    for root, dirs, files in os.walk(input_path):
        for file in files:
            filenames.append(file)
            all_files.append(os.path.join(root, file))

        for dir in dirs:
            all_dirs.append(os.path.join(root, dir))
    # print(filenames)
    # print(all_files)

    os.makedirs(output_path, exist_ok=True)
    scr = all_files
    it_scr = iter(scr)
    dst = [output_path / i for i in filenames]
    it_dst = iter(dst)

    for i in range(len(scr)):
        s = next(it_scr)
        d = next(it_dst)

        with open(s, 'rb') as rfp:
            with open(d, 'wb') as wfp:
                shutil.copyfileobj(rfp, wfp)

        print(f'-----{i/len(scr):.2f}-----')
        print('from:\n', s)
        print('to:\n', d)


if __name__ == '__main__':

    # # print sample name list
    # path = Path('/media/data/fuhaoyi/sequencing_data/20221207/combine')
    # sample_list = get_sample_list(path)

    # path = Path('/media/data/fuhaoyi/sequencing_data/20221211/combine')
    # sample_list = get_sample_list(path)

    # path = Path('/media/data/fuhaoyi/sequencing_data/20230103/combine')
    # sample_list = get_sample_list(path)


    # path = Path('/media/data/fuhaoyi/sequencing_data/20230108/'
    #             'PM-XS01KF2022030268-173KA-北京昌平实验室94个外来混合pooling97个子文库包1条lane测序不过滤任务单/'
    #             'ANNO_XS01KF2022030268_PM-XS01KF2022030268-173_2023-01-08_12-12-04_H3GGCDSX5/'
    #             'Rawdata')
    # outpath = Path('/media/data/fuhaoyi/sequencing_data/20230108/combine')
    # # if need
    # # move_files(path, outpath)
    # # renamefiles(outpath)
    # sample_list = get_sample_list(outpath)


    # path = Path('/media/data/fuhaoyi/sequencing_data/20230109/'
    #             'ANNO_XS01KF2022030268_PM-XS01KF2022030268-174_2023-01-09_17-05-28_H3GGCDSX5/'
    #             'Rawdata')
    # outpath = Path('/media/data/fuhaoyi/sequencing_data/20230109/combine')
    # # if need
    # # move_files(path, outpath)
    # # renamefiles(outpath)
    # sample_list = get_sample_list(outpath)


    path = Path('/media/data/fuhaoyi/sequencing_data/20230109_2/combine')
    #renamefiles(path)
    sample_list = get_sample_list(path)






