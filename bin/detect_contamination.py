from utils import *

if __name__ == '__main__':
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

