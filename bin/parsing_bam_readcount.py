import argparse
import pandas
import pandas as pd


def parse_args():
    Description = "Parsing bam-readcount output."
    Epilog = """Example usage: python3 parsing_bam_readcount.py <file_in> <file_out>"""

    parser = argparse.ArgumentParser(description=Description, epilog=Epilog)
    parser.add_argument("file_in", help="Input bam_readcount file.")
    parser.add_argument("file_out", help="Full path to output csv file.")

    # parser.add_argument(
    #     "-f",
    #     "--fasta",
    #     type=str,
    #     default=None,
    #     help="Fasta file used in mapping and variant calling for vcf header reference genome lenght info.",
    # )

    return parser.parse_args()



if __name__ == '__main__':

    # Process args
    args = parse_args()

    # Per-base/indel data fields
    base_fields = {
        'base': str,
        'count': int,
        'avg_mapping_quality': float,
        'avg_basequality': float,
        'avg_se_mapping_quality': float,
        'num_plus_strand': int,
        'num_minus_strand': int,
        'avg_pos_as_fraction': float,
        'avg_num_mismatches_as_fraction': float,
        'avg_sum_mismatch_qualities': float,
        'num_q2_containing_reads': int,
        'avg_distance_to_q2_start_in_q2_reads': float,
        'avg_clipped_length': float,
        'avg_distance_to_effective_3p_end': float
    }


    # Open the bam-readcount output file and read it line by line
    # Note that the output has no header, so we consume every line
    with open(args.file_in) as in_fh:

        readcount_table = []

        for line in in_fh:
            # Strip newline from end of line
            line = line.strip()
            # Fields are tab-separated, so split into a list on \t
            fields = line.split('\t')
            # The first four fields contain overall information about the position
            chrom = fields[0]  # Chromosome/reference name
            position = int(fields[1])  # Position (1-based)
            reference_base = fields[2]  # Reference base
            depth = int(fields[3])  # Depth of coverage
            # The remaining fields are data for each base or indel
            # Iterate over each base/indel

            readcount_list = []
            vaf_list = []
            for base_data_string in fields[4:]:
                # We will store per-base/indel data in a dict
                base_data = {}
                # Split the base/indel data on ':'
                base_values = base_data_string.split(':')
                # Iterate over each field of the base/indel data
                # Note that this relies on Python 3.5+ to keep the keys in order
                for i, base_field in enumerate(base_fields.keys()):
                    # Store each field of data, converting to the appropriate
                    # data type
                    base_data[base_field] = base_fields[base_field](base_values[i])

                # Do something with base_data, store it, filter it, etc.

                if base_data['base'] != '=':
                    base = base_data['base']
                    count = base_data['count']
                    if depth != 0:
                        vaf = base_data['count'] / depth
                    else:
                        vaf = 0
                    readcount_list.append(f'{base}:{count}')
                    vaf_list.append(f'{base}:{vaf}')

            readcount_str = ' '.join(readcount_list)
            vaf_str = ' '.join(vaf_list)

            readcount_table.append([chrom, position, reference_base,
                                    depth, readcount_str, vaf_str])

    readcount_df = pd.DataFrame(data=readcount_table, columns=['chrom', 'position',
                                                               'reference_base', 'depth',
                                                               'readcount_str', 'vaf'])

    readcount_df.to_csv(args.file_out, index=False)








