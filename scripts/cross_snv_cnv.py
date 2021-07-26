import datetime
import io
import math

import pandas as pd

l_mutation_id = []
l_sample_id = []
l_ref_counts = []
l_alt_counts = []
l_normal_cn = []
l_minor_cn = []
l_major_cn = []
l_tumour_content = []
l_error_rate = []

sample_id = 'TCGA-05-4244-01A-01D-1105-08'


def read_vcf(path):
    with open(path, 'r') as f:
        lines = [l for l in f if not l.startswith('##')]
    return pd.read_csv(
        io.StringIO(''.join(lines)),
        dtype={'#CHROM': str, 'POS': int, 'ID': str, 'REF': str, 'ALT': str,
               'QUAL': str, 'FILTER': str, 'INFO': str},
        sep='\t'
    ).rename(columns={'#CHROM': 'CHROM'})


def cross_rows(snv_row, cnv_row):
    # list of mutation_id, ref_counts, var_counts, normal_cn, minor_cn, major_cn
    aux_mutation_id = sample_id + ":" + str(snv_row["CHROM"]) + ":" + str(snv_row["POS"])
    if not l_mutation_id.__contains__(aux_mutation_id):
        list.append(l_mutation_id,
                    sample_id + ":" + str(snv_row["CHROM"]) + ":" + str(snv_row["POS"]))
        list.append(l_sample_id, sample_id)

        if snv_row["ID"] != 'Nan':
            composed_id = str.split(snv_row["ID"], sep="/")
            if len(composed_id) == 2:
                if composed_id[0] != 'Nan':
                    list.append(l_ref_counts, composed_id[0])
                    if composed_id[1] != 'Nan':
                        l_var_counts_composed = str.split(composed_id[1], sep=" ")
                        list.append(l_alt_counts, l_var_counts_composed[0])
            else:
                list.append(l_ref_counts, '0')
                list.append(l_alt_counts, '0')
            list.append(l_normal_cn, int(cnv_row['copy_number']))
            list.append(l_minor_cn, int(cnv_row['min_copy_number']))
            list.append(l_major_cn, int(cnv_row['max_copy_number']))
            list.append(l_tumour_content, 1.0)
            list.append(l_error_rate, 0.001)
        else:
            list.append(l_normal_cn, '0')
            list.append(l_minor_cn, '0')
            list.append(l_major_cn, '0')
            list.append(l_tumour_content, 1.0)
            list.append(l_error_rate, 0.001)


if __name__ == '__main__':

    snv_input = read_vcf('/home/cleon/extdata/tumor.vep_filtered.vcf')
    cnv_input = pd.read_csv('/home/cleon/extdata/TCGA-38-4629/CNVs/TCGA-38-4629_copy_number.tsv', sep='\t')

    # for each POS in SNV file
    for indexSNV, rowSNV in snv_input.iterrows():
        # find if it's within any CNV
        for indexCNV, rowCNV in cnv_input.iterrows():
            pos = rowSNV["POS"]
            # print(rowCNV.values)
            start = rowCNV['start']
            end = rowCNV['end']
            if rowSNV['CHROM'] != rowCNV['chromosome']:
                break
            if pos >= start:
                if pos <= end:
                    cross_rows(rowSNV, rowCNV)

    # dictionary of lists
    # mutation_id ref_counts var_counts normal_cn minor_cn major_cn
    dict2pyclone = {'mutation_id': l_mutation_id,
                    'sample_id': l_sample_id,
                    'ref_counts': l_ref_counts,
                    'alt_counts': l_alt_counts,
                    'normal_cn': l_normal_cn,
                    'minor_cn': l_minor_cn,
                    'major_cn': l_major_cn,
                    'tumour_content': l_tumour_content,
                    'error_rate': l_error_rate}

    df_output = pd.DataFrame(dict2pyclone)

    # saving the dataframe
    df_output.to_csv('/home/cleon/extdata/tumor_' + sample_id + '_pyclone.tsv',
                     sep='\t',
                     index=False)

