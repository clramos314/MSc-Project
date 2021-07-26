import io
import numpy as np
import pandas as pd

sample_id = 'TCGA-05-4244-01A-01D-1105-08'
frames = []


def read_vcf(path):
    d_vcf = {}
    headers = []
    lines = []
    with open(path, 'r') as f:
        for line in f:
            if not line.startswith('##'):
                lines.append(line)
            else:
                headers.append(line.replace('\n', ''))

    d_vcf['headers'] = pd.DataFrame(data=headers)

    d_vcf['body'] = pd.read_csv(io.StringIO(''.join(lines)),
                                dtype={'#CHROM': str, 'POS': int, 'ID': str, 'REF': str, 'ALT': str,
                                       'QUAL': str, 'FILTER': str, 'INFO': str},
                                sep='\t'
                                )
    return d_vcf


if __name__ == '__main__':
    dict_final = {}
    dict_vcf = read_vcf('/home/cleon/extdata/tumor.vep_filtered.vcf')
    header_vcf = dict_vcf['headers']
    df_filtered = dict_vcf['body']

    df_clusters = pd.read_csv(f'/home/cleon/extdata/tumor_{sample_id}_clustered.tsv',
                              sep="\t",
                              low_memory=False)

    arr = np.array(df_clusters['cluster_id'].unique())
    for idx, x in np.ndenumerate(arr):
        df_cluster = pd.DataFrame(df_clusters.loc[df_clusters['cluster_id'] == 0])

        for index, row in df_cluster.iterrows():
            mutation_id = row['mutation_id']
            chro = mutation_id.split(':')[1]
            pos = mutation_id.split(':')[2]

            rslt_df = df_filtered.loc[df_filtered['POS'] == int(pos)]
            if pd.DataFrame(rslt_df).size > 0:
                frames.append(rslt_df)

    header_vcf.to_csv(f'/home/cleon/extdata/input_pandrugs_{sample_id}.vcf',
                      sep='\t',
                      index=False,
                      header=False)

    df_frames = pd.concat(frames)
    df_frames.to_csv(f'/home/cleon/extdata/input_pandrugs_{sample_id}.vcf',
                     sep='\t',
                     index=False,
                     header=True,
                     mode='a')
