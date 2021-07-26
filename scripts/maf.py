# Copyright (c) 2016-2017. Mount Sinai School of Medicine
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

from __future__ import print_function, division, absolute_import
import logging

import pandas as pd
from typechecks import require_string

TCGA_PATIENT_ID_LENGTH = 12

MAF_COLUMN_NAMES = [
    'Hugo_Symbol',
    'Entrez_Gene_Id',
    'Center',
    'NCBI_Build',
    'Chromosome',
    'Start_Position',
    'End_Position',
    'Strand',
    'Variant_Classification',
    'Variant_Type',
    'Reference_Allele',
    'Tumor_Seq_Allele1',
    'Tumor_Seq_Allele2',
    'dbSNP_RS',
    'dbSNP_Val_Status',
    'Tumor_Sample_Barcode',
    'Matched_Norm_Sample_Barcode',
    'Match_Norm_Seq_Allele1',
    'Match_Norm_Seq_Allele2',
]


def load_maf_dataframe(path, nrows=None, raise_on_error=True, encoding=None):
    """
    Load the guaranteed columns of a TCGA MAF file into a DataFrame

    Parameters
    ----------
    path : str
        Path to MAF file

    nrows : int
        Optional limit to number of rows loaded

    raise_on_error : bool
        Raise an exception upon encountering an error or log an error

    encoding : str, optional
        Encoding to use for UTF when reading MAF file.
    """
    require_string(path, "Path to MAF")

    n_basic_columns = len(MAF_COLUMN_NAMES)

    # pylint: disable=no-member
    # pylint gets confused by read_csv
    df = pd.read_csv(
        path,
        comment="#",
        sep="\t",
        low_memory=False,
        skip_blank_lines=True,
        header=0,
        nrows=nrows,
        encoding=encoding)

    if len(df.columns) < n_basic_columns:
        error_message = (
            "Too few columns in MAF file %s, expected %d but got  %d : %s" % (
                path, n_basic_columns, len(df.columns), df.columns))
        if raise_on_error:
            raise ValueError(error_message)
        else:
            logging.warn(error_message)

    # check each pair of expected/actual column names to make sure they match
    for expected, actual in zip(MAF_COLUMN_NAMES, df.columns):
        if expected != actual:
            # MAFs in the wild have capitalization differences in their
            # column names, normalize them to always use the names above
            if expected.lower() == actual.lower():
                # using DataFrame.rename in Python 2.7.x doesn't seem to
                # work for some files, possibly because Pandas treats
                # unicode vs. str columns as different?
                df[expected] = df[actual]
                del df[actual]
            else:
                error_message = (
                    "Expected column %s but got %s" % (expected, actual))
                if raise_on_error:
                    raise ValueError(error_message)
                else:
                    logging.warn(error_message)

    return df


if __name__ == '__main__':
    chromosome = []
    start = []
    end = []
    allele = []
    strand = []
    identifier = []

    maf_df = load_maf_dataframe("/home/cleon/extdata/tumor.tsv")
    # iterate through each row and select
    for index, row in maf_df.iterrows():
        chromosome.append(str(row["Chromosome"]))
        start.append(row["Start_Position"])
        end.append(row["End_Position"])
        allele_tmp = row["Reference_Allele"]
        if row["Reference_Allele"] != row["Tumor_Seq_Allele1"]:
            allele_tmp = allele_tmp + "/" + str(row["Tumor_Seq_Allele1"])
        if row["Reference_Allele"] != row["Tumor_Seq_Allele2"]:
            allele_tmp = allele_tmp + "/" + str(row["Tumor_Seq_Allele2"])
        allele.append(allele_tmp)
        strand.append(row["Strand"])
        # identifier puedes añadir la info de profundidad, número de lecturas de referencia y variante, ...
        identifier.append(str(row["t_ref_count"]) + " " + str(row["t_alt_count"]))

    df_output = pd.DataFrame({'chromosome': chromosome,
                              'start': start,
                              'end': end,
                              'allele': allele,
                              'strand': strand,
                              'identifier': identifier})

    # saving the dataframe
    df_output.to_csv('/home/cleon/extdata/subset_pyclone.tsv',
                     header=False,
                     sep='\t',
                     index=False)
    #print(maf_df)
