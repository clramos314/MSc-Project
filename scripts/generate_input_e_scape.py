import os
import re
import pdb
import numpy as np
import pandas as pd
import tables

dire = "/home/cleon/workspace_mapscape/mapscape/data/input/"
lists = ["citup_output_20ID00738.h5"]


if __name__ == '__main__':

    print('start')

    for filename in lists:
        print(filename)
        patient = filename.split(".")[0].split("_")[2]
        samples = pd.read_csv(dire+"samples_"+patient+".tsv", sep="\t", header=0, low_memory=False)

        inputf = open(dire + filename,'r')

        f = tables.open_file(dire + filename)
        table1 = f.root.results.optimal
        #seleccionar el mejor arbol
        best = list(np.array(table1.index))[list(np.array(table1.values)).index(True)]
        print(best)
        array = np.array(f.root.trees._f_get_child(str(best)).adjacency_list.block0_values)
        np.savetxt(dire + filename.replace(".h5", "_tree.tsv"), array, delimiter='\t', newline='\n', header='source\ttarget', fmt='%d', comments='')

        #para mapscape
        clones = np.array(f.root.trees._f_get_child(str(best)).clone_freq.block0_items)
        freq = np.array(f.root.trees._f_get_child(str(best)).clone_freq.block0_values)
        #sample_list = [samples[samples["index"] == x]["values"].tolist()[0] for x in np.array(f.root.trees._f_get_child(str(best)).clone_freq.axis1)]

        sample_list = []
        for x in np.array(freq):
            sample_list.append(x)

        clone_str = pd.DataFrame(
            {'sample_id': np.repeat(np.array(f.root.trees._f_get_child(str(best)).clone_freq.axis1), len(clones)),
             'clone_id': np.tile(clones, len(sample_list)),
             'clonal_prev': np.concatenate(freq, axis=0)}, columns=['sample_id', 'clone_id', 'clonal_prev'])

        ##clone_str = pd.DataFrame({'timepoint': np.repeat(sample_list, len(clones)), 'clone_id': np.tile(clones, len(samples)), 'clonal_prev': np.concatenate(freq, axis=0)}, columns=['timepoint', 'clone_id','clonal_prev'])

        clone_str.to_csv(dire + filename.replace(".h5", "_clonal_prev.tsv"), index=False, header=True, sep="\t")
        #para obtener las mutaciones tanto para timescape como mapscape
        variant_assignment = pd.DataFrame(
            {'index': np.array(f.root.trees._f_get_child(str(best)).variant_assignment.index), 
            'values': np.array(f.root.trees._f_get_child(str(best)).variant_assignment.values)}, columns=['index', 'values'])
            
        variant_assignment.to_csv(dire + filename.replace(".h5", "_variant_assignment.tsv"), index=False, header=True, sep="\t")

        va_values = np.array(f.root.trees._f_get_child(str(best)).variant_assignment.values)

        #mutations = pd.read_csv("/home/epineiro/Analysis/Paradiff/clonality/leia/pyclone_analysis/"+patient+"/tables/loci.tsv", sep ="\t", header=0, low_memory=False)
        mutations = pd.read_csv("/home/cleon/workspace_mapscape/mapscape/data/input/loci.tsv", sep="\t", header=0, low_memory=False)
        for i in variant_assignment["index"].tolist():
            mutations.loc[mutations['cluster_id'] == i, 'cluster_citup'] = str(variant_assignment[variant_assignment["index"] == i]["values"].tolist()[0])
        mutations = pd.DataFrame(
            {'chrom': mutations["mutation_id"].str.split(':', expand=True)[1],
            'coord': mutations["mutation_id"].str.split(':', expand=True)[2],
            'clone_id': mutations["cluster_citup"],
            'sample_id': mutations["sample_id"],
            'VAF': mutations["variant_allele_frequency"]}, columns=['chrom', 'coord', 'clone_id', 'sample_id', 'VAF'])
    
        mutations.to_csv(dire + filename.replace(".h5", "_mutations.tsv"), index=False, header=True, sep="\t")

    f.close()
    inputf.close()
    print('end')
