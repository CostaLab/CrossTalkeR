import pandas as pd
def correct_lr(data):
    '''
    Invert the RL to LR and R1R2 to r2>r1
    '''
    import pandas as pd
    def swap(a,b): return b,a
    data = data.to_dict('index')
    for k,v in data.items():
        if v['isReceptor_fst'] and v['isReceptor_scn']:
            v['isReceptor_fst'],v['isReceptor_scn'] = swap(v['isReceptor_fst'],v['isReceptor_scn'])
            v['Ligand'],v['Receptor'] = swap(v['Ligand'],v['Receptor'])
            v['Ligand.Cluster'],v['Receptor.Cluster'] = swap(v['Ligand.Cluster'],v['Receptor.Cluster'])
        elif v['isReceptor_fst'] and not v['isReceptor_scn']:
            v['isReceptor_fst'],v['isReceptor_scn'] = swap(v['isReceptor_fst'],v['isReceptor_scn'])
            v['Ligand'],v['Receptor'] = swap(v['Ligand'],v['Receptor'])
            v['Ligand.Cluster'],v['Receptor.Cluster'] = swap(v['Ligand.Cluster'],v['Receptor.Cluster'])
    res_df = pd.DataFrame.from_dict(data,orient='index')
    return (res_df)
def cpdb2df(data):
    data = data.fillna(0)
    df_data = {}
    df_data['Ligand'] = []
    df_data['Receptor'] = []
    df_data['Ligand.Cluster'] = []
    df_data['Receptor.Cluster'] = []
    df_data['isReceptor_fst'] = []
    df_data['isReceptor_scn'] = []
    df_data['MeanLR'] = []
    for i in range(data.shape[0]):
        pair = list(data['interacting_pair'])[i].split('_')
        for j in range(data.iloc[:,12:].shape[1]):
            c_pair = list(data.columns)[j+12].split('|')
            if int(data.iloc[i,j+12]) !=0.0:
                df_data['Ligand'].append(pair[0])
                df_data['Receptor'].append(pair[1])
                df_data['Ligand.Cluster'].append(c_pair[0])
                df_data['Receptor.Cluster'].append(c_pair[1])
                df_data['isReceptor_fst'].append(list(data['receptor_a'])[i])
                df_data['isReceptor_scn'].append(list(data['receptor_b'])[i])
                df_data['MeanLR'].append(data.iloc[i,j+12])
    data_final = pd.DataFrame.from_dict(df_data)
    return(data_final)

import os
os.chdir('/home/nagai/Documents/sarscov/LR')
s1 = pd.read_csv('./CTR_filtered/significant_means.txt',sep='\t')
s2 = pd.read_csv('./COVID_filtered/significant_means.txt',sep='\t')
#dict with the mapping

s1_filtered = cpdb2df(s1)
s2_filtered = cpdb2df(s2)
s1_filtered = correct_lr(s1_filtered)
s2_filtered = correct_lr(s2_filtered)

s1_filtered.to_csv('s1_filtered_corrected.csv')
s2_filtered.to_csv('s2_filtered_corrected.csv')
