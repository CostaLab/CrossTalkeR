## CellPhoneDB CookBook

*by James Nagai*

### 0- Installation



```pip3 install cellphonedb```



### 1 - Step Extract data from Seurat Object

```R
require(EWCE)
require(tibble)
require(biomaRt)
require(tidyr)
require(dplyr)
# Basic function to convert mouse to human gene names
data_ <- readRDS('path_to_data')
alldata <- NormalizeData(alldata, normalization.method = "LogNormalize", scale.factor = 10000)
allgenes <- rownames(alldata)
matrix1 <- as.data.frame(alldata@assays$RNA@data)
matrix1 <- matrix1[rowSums(matrix1[,2:dim(matrix1)[2]])!=0,]

### If you are using a mouse data, then its needed to convert the gene names to human orthologs
human = useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl")
mouse = useEnsembl("ensembl", dataset = "mmusculus_gene_ensembl")
genesV2 = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = rownames(alldata@assays$RNA@data) , mart = mouse, attributesL = c("hgnc_symbol","hgnc_id",'ensembl_gene_id'), martL = human, uniqueRows=T)
print(head(genesV2))
matrix1 <- matrix1[match(genesV2$MGI.symbol,rownames(alldata),nomatch=F),]
matrix1$gene <- genesV2$Gene.stable.ID
#rownames(matrix1) <-  genesV2$Gene.stable.ID

### Subseting the matrix
s1 <- grepl('state1',alldata@meta.data$cond)
s2 <- grepl('state2',alldata@meta.data$cond)
s1[match('gene',colnames(matrix1))] <- TRUE
s2[match('gene',colnames(matrix1))] <- TRUE

## Checking the dimensions
print(dim(matrix1[,s1]))
print(dim(matrix1[,s2]))


## If the cluster names are categorical, you will need to convert it to numerical
alldata@meta.data$sub_cell_type <- as.factor(alldata@meta.data$sub_cell_type)
print(levels(alldata@meta.data$sub_cell_type))
levels(alldata@meta.data$sub_cell_type) <- 1:length(levels(alldata@meta.data$sub_cell_type))
print(1:length(levels(alldata@meta.data$sub_cell_type)))
alldata@meta.data$sub_cell_type <- as.numeric(alldata@meta.data$sub_cell_type)


write.table(matrix1[,s1], 's1_filtered_hcount.csv',row.names=T,sep=',')
write.table(matrix1[,s2], 's2_filtered_hcount.csv',row.names=T,sep=',')
metadata <- data.frame(cells=rownames(alldata@meta.data[grepl('state1',alldata@meta.data$stim),]),cluster=alldata@meta.data$sub_cell_type[grepl('state1',alldata@meta.data$stim)])
metadata_s2 <- data.frame(cells=rownames(alldata@meta.data[!grepl('state1',alldata@meta.data$stim),]),cluster=alldata@meta.data$sub_cell_type[!grepl('state1',alldata@meta.data$stim)]) ## Just negate grepl('state1',alldata@meta.data$stim),]
print('Writing Metadata')
write.csv(metadata, 's1_filtered_meta.csv', row.names=FALSE)
write.csv(metadata_tac, 's2_filtered_meta.csv', row.names=FALSE)```
```

**Note that the gene ID needs to be unique, you will need to use an approach to combine multiple mapped orthologs (sum, max or bitwiseor )**

### Run CellPhoneDB

** The parameters list is available at https://github.com/Teichlab/cellphonedb

#### Ensembl based ID

```
#! /bin/bash
mkdir s1 s2  # creating the output folders

cellphonedb method statistical_analysis s1_filtered_meta.csv  s1_filtered_hcount.csv  --threads 30 --output-path s1/

cellphonedb method statistical_analysis s2_filtered_meta.csv  s2_filtered_hcount.csv  --threads 30 --output-path s2/

```

#### HUGO based ID

```
#! /bin/bash
mkdir s1 s2  # creating the output folders

cellphonedb method statistical_analysis s1_filtered_meta.csv --counts-data hgnc_symbol s1_filtered_hcount.csv  --threads 30 --output-path s1/ 

cellphonedb method statistical_analysis s2_filtered_meta.csv --counts-data hgnc_symbol s2_filtered_hcount.csv  --threads 30 --output-path s2/ 

```

### Extracting LR

```python
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
def cpdb2df(data,clsmapping):
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
            if float(data.iloc[i,j+12]) !=0.0:
                df_data['Ligand'].append(pair[0])
                df_data['Receptor'].append(pair[1])
                df_data['Ligand.Cluster'].append(clsmapping[c_pair[0]])
                df_data['Receptor.Cluster'].append(clsmapping[c_pair[1]])
                df_data['isReceptor_fst'].append(list(data['receptor_a'])[i])
                df_data['isReceptor_scn'].append(list(data['receptor_b'])[i])
                df_data['MeanLR'].append(data.iloc[i,j+12])
    data_final = pd.DataFrame.from_dict(df_data)s1_filtered_corrected.csv
    return(data_final)
            
def cpdb2df_nocls(data):
    '''
   		When the cluster name is used on CPDB
    '''
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
            if float(data.iloc[i,j+12]) !=0.0:
                df_data['Ligand'].append(pair[0])
                df_data['Receptor'].append(pair[1])
                df_data['Ligand.Cluster'].append(c_pair[0])
                df_data['Receptor.Cluster'].append(c_pair[1])
                df_data['isReceptor_fst'].append(list(data['receptor_a'])[i])
                df_data['isReceptor_scn'].append(list(data['receptor_b'])[i])
                df_data['MeanLR'].append(data.iloc[i,j+12])
    data_final = pd.DataFrame.from_dict(df_data)s1_filtered_corrected.csv
    return(data_final)
            
s1 = pd.read_csv('./s1/significant_means.txt',sep='\t')
s2 = pd.read_csv('./s2/ significant_means.txt',sep='\t')
#dict with the mapping
num_to_clust = {'1':'Cluters 1',
                '2':'Cluters 2',
                '3':'Cluters 3',
                '4':'Cluters 4',
                '5':'Cluters 5'}

s1_filtered = cpdb2df(s1,num_to_clust)
s2_filtered = cpdb2df(s2,num_to_clust)
s1_filtered = correct_lr(s1_filtered)
s2_filtered = correct_lr(s2_filtered)

s1_filtered.to_csv('s1_filtered_corrected.csv')
s2_filtered.to_csv('s2_filtered_corrected.csv')

```

### CrossTalkeR

````R
library('CrossTalkeR')

# the method always consider the first path as control: the multiple control case will be handle soon
paths <- c('CTR' = 's1_filtered_corrected.csv', 
           'EXP' = 's1_filtered_corrected.csv',
           'EXP1' = 's1_filtered_corrected.csv',)
# Selected gene list     
genes <- c('TGFB1', 'PF4')

# Generating the report and the object

data <- generate_report(paths=paths, # paths list
						selected_genes=genes, # Selected list
						output='/home/nagai/Documents/', # output path
						threshold = 0, # threshold of prune edges 0=keep all
						out_file='All_Nils.html' report name
						)

````



