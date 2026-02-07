import anndata
import random
import numpy as np
import pandas as pd
from tqdm import tqdm
from numpy.random import choice
from scipy import sparse as sp

def generate_simulated_data_scaden(sc_data, outname=None,
                                   d_prior=None,
                                   n=500, samplenum=5000,
                                   random_state=None, sparse=True, sparse_prob=1,
                                   rare=False, rare_percentage=0.4, id_col = "patient_id",
                                   cell_type = None, add_realdat = False, real_dat = None):
    # sc_data should be a cell*gene matrix, no null value, txt file, sep='\t'
    # index should be cell names
    # columns should be gene labels
    print('Reading single-cell dataset, this may take 1 min')
    if isinstance(sc_data, str): 
        if '.h5ad' in sc_data:
            print('You are using H5AD format data, please make sure "CellType" occurs in the adata.obs')
            sc_data = anndata.read_h5ad(sc_data)

            if cell_type not in sc_data.obs.columns.tolist():
                raise Exception('Please make sure ' + cell_type + ' column present in the adata.obs.')
            
            if id_col not in sc_data.obs.columns.tolist():
                raise Exception('Please make sure ' + id_col + ' column present in the adata.obs.')

            patientid = sc_data.obs[id_col].unique()
            if isinstance(sc_data.X, np.ndarray):
                pass
            else:
                sc_data.X = sc_data.X.todense()
    elif isinstance(sc_data, anndata.AnnData):
        print('You are using H5AD format data, please make sure "CellType" occurs in the adata.obs')
        
        if cell_type not in sc_data.obs.columns.tolist():
            raise Exception('Please make sure ' + cell_type + ' column present in the adata.obs.')
        
        if id_col not in sc_data.obs.columns.tolist():
                raise Exception('Please make sure ' + id_col + ' column present in the adata.obs.')
        
        patientid = sc_data.obs[id_col].unique()
        if isinstance(sc_data.X, np.ndarray):
            pass
        else:
            sc_data.X = sc_data.X.todense()
    else:
        raise Exception("Please check the format of single-cell data!")
    print('Reading dataset is done')

    # set min included cell type num based on data.
    min_count = []
    for i in patientid:
        sc_data_tmp = sc_data[sc_data.obs[id_col] == i]
        min_count.append(len(sc_data_tmp.obs[cell_type].unique()))
    min_ct_num = min(min_count)

    # for pid in patientid:
    first_p = True

    for pid in patientid:

        sc_data_patient = sc_data[sc_data.obs[id_col] == pid]

        # print(sc_data_patient.X)
        # print(isinstance(sc_data_patient.X, sp.csc_matrix))
        if isinstance(sc_data_patient.X, sp.csc_matrix):
            sc_data_patient = pd.DataFrame(sc_data_patient.X.todense(), index=sc_data_patient.obs[cell_type], columns=sc_data_patient.var.index) # sc_data_patient.X.todense()
        else:
            sc_data_patient = pd.DataFrame(sc_data_patient.X, index=sc_data_patient.obs[cell_type], columns=sc_data_patient.var.index) # sc_data_patient.X.todense()
        # sc_data_patient = pd.DataFrame(sc_data_patient.X.todense(), index=sc_data_patient.obs[cell_type], columns=sc_data_patient.var.index)
        
        sc_data_patient.dropna(inplace=True)
        sc_data_patient[cell_type] = sc_data_patient.index
        sc_data_patient.index = range(len(sc_data_patient))

        num_celltype = len(sc_data_patient[cell_type].value_counts())
        genename = sc_data_patient.columns[:-1]

        celltype_groups = sc_data_patient.groupby(cell_type).groups
        sc_data_patient.drop(columns=cell_type, inplace=True)

        sc_data_patient = anndata.AnnData(sc_data_patient)
        # sc.pp.normalize_total(sc_data_patient, target_sum=1e4)

        # use ndarray to accelerate
        # change to C_CONTIGUOUS, 10x faster
        sc_data_patient = sc_data_patient.X
        sc_data_patient = np.ascontiguousarray(sc_data_patient, dtype=np.float32)

        # make random cell proportions
        if random_state is not None and isinstance(random_state, int):
            print('You specified a random state, which will improve the reproducibility.')

        # make the dictionary
        for key, value in celltype_groups.items():
            celltype_groups[key] = np.array(value)

        prop = np.zeros((samplenum, num_celltype))

        for i in range(int(prop.shape[0] * sparse_prob)):
    
            non_zero_cell = random.randint(min_ct_num, num_celltype) # we choose a random number of cell types between 5 and maximum available cell types for a patient was chosen to included
            indices = np.random.choice(np.arange(prop.shape[1]), replace=False, size=non_zero_cell) # a random fraction will be assigned to a random cell type.
            selected_cell_type_C = np.random.choice(indices)
            selected_cell_type_C_prop = random.uniform(0,1)

            prop[i, selected_cell_type_C] = selected_cell_type_C_prop
            indices_withoutC = np.setdiff1d(indices, selected_cell_type_C)
            prop[i, indices_withoutC] = np.random.uniform(0,1, len(indices_withoutC)) # randomly generated fractions for other cell types, normalised them between [0,1]
            prop[i, indices_withoutC] = normalize(prop[i, indices_withoutC], 1-prop[i, selected_cell_type_C]) # normalised them againbetween [0, 1- fc] to ensure cell-type proportions sum to 1

        # precise number for each celltype
        cell_num = np.floor(n * prop)

        # precise proportion based on cell_num
        prop = cell_num / np.sum(cell_num, axis=1).reshape(-1, 1)

        # start sampling
        sample = np.zeros((prop.shape[0], sc_data_patient.shape[1]))
        allcellname = celltype_groups.keys()
        print('Sampling cells to compose pseudo-bulk data')
        for i, sample_prop in tqdm(enumerate(cell_num)):
            for j, cellname in enumerate(allcellname):
                select_index = choice(celltype_groups[cellname], size=int(sample_prop[j]), replace=True)
                sample[i] += sc_data_patient[select_index].sum(axis=0)
                                            
        prop = pd.DataFrame(prop, columns=celltype_groups.keys())

        if first_p:
            simudata = anndata.AnnData(X=sample,
                                       obs=prop,
                                       var=pd.DataFrame(index=genename))
            first_p = False
        else:
            simudata_tmp = anndata.AnnData(X=sample,
                                           obs=prop,
                                           var=pd.DataFrame(index=genename))
            simudata = anndata.concat([simudata, simudata_tmp], join='outer')

    if add_realdat == True:
        real_dat_h5ad = anndata.read_h5ad(real_dat)
        simudata = anndata.concat([simudata, real_dat_h5ad])

    print('Sampling is done')
    if outname is not None:
        simudata.write_h5ad(outname + '.h5ad')
    simudata.obs.fillna(0, inplace=True)

    return simudata

def normalize(arr, target_sum):
    return arr / np.sum(arr) * target_sum