import anndata
import numpy as np
import pandas as pd
from tqdm import tqdm
from numpy.random import choice

def generate_simulated_data_tape(sc_data, outname=None,
                                 d_prior=None,
                                 n=500, samplenum=5000,
                                 random_state=None, sparse=True, sparse_prob=0.5,
                                 rare=False, rare_percentage=0.4,
                                 cell_type = None, add_noise = False, noise_level=0.1,
                                 return_noise_free = False, add_realdat = False, real_dat = None):
    # sc_data should be a cell*gene matrix, no null value, txt file, sep='\t'
    # index should be cell names
    # columns should be gene labels
    print('Reading single-cell dataset, this may take 1 min')
    if isinstance(sc_data, str):
        if '.txt' in sc_data:
            sc_data = pd.read_csv(sc_data, index_col=0, sep='\t')
            sc_data.dropna(inplace=True)
            sc_data['celltype'] = sc_data.index
            sc_data.index = range(len(sc_data))
        elif '.h5ad' in sc_data:
            print('You are using H5AD format data, please make sure' + cell_type + 'occurs in the adata.obs')
            sc_data = anndata.read_h5ad(sc_data)
            if isinstance(sc_data.X, np.ndarray):
                pass
            else:
                sc_data.X = sc_data.X.toarray()

            sc_data = pd.DataFrame(sc_data.X, index=sc_data.obs[cell_type], columns=sc_data.var.index)
            sc_data.dropna(inplace=True)
            sc_data['celltype'] = sc_data.index
            sc_data.index = range(len(sc_data))
            print("Finished Reading.")
    elif type(sc_data) is pd.DataFrame:
        sc_data.dropna(inplace=True)
        sc_data['celltype'] = sc_data.index
        sc_data.index = range(len(sc_data))
    elif isinstance(sc_data, anndata.AnnData):
        print('You are using H5AD format data, please make sure' + cell_type + 'occurs in the adata.obs')
        if isinstance(sc_data.X, np.ndarray):
            pass
        else:
            sc_data.X = sc_data.X.toarray()

        if isinstance(sc_data.X, np.ndarray):
            sc_data = pd.DataFrame(sc_data.X, index=sc_data.obs[cell_type], columns=sc_data.var.index)
        else:
            sc_data = pd.DataFrame(sc_data.X.todense(), index=sc_data.obs[cell_type], columns=sc_data.var.index)
        sc_data.dropna(inplace=True)
        sc_data['celltype'] = sc_data.index
        sc_data.index = range(len(sc_data))
    else:
        raise Exception("Please check the format of single-cell data!")
    print('Reading dataset is done')

    num_celltype = len(sc_data['celltype'].value_counts())
    genename = sc_data.columns[:-1]

    celltype_groups = sc_data.groupby('celltype').groups
    sc_data.drop(columns='celltype', inplace=True)

    ### normalize with scanpy
    print('Normalizing raw single cell data with scanpy.pp.normalize_total')
    sc_data = anndata.AnnData(sc_data)
    # sc.pp.normalize_total(sc_data, target_sum=1e4)

    # use ndarray to accelerate
    # change to C_CONTIGUOUS, 10x faster
    sc_data = sc_data.X
    sc_data = np.ascontiguousarray(sc_data, dtype=np.float32)
    # make random cell proportions

    if random_state is not None and isinstance(random_state, int):
        print('You specified a random state, which will improve the reproducibility.')

    if d_prior is None:
        print('Generating cell fractions using Dirichlet distribution without prior info (actually random)')
        if isinstance(random_state, int):
            np.random.seed(random_state)
        prop = np.random.dirichlet(np.ones(num_celltype), samplenum)
        print('RANDOM cell fractions is generated')
    elif d_prior is not None:
        print('Using prior info to generate cell fractions in Dirichlet distribution')
        assert len(d_prior) == num_celltype, 'dirichlet prior is a vector, its length should equals ' \
                                             'to the number of cell types'
        if isinstance(random_state, int):
            np.random.seed(random_state)
        prop = np.random.dirichlet(d_prior, samplenum)
        print('Dirichlet cell fractions is generated')

    # make the dictionary
    for key, value in celltype_groups.items():
        celltype_groups[key] = np.array(value)

    prop = prop / np.sum(prop, axis=1).reshape(-1, 1)
    # sparse cell fractions
    if sparse:
        print("You set sparse as True, some cell's fraction will be zero, the probability is", sparse_prob)
        ## Only partial simulated data is composed of sparse celltype distribution
        for i in range(int(prop.shape[0] * sparse_prob)):
            indices = np.random.choice(np.arange(prop.shape[1]), replace=False, size=int(prop.shape[1] * sparse_prob))
            prop[i, indices] = 0

        prop = prop / np.sum(prop, axis=1).reshape(-1, 1)

    if rare:
        print(
            'You will set some cell type fractions are very small (<3%), '
            'these celltype is randomly chosen by percentage you set before.')
        ## choose celltype
        np.random.seed(0)
        indices = np.random.choice(np.arange(prop.shape[1]), replace=False, size=int(prop.shape[1] * rare_percentage))
        prop = prop / np.sum(prop, axis=1).reshape(-1, 1)

        for i in range(int(0.5 * prop.shape[0]) + int(int(rare_percentage * 0.5 * prop.shape[0]))):
            prop[i, indices] = np.random.uniform(0, 0.03, len(indices))
            buf = prop[i, indices].copy()
            prop[i, indices] = 0
            prop[i] = (1 - np.sum(buf)) * prop[i] / np.sum(prop[i])
            prop[i, indices] = buf

    # precise number for each celltype
    cell_num = np.floor(n * prop)

    # precise proportion based on cell_num
    prop = cell_num / np.sum(cell_num, axis=1).reshape(-1, 1)

    # start sampling
    sample = np.zeros((prop.shape[0], sc_data.shape[1]))
    sample_pure = np.zeros((prop.shape[0], sc_data.shape[1]))
    allcellname = celltype_groups.keys()
    print('Sampling cells to compose pseudo-bulk data')
    if add_noise:
        for i, sample_prop in tqdm(enumerate(cell_num)):
            tmp_index = []
            for j, cellname in enumerate(allcellname):
                select_index = choice(celltype_groups[cellname], size=int(sample_prop[j]), replace=True)
                sample[i] += sc_data[select_index].sum(axis=0)
                tmp_index = tmp_index + list(select_index) # get the full tmp index

                sample_pure[i] += sc_data[select_index].sum(axis=0)

            if add_noise:
                means = sc_data[tmp_index].mean(axis=0)
                stds = sc_data[tmp_index].std(axis=0)
                noises = np.random.normal(loc=means, scale=stds)
                noises_added = np.abs(noise_level*noises)
                sample[i] += noises_added
    else:
        for i, sample_prop in tqdm(enumerate(cell_num)):
            for j, cellname in enumerate(allcellname):
                select_index = choice(celltype_groups[cellname], size=int(sample_prop[j]), replace=True)
                sample[i] += sc_data[select_index].sum(axis=0)

    prop = pd.DataFrame(prop, columns=celltype_groups.keys())
    simudata = anndata.AnnData(X=sample,
                               obs=prop,
                               var=pd.DataFrame(index=genename))
    
    if return_noise_free:
        simudata_nonoise = anndata.AnnData(X=sample_pure,
                                            obs=prop,
                                            var=pd.DataFrame(index=genename))
        
    if add_realdat == True:
        real_dat_h5ad = anndata.read_h5ad(real_dat)
        simudata = anndata.concat([simudata, real_dat_h5ad])
    
    print('Sampling is done')
    if outname is not None:
        simudata.write_h5ad(outname + '.h5ad')
        
    if return_noise_free:
        return simudata, simudata_nonoise
    else:
        return simudata