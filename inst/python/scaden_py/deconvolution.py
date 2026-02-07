import anndata
import pandas as pd
# from scaden_py import generate_simulated_data, ProcessInputData, scaden
from scaden_py.simulation_tape import generate_simulated_data_tape
from scaden_py.simulation_scaden import generate_simulated_data_scaden
from scaden_py.utils import ProcessInputData
from scaden_py.model import scaden

def ScadenDeconvolution(necessary_data, real_bulk, sep='\t', sparse=True,
                        batch_size=128, epochs=128, samplenum = 5000, n = 500, d_prior = None,
                        cell_type="None", generate_sim_method = "scaden", cut_variance = True, id_col = "patient_id",
                        add_noise = False, noise_level = 0.01, add_realdat = False, real_dat = None, markergene = None, quan_norm = False):
    if type(necessary_data) is str:
        postfix = necessary_data.split('.')[-1]
        if postfix == 'txt' or postfix == 'gz':
            if generate_sim_method == "tape":
                simudata = generate_simulated_data_tape(sc_data=necessary_data, samplenum = samplenum, sparse=sparse, n = n, d_prior = d_prior,
                                                        add_noise = add_noise, noise_level = noise_level, return_noise_free = False,
                                                        add_realdat = add_realdat, real_dat = real_dat, cell_type=cell_type)
            elif generate_sim_method == "scaden":
                simudata = generate_simulated_data_scaden(sc_data=necessary_data, samplenum = samplenum, sparse=sparse, n = n,
                                                          id_col=id_col, cell_type=cell_type)  
            else:
                raise Exception('We only provide two approach to generate simulated training data: scaden and tape.')  

        elif postfix == 'h5ad':
            # simudata = anndata.read_h5ad(necessary_data)
            # simudata = generate_simulated_data(sc_data=necessary_data, samplenum=5000, sparse=sparse, cell_type=cell_type)
            if generate_sim_method == "tape":
                simudata = generate_simulated_data_tape(sc_data=necessary_data, samplenum=samplenum, sparse=sparse, n = n, d_prior = d_prior,
                                                        add_noise = add_noise, noise_level = noise_level, return_noise_free = False,
                                                        add_realdat = add_realdat, real_dat = real_dat, cell_type=cell_type)
            elif generate_sim_method == "scaden":
                simudata = generate_simulated_data_scaden(sc_data=necessary_data, samplenum=samplenum, sparse=sparse, n = n, 
                                                          id_col=id_col, cell_type=cell_type)
            else:
                raise Exception('We only provide two approach to generate simulated training data: scaden and tape.')   

        elif postfix == 'pth':
            raise Exception('Do not accept a model as input')
        else:
            raise Exception('Please give the correct input')
    else:
        if generate_sim_method == "tape":
            if type(necessary_data) is pd.DataFrame:
                simudata = generate_simulated_data_tape(sc_data=necessary_data, samplenum=samplenum, sparse=sparse, n = n, d_prior = d_prior,
                                                        add_noise = add_noise, noise_level = noise_level, return_noise_free = False,
                                                        add_realdat = add_realdat, real_dat = real_dat, cell_type=cell_type)
            elif type(necessary_data) is anndata.AnnData:
                simudata = generate_simulated_data_tape(sc_data=necessary_data, samplenum=samplenum, sparse=sparse, n = n, d_prior = d_prior,
                                                        add_noise = add_noise, noise_level = noise_level, return_noise_free = False,
                                                        add_realdat = add_realdat, real_dat = real_dat, cell_type=cell_type)
        else:
            simudata = generate_simulated_data_scaden(sc_data=necessary_data, samplenum=samplenum, sparse=sparse, n = n, 
                                                      id_col=id_col, cell_type=cell_type)

    train_x, train_y, test_x, genename, celltypes, samplename = \
        ProcessInputData(simudata, real_bulk, sep=sep, datatype='counts', cut_variance = cut_variance, markergene = markergene, quan_norm = quan_norm) # TPM counts
    print('training data shape is ', train_x.shape, '\ntest data shape is ', test_x.shape)
    pred = test_scaden(train_x,train_y,test_x,batch_size=batch_size,epochs=epochs)
    pred = pd.DataFrame(pred, columns=celltypes, index=samplename)
    return pred
            
    # test_x = ProcessInputData(simudata, real_bulk, sep=sep, datatype='counts', cut_variance = cut_variance, markergene = markergene)
    # return test_x

def test_scaden(train_x,train_y,test_x,batch_size=128,epochs=128):
    architectures = {'m256': ([256, 128, 64, 32], [0, 0, 0, 0]),
                     'm512': ([512, 256, 128, 64], [0, 0.3, 0.2, 0.1]),
                     'm1024': ([1024, 512, 256, 128], [0, 0.6, 0.3, 0.1])}
    model = scaden(architectures, train_x, train_y, batch_size=batch_size, epochs=epochs)
    model.train()
    pred = model.predict(test_x)
    return pred