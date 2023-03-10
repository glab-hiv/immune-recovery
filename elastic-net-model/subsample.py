import os
import numpy as np
from tqdm import tqdm
import scanpy as sc
from utils import *
import sys
from kh import *

directory = 'data'

print('reading in data')
sample_id = 'A5248' #'LT-ART'
adata = sc.read_h5ad(os.path.join(directory, f'{sample_id}_preprocessed.h5ad'))    
       
# KH subsampling: please see repo: https://github.com/CompCy-lab/SketchKH
fcs_files = list(adata.obs['fcs_file'].cat.categories)

gamma = 1.
num_samples = 2500
num_sample_sets = len(fcs_files)
subsample_idx = []
for i in tqdm(range(num_sample_sets)):
    sys.stdout.write(str(i) + '\n')
    fcs_file = adata.obs['fcs_file'].cat.categories[i]
    fcs_data = adata[adata.obs['fcs_file'] == fcs_file]
    phi = random_feats(fcs_data.X, gamma=gamma, frequency_seed=0)
    kh_indices = kernel_herding(phi, num_samples)
    subsample_idx.append(kh_indices)

np.save(os.path.join(directory, f'{sample_id}_kh_subsample_{num_samples}.npy'), np.asarray(subsample_idx))