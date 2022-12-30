from utils import *
import os

markers_dict = {'89Y_CD45':'CD45', 
                '115In_CD3':'CD3',
                '141Pr_CD196_CCR6':'CCR6',
                '142Nd_CD19_and_Bead2':'CD19',
                '169Tm_CD25': 'CD25',
                '143Nd_CD45RA':'CD45RA',
                '144Nd_CD31':'CD31',
                '145Nd_CD4':'CD4',
                '146Nd_CD11c':'CD11c',
                '147Sm_CD99':'CD99',
                '149Sm_CD194_CCR4':'CCR4',
                '151Eu_CD14_and_Bead3':'CD14',
                '152Sm_Bcl_2':'Bcl-2',
                '153Eu_CD185_CXCR5_and_Bead4':'CXCR5',
                '155Gd_CD279_PD_1':'PD-1',
                '156Gd_CD183_CXCR3':'CXCR3',
                '158Gd_CD27':'CD27',
                '159Tb_FoxP3':'FoxP3',
                '160Gd_CD28':'CD28',
                '162Dy_Ki_67':'Ki-67',
                '164Dy_CD95_FAS':'CD95',
                '165Ho_CD127_and_Bead5':'CD127',
                '166Er_CD57':'CD57',
                '167Er_CD197_CCR7':'CCR7',
                '168Er_CD8':'CD8',
                '170Er_HLA_DR':'HLA-DR',
                '171Yb_CD195_CCR5':'CCR5',
                '172Yb_CD38':'CD38',
                '173Yb_TCRgd':'TCRgd',
                '175Lu_CCR10_and_Bead6':'CCR10',
                '176Yb_CD56_and_Bead7':'CD56',
                '209Bi_CD16':'CD16'}

#create csv of metadata files
data_directory = 'data'
sample_folder = 'A5248'

dir = os.path.join(data_directory, sample_folder)
metadata_A5248 = obtain_metadata(dir)
metadata_A5248.to_csv(os.path.join(data_directory, '{}_metadata.csv'.format(sample_folder)))
adata_A5248 = run_preprocess(data_directory = data_directory, sample_folder = sample_folder, markers_dict = markers_dict)

sample_folder = 'LT-ART'
dir = os.path.join(data_directory, sample_folder)
metadata_LT = obtain_metadata(dir)
metadata_LT.to_csv(os.path.join(data_directory, '{}_metadata.csv'.format(sample_folder)))
adata_LT = run_preprocess(data_directory = data_directory, sample_folder = sample_folder, markers_dict = markers_dict)