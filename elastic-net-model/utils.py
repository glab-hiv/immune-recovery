import flowkit as fk
import os
import glob
import anndata
import re
import numpy as np
import pandas as pd
import phenograph
from sklearn.cluster import *
from sklearn.metrics import calinski_harabasz_score, silhouette_score, davies_bouldin_score
from warnings import warn
import seaborn as sns
from sklearn.preprocessing import StandardScaler
import matplotlib.pyplot as plt
from sklearn.linear_model import ElasticNet
from sklearn.model_selection import GridSearchCV, LeaveOneGroupOut
from sklearn.metrics import r2_score, mean_squared_error
from tqdm import tqdm_notebook, tqdm
from IPython import get_ipython
import logging

def obtain_metadata(dir):
    """
    Parses metadata into dataframe.
    Parameters
    ----------
    dir: str
        string referring to the directory where metadata csv files are found
    Returns
    ----------
    metadata: pandas.DataFrame
        dataframe containing fcs_file id, patient id, and timepoint information
    """
    files = glob.glob(os.path.join(dir, "*"))
    fcs_file = [i.split('\\')[-1].split('.fcs')[0] for i in files]
    id = [str(i.split('_')[1]) for i in fcs_file]
    timepoint = [str(re.findall('\d+', i.split('_')[2])[0]) for i in fcs_file]

    metadata = pd.DataFrame({'fcs_file': fcs_file,
                            'patient_id': id,
                            'timepoint': timepoint})
    return metadata

def preprocess_fcs(directory = None,
                    fcs_filename = None,
                    markers_dict = None,
                    cofactor = 5):
    """
    Arc-sinh transform data from FCS file of interest
    Parameters
    ----------
    directory: str
        string referring to the directory where fcs files are found
    fcs_filename: str
        string referring to the fcs filename of interest
    markers_dict: dict
        dictonary containing marker information and relabeling
    cofactor: int:
        integer for arc-sinh transformation
    Returns
    ----------
    df: pandas.DataFrame
        arc-sinh transformed data
    """
    sample = fk.Sample(os.path.join(directory, fcs_filename + '.fcs'), ignore_offset_error = True) #214346 events, 71 channels
    events = sample._get_raw_events()
    channel_names = np.array(sample.pns_labels)  
    markers = list(markers_dict.keys())
    idx = np.where(np.isin(channel_names, markers) == True)[0]
    events_trans = np.arcsinh(1./cofactor * events[:,idx])
    channel_names = channel_names[idx]

    df = pd.DataFrame(events_trans, columns = channel_names)
    df.rename(columns = markers_dict, inplace = True) 

    return df

def run_preprocess(data_directory = None, sample_folder = None, markers_dict = None, cofactor = 5):
    """
    Preprocessing data from a series of fcs files
    Parameters
    ----------
    data_directory: str
        string referring to the directory of fcs folder containing fcs files
    sample_folder: str
        string referring to the folder containing fcs files
    markers_dict: dict
        dictonary containing marker information and relabeling
    cofactor: int:
        integer for arc-sinh transformation
    Returns
    ----------
    df: pandas.DataFrame
        arc-sinh transformed data
    """
    metadata_all = pd.read_csv(os.path.join(data_directory, sample_folder+'_metadata.csv'), index_col = 0)
    
    dir = os.path.join(data_directory, sample_folder)
    file_list = glob.glob(os.path.join(dir, "*"))
    adata_ = []
    for file in file_list:
        fcs_filename = file.split('\\')[-1].split('.fcs')[0]
        df = preprocess_fcs(directory = dir, fcs_filename = fcs_filename, markers_dict = markers_dict, cofactor = cofactor)
        adata = anndata.AnnData(df)
        metadata = metadata_all[metadata_all['fcs_file'] == fcs_filename]
        adata.obs = pd.DataFrame(np.tile(metadata.values, (len(adata.obs),1)), columns = metadata.columns).astype('str')
        adata_.append(adata)

    adata = anndata.concat(adata_)
    adata.obs_names = [str(i) for i in range(0, len(adata.obs_names))]
    adata.write(os.path.join(data_directory, sample_folder+"_preprocessed.h5ad"))

    return adata

def phenograph_clustering(data: pd.DataFrame,
                          features: list,
                          verbose: bool,
                          print_performance_metrics: bool = True,
                          **kwargs):
    """
    Perform clustering of single cell data using PhenoGraph algorithm (https://github.com/dpeerlab/PhenoGraph)
    Code accessed and modified from flow and mass cytometry analysis toolkit CytoPy: https://github.com/burtonrj/CytoPy
    
    Parameters
    ----------
    data: pandas.DataFrame
        dataframe with dimensions samples x markers
    features: list
        columns to peform clustering on (e.g. phenotypic markers)
    verbose: bool
        if True, provides a progress bar when global_clustering is False
    print_performance_metrics: bool = True
        Print Calinski-Harabasz Index, Silhouette Coefficient, and Davies-Bouldin Index
        (see https://scikit-learn.org/stable/modules/clustering.html#clustering-performance-evaluation)
    kwargs:
        Additional keyword arguments passed when calling phenograph.cluster
    Returns
    -------
    Pandas.DataFrame, scipy.sparse.base.spmatrix, float
        modified dataframe with clustering IDs assigned to the column 'cluster_label', sparse graph
        matrix, and modularity score for communities (Q)
    """
    _print = vprint(verbose=verbose)
    data["cluster_label"] = None
    graphs = dict()
    q = dict()
    for _id, df in data.groupby("sample_id"):
        _print(f"----- Clustering {_id} -----")
        communities, graph, q_ = phenograph.cluster(df[features], clustering_algo = 'leiden', **kwargs)
        graphs[_id], q[_id] = graph, q_
        df["cluster_label"] = communities
        data.loc[df.index, ["cluster_label"]] = df.cluster_label
        if print_performance_metrics:
            clustering_performance(df[features], df["cluster_label"].values)
        _print("-----------------------------")
        _print("\n")
    return data, graphs, q

def phenograph_metaclustering(data: pd.DataFrame,
                              features: list,
                              verbose: bool = True,
                              print_performance_metrics: bool = True,
                              **kwargs):
    """
    Perform meta-clustering of single cell data using PhenoGraph algorithm (https://github.com/dpeerlab/PhenoGraph)
    Code accessed and modified from flow and mass cytometry analysis toolkit CytoPy: https://github.com/burtonrj/CytoPy
    Parameters
    ----------
    data: pandas.DataFrame
        clustered data with columns for sample_id and cluster_label
    features: list
       columns to peform clustering on (e.g. phenotypic markers)
    print_performance_metrics: bool = True
        print Calinski-Harabasz Index, Silhouette Coefficient, and Davies-Bouldin Index
        (see https://scikit-learn.org/stable/modules/clustering.html#clustering-performance-evaluation)
    verbose: bool (default=True)
        whether to provide feedback to stdout
    kwargs:
        keyword arguments passed to phenograph.cluster
    Returns
    -------
    pandas.DataFrame
        Updated dataframe with a new column named 'meta_label' with the meta-clustering
        associations
    """
    vprint_ = vprint(verbose)
    vprint_("----- Phenograph meta-clustering ------")
    metadata = summarise_clusters(data, features, None, None, 'median')
    vprint_("...summarising clusters")
    vprint_("...clustering the clusters")
    communities, graph, q = phenograph.cluster(metadata[features].values, clustering_algo = 'leiden', **kwargs)
    metadata["meta_label"] = None
    metadata["meta_label"] = communities
    if print_performance_metrics:
        clustering_performance(metadata[features], metadata["meta_label"].values)
    vprint_("...assigning meta-labels")
    data = _assign_metalabels(data, metadata)
    vprint_("------ Complete ------")
    return data, graph, q

def summarise_clusters(data: pd.DataFrame,
                       features: list,
                       scale: str or None = None,
                       scale_kwargs: dict or None = None,
                       summary_method: str = "median"):
    """
    Code accessed and modified from flow and mass cytometry analysis toolkit CytoPy: https://github.com/burtonrj/CytoPy
    Average cluster parameters along columns average to generated a centroid for
    meta-clustering
    Parameters
    ----------
    data: Pandas.DataFrame
        Clustering results to average
    features: list
        List of features to use when generating centroid
    summary_method: str (default='median')
        Average method, should be mean or median
    scale: str, optional
        Perform scaling of centroids; see cytopy.transform.Scaler
    scale_kwargs: dict, optional
        Additional keyword arguments passed to Scaler
    Returns
    -------
    Pandas.DataFrame
    Raises
    ------
    ValueError
        If invalid method provided
    """
    if summary_method == "median":
        data = data.groupby(["sample_id", "cluster_label"])[features].median().fillna(0).reset_index()
    elif summary_method == "mean":
        data = data.groupby(["sample_id", "cluster_label"])[features].mean().fillna(0).reset_index()
    else:
        raise ValueError("summary_method should be 'mean' or 'median'")
    scale_kwargs = scale_kwargs or {}
    if scale is not None:
        scaler = Scaler(method=scale, **scale_kwargs)
        data = scaler(data=data, features=features)
    return data


def clustering_performance(data: pd.DataFrame,
                           labels: list):
    'Code accessed and modified from flow and mass cytometry analysis toolkit CytoPy: https://github.com/burtonrj/CytoPy'
    print("Clustering performance...")
    print(f"Silhouette coefficient: {silhouette_score(data.values, labels, metric='euclidean')}")
    print(f"Calinski-Harabasz index: {calinski_harabasz_score(data.values, labels)}")
    print(f"Davies-Bouldin index: {davies_bouldin_score(data.values, labels)}")

def frequency_feats(data = None, meta_label = 'meta_label'):    
    """
    Computes frequency features (proportion of cells assigned to a metacluster within a sample)

    Parameters
    ----------
    data: pandas.DataFrame
    meta_label: str
        string referencing column containing metacluster label ID

    Returns
    ----------
    freq_feats: pandas.DataFrame
        dataframe containing frequency features for each sample
    """
    num = data.groupby(['sample_id', meta_label]).size().unstack().fillna(0)
    denom = num.sum(1)
    freq_feats = num.divide(denom, axis = 0)
    freq_feats.columns = ['freq_{}'.format(i) for i in freq_feats.columns]

    return freq_feats

def functional_feats(data = None, meta_label = 'meta_label', functional_markers = None):
    """
    Computes functional features (median expression of metaclusters within a sample)

    Parameters
    ----------
    data: pandas.DataFrame
    meta_label: str
        string referencing column containing metacluster label ID
    functional_markers: list
        list containing protein markers that should be used for computing functional markers (e.g. median expression per cluster)
    Returns
    ----------
    func_feats: pandas.DataFrame
        dataframe containing functional features for each sample
    """
    func_feats = data.groupby(['sample_id', meta_label])[functional_markers].median().unstack().fillna(0)
    func_feats.columns = ['{}_exp_{}'.format(i, j) for i, j in func_feats.columns]
    return func_feats

def delta_abundance(data = None, meta_label = 'meta_label', time_label = None, patient_label = None):    
    """
    Computes change in total abundance of cells with respect to time 0 for all samples for all metaclusters

    Parameters
    ----------
    data: pandas.DataFrame
    meta_label: str
        string referencing column containing metacluster label ID
    time_label: str
        string referencing timepoint label: e.g. ACTG: 'day', CHI: 'month'
    patient_label: str
        string referencing column containing patient label IDs
    Returns
    ----------
    func_feats: pandas.DataFrame
        dataframe containing functional features for each sample
    """
    abund_df = data.groupby(['sample_id', meta_label]).size().unstack().fillna(0)
    abund_df.columns = abund_df.columns.astype('str')
    abund_df[patient_label] = [i.split('_')[1].split('_')[0] for i in abund_df.index]

    delta_abund_feats = []
    for i, j in abund_df.groupby(patient_label):
        j = j.drop(columns = patient_label)
        d0_ix = j.index[j.index.str.contains(f'{time_label}0')][0]
        delta_abund_feats.append(j.sub(j.loc[d0_ix]))

    delta_abund_feats = pd.concat(delta_abund_feats)
    delta_abund_feats.columns = ['diff_abund_{}'.format(i) for i in delta_abund_feats.columns]

    return delta_abund_feats

def delta_functional_feats(data = None, meta_label = 'meta_label', functional_markers = None,  time_label = None, patient_label = None):
    """
    Computes change in functional marker expression of a metacluster with respect to time 0 for all samples for all metaclusters

    Parameters
    ----------
    data: pandas.DataFrame
    meta_label: str
        string referencing column containing metacluster label ID
    functional_markers: list
        list containing protein markers that should be used for computing functional markers (e.g. median expression per cluster)
    time_label: str
        string referencing timepoint label: e.g. ACTG: 'day', CHI: 'month'
    patient_label: str
        string referencing column containing patient label IDs
    Returns
    ----------
    func_feats: pandas.DataFrame
        dataframe containing functional features for each sample
    """
    functional_sub_df = functional_feats(data = data, meta_label =meta_label, functional_markers = functional_markers)
    functional_sub_df[patient_label] = [i.split('_')[1].split('_')[0] for i in functional_sub_df.index]

    functional_sub_feats = []
    for i, j in functional_sub_df.groupby(patient_label):
        j = j.drop(columns = patient_label)
        d0_ix = j.index[j.index.str.contains(f'{time_label}0')][0]
        functional_sub_feats.append(j.sub(j.loc[d0_ix]))

    functional_sub_feats = pd.concat(functional_sub_feats)
    functional_sub_feats.columns = ['diff_{}'.format(i) for i in functional_sub_feats.columns]

    return functional_sub_feats

def engineer_features(data = None, meta_label = 'meta_label', functional_markers = None, time_label = None, patient_label = None):
    """
    Derives frequency, functional, abundance, and change in abundance features from metacluster matrix  

    Parameters
    ----------
    data: pandas.DataFrame
    meta_label: str
        string referencing column containing metacluster label ID
    functional_markers: list
        list containing protein markers that should be used for computing functional markers (e.g. median expression per cluster)
    time_label: str
        string referencing timepoint label: e.g. ACTG: 'day', CHI: 'month'
    patient_label: str
        string referencing column containing patient label IDs

    Returns
    ----------
    engineered_feats: pandas.DataFrame
        dataframe containing engineered features for each sample: dimensions: samples x engineered features 
    scores: pandas.DataFrame
        dataframe containing the R^2 and MSE scores across LOO cross validation folds
    """
    freq_feats = frequency_feats(data, meta_label = meta_label)
    func_feats = functional_feats(data, meta_label = meta_label, functional_markers = functional_markers)

    delta_abund_feats = delta_abundance(data, meta_label =meta_label, time_label = time_label, patient_label = patient_label)
    delta_func_feats = delta_functional_feats(data, meta_label = meta_label, functional_markers = functional_markers,  time_label = time_label, patient_label = patient_label)

    engineered_feats = pd.concat([freq_feats, func_feats, delta_abund_feats, delta_func_feats], axis = 1)
    engineered_feats = engineered_feats.loc[:, ~(engineered_feats==0).all(axis=0)]

    return engineered_feats

def elastic_net(df, timepoint_label):
    """
    Performs elastic net regression on metacluster-derived immune features 
    Parameters
    ----------
    df: pandas.DataFrame
        dataframe containing engineered immune features: dimensions samples x derived immune features 
    timepoint_label: str
        string referencing timepoint label for response vector. For ACTG: 'day, For CHI: 'month

    Returns
    ----------
    prediction_coef: pandas.DataFrame
        dataframe containing the coefficients for each feature across LOO cross validation folds
    scores: pandas.DataFrame
        dataframe containing the R^2 and MSE scores across LOO cross validation folds
    """
    scores = pd.DataFrame()
    prediction_coef = pd.DataFrame()
    
    param_grid = {'l1_ratio': np.arange(0, 1, 0.01), 'alpha': [1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 0.0, 1.0, 10.0, 100.0]}

    scaler = StandardScaler()
    X = scaler.fit_transform(df)

    y = np.asarray([int(i.split(timepoint_label)[1]) for i in df.index])
    patients = np.asarray([i.split('_')[1].split('_')[0] for i in df.index])
    feature_names = df.columns

    logo = LeaveOneGroupOut()

    for train_ix, test_ix in logo.split(X, y, patients):
        X_train, X_test = X[train_ix, :], X[test_ix, :]
        y_train, y_test = y[train_ix], y[test_ix]
        patients_train, patients_test = patients[train_ix], patients[test_ix]
        model = ElasticNet(max_iter = 2000)
        gsearch = GridSearchCV(model, param_grid, scoring = 'neg_mean_absolute_error', n_jobs = -1, cv = logo, refit = True)
        
        result = gsearch.fit(X_train, y_train, groups = patients_train)
        best_model = result.best_estimator_
        y_pred = best_model.predict(X_test)

        scores_ = pd.DataFrame({'R2': [r2_score(y_test, y_pred)],
                                'RMSE': [np.sqrt(mean_squared_error(y_test, y_pred))]})
        scores = pd.concat([scores_, scores], axis = 0)
        prediction_coef_ = pd.DataFrame(best_model.coef_, index = feature_names)
        prediction_coef = pd.concat([prediction_coef_, prediction_coef], axis= 1)

    return prediction_coef, scores

def plot_features(feature_list = None, feats = None, hue = None, filename_save = None, time = None, colors = None):
    df_lineplot = feats.copy()
    df_lineplot['timepoint'] = [int(i.split(time)[1]) for i in df_lineplot.index]
    df_lineplot['patient'] = [i.split('_')[1].split('_')[0] for i in df_lineplot.index]

    sns.set_style('ticks')
    fig, axes = plt.subplots(1, 6, figsize = (28,3.75), gridspec_kw={'hspace': 0.3, 'bottom':0.15, 'wspace':0.35}) #4,5 #26,16
    for i, ax in zip(range(len(feature_list)), axes.flat):
        if hue is None:
            g = sns.lineplot(y = feature_list[i], x = 'timepoint', data = df_lineplot, err_style='band', ax = ax, color = '#675068', marker = 'o',ci = 'sd')
            ax.tick_params(labelsize=16)
            ax.set_title(feature_list[i], fontsize = 14)
            ax.set_xlabel(time, fontsize = 16)
            ax.set_ylabel('', fontsize = 16)
        else:
            g = sns.lineplot(y = feature_list[i], x = 'timepoint', data = df_lineplot, hue = hue, ax = ax, marker = 'o', palette = colors, ci = 'sd')
            ax.tick_params(labelsize=16)
            ax.set_title(feature_list[i], fontsize = 14)
            ax.set_xlabel(time, fontsize = 16)
            ax.set_ylabel('', fontsize = 16)
            if (i == 4) or (i == 9):
                ax.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0., title = hue, prop={'size':12}, title_fontsize=14)
            else:
                ax.legend([],[], frameon=False)

        plt.savefig(filename_save+'_features.pdf', bbox_inches="tight", dpi = 150)

def plot_R2(scores, filename):
    sns.set_style('ticks')
    fig, axes = plt.subplots(1, 1, figsize = (2.5,3), gridspec_kw={'hspace': 0.1, 'bottom':0.15, 'wspace':0.1})
    sns.boxplot(data = scores['R2'].values, width=0.4, dodge=True, color='#675068', fliersize=0, linewidth = 1, ax = axes, boxprops=dict(alpha=0.6))
    sns.stripplot(data = scores['R2'].values, color='#675068', linewidth = 0.8, size=5, edgecolor="black", split=True, jitter = True, dodge = False)
    axes.tick_params(labelsize=16)
    axes.set(xticklabels=[])
    axes.set_ylabel('$R^2$', fontsize = 16)
    axes.set_ylim(0, 0.85)
    plt.savefig(filename+'_R2_score.pdf', bbox_inches="tight")

def setup_standard_logger(name: str,
                          default_level: int or None = None,
                          log: str or None = None) -> logging.Logger:
    """
    Code accessed from flow and mass cytometry analysis toolkit CytoPy: https://github.com/burtonrj/CytoPy
    Convenience function for setting up logging.
    Parameters
    ----------
    name: str
        Name of the logger
    default_level: int
        Default level at which data is logged (defaults to INFO)
    log: str, optional
        Optional filepath to print logs too, if not provided, logging is printed to stdout
    Returns
    -------
    logging.Logger
    """
    default_level = default_level or logging.INFO
    logger = logging.getLogger(name)
    logger.setLevel(default_level)
    if log is not None:
        handler = logging.FileHandler(filename=log)
    else:
        handler = logging.StreamHandler()
    handler.setFormatter(logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s'))
    logger.addHandler(handler)
    return logger


def progress_bar(x: iter,
                 verbose: bool = True,
                 **kwargs) -> callable:
    """
    Code accessed from flow and mass cytometry analysis toolkit CytoPy: https://github.com/burtonrj/CytoPy

    Generate a progress bar using the tqdm library. If execution environment is Jupyter, return tqdm_notebook
    otherwise used tqdm.
    Parameters
    -----------
    x: iterable
        some iterable to pass to tqdm function
    verbose: bool, (default=True)
        Provide feedback (if False, no progress bar produced)
    kwargs:
        additional keyword arguments for tqdm
    :return: tqdm or tqdm_notebook, depending on environment
    """
    if not verbose:
        return x
    if which_environment() == 'jupyter':
        return tqdm_notebook(x, **kwargs)
    return tqdm(x, **kwargs)


def which_environment() -> str:
    """
    Code accessed from flow and mass cytometry analysis toolkit CytoPy: https://github.com/burtonrj/CytoPy

    Test if module is being executed in the Jupyter environment.
    Returns
    -------
    str
        'jupyter', 'ipython' or 'terminal'
    """
    try:
        ipy_str = str(type(get_ipython()))
        if 'zmqshell' in ipy_str:
            return 'jupyter'
        if 'terminal' in ipy_str:
            return 'ipython'
    except:
        return 'terminal'


def vprint(verbose: bool):
    """
    Code accessed from flow and mass cytometry analysis toolkit CytoPy: https://github.com/burtonrj/CytoPy 

    Utility function for optional printing.
    Parameters
    ----------
    verbose: bool
        If True, returns print function, else False
    Returns
    -------
    callable
    """
    return print if verbose else lambda *a, **k: None

def sample_n(data: pd.DataFrame,
              sample_label: str):
    'Code accessed from flow and mass cytometry analysis toolkit CytoPy: https://github.com/burtonrj/CytoPy'
    sample_size = data[sample_label].value_counts()
    sample_size.name = "sample_n"
    return pd.DataFrame(sample_size).reset_index().rename({"index": sample_label}, axis=1)
    
def cluster_n(data: pd.DataFrame,
               cluster_label: str,
               sample_label: str):
    'Code accessed from flow and mass cytometry analysis toolkit CytoPy: https://github.com/burtonrj/CytoPy'
    sample_cluster_counts = data.groupby(sample_label)[cluster_label].value_counts()
    sample_cluster_counts.name = "cluster_n"
    return pd.DataFrame(sample_cluster_counts).reset_index()


def cluster_size(sample_n: pd.DataFrame,
                  cluster_n: pd.DataFrame):
    'Code accessed from flow and mass cytometry analysis toolkit CytoPy: https://github.com/burtonrj/CytoPy'
    cluster_size = cluster_n.merge(sample_n, on="sample_id")
    cluster_size["cluster_size"] = cluster_size["cluster_n"] / cluster_size["sample_n"]
    return cluster_size


def label_centroids(data: pd.DataFrame,
                     centroids: pd.DataFrame,
                     sample_label: str,
                     cluster_label: str,
                     target_label: str):
    'Code accessed from flow and mass cytometry analysis toolkit CytoPy: https://github.com/burtonrj/CytoPy'
    data = data[[sample_label, cluster_label, target_label]].drop_duplicates()
    return centroids.merge(data, on=[sample_label, cluster_label])

def _assign_metalabels(data: pd.DataFrame,
                       metadata: pd.DataFrame):
    """
    Code accessed and modified from flow and mass cytometry analysis toolkit CytoPy: https://github.com/burtonrj/CytoPy
    Assign the meta-cluster labels to the original data and return the modified dataframe with the meta cluster
    labels in a new column called 'meta_label'
    Parameters
    ----------
    data: Pandas.DataFrame
    metadata: Pandas.DataFrame
    Returns
    -------
    Pandas.DataFrame
    """
    return pd.merge(data, metadata[["sample_id", "cluster_label", "meta_label"]], on=["sample_id", "cluster_label"])

