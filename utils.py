from ast import arg
import os
import csv
import numpy as np
import pandas as pd
import pathlib
#from numba import jit


# special symbols
DIR_SEP = '/'

# paths
# assuming the structure PROJ_DIR/dataproc/utils.py
PROJ_DIR = str(pathlib.Path(__file__).parent.resolve())
#PROJ_DIR = os.path.join(PROJ_DIR, 'GRiPE')
#PROJ_DIR = os.path.abspath(os.pardir) 
DATAPROC_DIR = 'dataproc'
BIODATA_DIR = 'biodata'
# default directory with results
RES_DIR  = 'results'
FIG_DIR = os.path.join(DATAPROC_DIR, 'figures')
ID_DIR = os.path.join(DATAPROC_DIR, 'experiment_id_lists')
JAVA_CLASS = 'simulator.SimulatorCLI'
JAVA_CLASS_PATH = os.path.join(PROJ_DIR, 'out', 'production', 'GRiPE')
SITE_OCCUPANCY_SCRIPT = 'site_occupancy.py'

# parts of file names
OCCUPANCY = '_occupancy_'
PARAMS = '_params_'
AFFINITY = '_affinity_landscape_'
TF_SPECIES = '_TF_species_'
STATUS = '_status_'
REPRESSED_LEN = '_repressed_lengths_'
TARGET_SITE_FOLLOW = '_target_site_follow_'

# file fields
SPEC_WAITING_TIME = 'SPECIFICWAITINGTIME'
TOTAL_SIM_TIME = 'STOP_TIME'
ENSEMBLE_SIZE = 'ENSEMBLE_SIZE'
TS_FILE = 'TS_FILE'
BTRACK = 'DNA_AVAILABILITY_FILE'
ASSOC_RATE = 'ASSOCRATE'
UNBIND_PROB = 'UNBINDINGPROBABILITY'
JUMP_PROB = 'JUMPINGPROBABILITY'
SIZE_LEFT = 'SIZELEFT'
SIZE_RIGHT = 'SIZERIGHT'
MOTIF = 'DBD'
PWM = 'PWM:'
REPRESSION_RATE = 'REPRESSIONRATE'
DEREPRESSION_ATTENUATION_FACTOR = 'DEREPRESSIONATTENUATIONFACTOR'
COPY_NUMBER = 'COPYNUMBER'
SPEC_ENERGY_THRESH = 'SPECIFICENERGYTHRESHOLD'
SCORE_TO_ENERGY = 'ES'
OUTPUT_REPRESSED_LENGTHS = 'OUTPUT_REPRESSED_LENGTHS'
DEBUG_MODE = 'DEBUG_MODE'
DNA_DEREPRESSION_RATE = 'DNA_DEREPRESSION_RATE'
REPR_LEN_LEFT = 'REPRLENLEFT'
REPR_LEN_RIGHT = 'REPRLENRIGHT'
OUTPUT_TF_POINTS = 'OUTPUT_TF_POINTS'
FOLLOW_TS = 'FOLLOW_TS'
OUTPUT_FOLDER = 'OUTPUT_FOLDER'
OUTPUT_FILENAME = 'OUTPUT_FILENAME'

# output file names
PLOT_PROB_SITE = 'sites_probabilities_'
PLOT_PROB_ENERGY = 'energy_probabilities_'

# desired types of dataframe columns
TF_FIELDS_TYPES_DICT = {REPRESSION_RATE: float, 
                        SPEC_WAITING_TIME: float, 
                        SPEC_ENERGY_THRESH: float, 
                        SCORE_TO_ENERGY: float, 
                        COPY_NUMBER: int, 
                        DEREPRESSION_ATTENUATION_FACTOR: float,
                        SIZE_LEFT: int, SIZE_RIGHT: int}

# order of TFs
TF_ORDER_LIST_DROSOPHILA = ['Kr', 'Hb', 'Gt', 'Kni', 'Bcd', 'Cad', 'Tll', 'Hkb']



# memoizer
def memoize(f):
    results = {}
    def helper(*args):
        if args not in results:
            results[args] = f(*args)
        return results[args]
    return helper 

#
# TF-DNA functions
#
def get_affinity_from_df(df_site_affinity, tf_size, btrack, affinity_thresh=None):
    affinity = df_site_affinity.to_numpy().transpose()[1:]
    for tf_i in range(len(tf_size)):
        thresh = -np.inf
        if isinstance(affinity_thresh, (list, np.ndarray)):
            thresh = affinity_thresh[tf_i]
        for bp_i in range(affinity.shape[1] - tf_size[tf_i] + 1):
            #if (not btrack[bp_i]) or (not btrack[bp_i + tf_size[tf_i] - 1]):
            #    affinity[tf_i, bp_i: bp_i+tf_size[tf_i]-1] = -np.inf
            if not is_opened(bp_i, tf_size[tf_i], btrack):
                affinity[tf_i, bp_i] = -np.inf
            affinity[tf_i, bp_i] = max(affinity[tf_i, bp_i], thresh)
        # last positions which TFs cannot reach
        affinity[tf_i,-tf_size[tf_i]+1:] = -np.inf
    return affinity

def get_occupancy_from_df(df_occupancy, ensemble_size=1):
    return df_occupancy.to_numpy().transpose()[2:] / ensemble_size

#@memoize
def get_tf_size(df_tf_data):
    tf_list = np.array(df_tf_data.index, dtype=str)
    tf_size = np.array([get_motif_size(df_tf_data.at[tf, MOTIF])
                        + int(df_tf_data.at[tf, SIZE_LEFT])
                        + int(df_tf_data.at[tf, SIZE_RIGHT]) for tf in tf_list], dtype=int)
    return tf_size
    
def get_motif_size(str_pwm):
    if str_pwm == '':
        return 0
    if not str_pwm.startswith('PWM:'):
        raise ValueError('Non-PWM string is not supported.')
    # pwm
    return len(str_pwm.split(';')[0].split(','))


def is_opened(bp_i, motif_len, btrack):
    if np.sum(btrack[bp_i: bp_i+motif_len]) < motif_len:
        return False
    return True
        
def bp_to_idx(bp):
    bps = list('acgtn')
    return bps.index(bp)

def complementary(bp):
    compl_bps = list('tgcan')
    return compl_bps[bp_to_idx(bp)]

def complementary_idx(idx):
    return 3 - idx

# get dict with PWMs for each TF from dataframe
def get_pwms(df):
    d = dict()
    for tf in df.index:
        # prepare string like 'PWM:A=[-0.1, 0.4, 1.02]; C=[...'
        str_pwm = df.at[tf, MOTIF].replace(PWM, '').replace(' ', '')
        # get rows
        str_pwm_rows = str_pwm.split(';')
        pwm = []
        # parse rows
        for srow in str_pwm_rows:
            srow = srow.split('=')[-1][1:-1].split(',')
            row = [float(el) for el in srow]
            pwm.append(row)
        # add the resulting pwm to the dict
        d[tf.lower()] = np.array(pwm)
    return d

def pwm_string_to_pwm(pwm_string):
    pwm_str_nucl = pwm_string.replace('PWM: ', '').split('; ')
    pwm = []
    for str_nucl_scores in pwm_str_nucl:
        str_nucl_scores = str_nucl_scores.split(' = ')[-1]
        str_nucl_scores = str_nucl_scores.replace('[', '').replace(']', '')
        pwm.append([float(s) for s in str_nucl_scores.split(',')])
    return np.array(pwm)

#@memoize
def pwm_score(pwm, seq, norm=True):
    score = 0
    max_score = 0
    for i, bp in enumerate(seq):
        score += pwm[bp_to_idx(bp), i]
        max_score += pwm[:, i].max()
    if norm: score -= max_score
    return score
"""
@jit(nopython=True)
def pwm_score_jit(pwm, seq_idx, max_score_bp, norm=True):
    score = 0
    max_score = 0
    #max_score_bp = np.max(pwm, axis=0)
    for i in np.arange(seq_idx.shape[0]):
        score = score + pwm[seq_idx[i], i]
        max_score += max_score_bp[i]
    if norm:
        return score - max_score
    return score

@jit(nopython=True)
def pwm_score_jit(pwm, seq_idx, max_score_bp, norm=True):
    score = 0
    max_score = 0
    #max_score_bp = np.max(pwm, axis=0)
    for i in np.arange(seq_idx.shape[0]):
        score = score + pwm[seq_idx[i], i]
        max_score += max_score_bp[i]
    if norm:
        return score - max_score
    return score
"""
#
# Helpful utility functions
#
def sum_neighbours(array):
    n = len(array)
    if not n % 2 == 0:
        raise ValueError('Odd array length')
    return np.array([array[i] + array[i+1] for i in range(0, n, 2)])

def double_each_element(array):
    n = len(array)
    return np.array([[el,]*2 for el in array]).reshape(2*n)
    

#
# Results selecting and processing functions
#
#@memoize
def process_results(path, prefix, id, is_energy=False):
    # load occupancy (occupancy.wig)
    print(path, prefix, id)
    df_occupancy = read_wig_to_dataframe(get_name(path, prefix, OCCUPANCY, id, 'wig'))
    df_occupancy['collisionsCount'] = np.array(df_occupancy['collisionsCount'], dtype=int)
    # load energies (PFM-based)
    df_site_affinity = read_wig_to_dataframe(get_name(path, prefix, AFFINITY, id, 'wig'))
    if is_energy:
        df_site_affinity.iloc[:,1:] = -df_site_affinity.iloc[:,1:]
    # load tf data
    df_tf_data = read_csv_to_dataframe(get_name(path, prefix, TF_SPECIES, id, 'csv'))
    is_two_state = list(map(lambda ele: ele.strip().lower().capitalize() == "True", df_tf_data['ISTWOSTATERANDOMWALK'].to_list()))
    affinity_thresh = df_tf_data[SPEC_ENERGY_THRESH].to_numpy() / df_tf_data['ES'].to_numpy()
    for i, b in enumerate(is_two_state):
        if not b:
            affinity_thresh[i] = -np.inf
    # load params
    dict_params = read_params_to_dict(get_name(path, prefix, PARAMS, id, 'grp'))
    for key in dict_params.keys():
        dict_params[key] = dict_params[key].replace('"','') # remove ""
    # both directions or not (two strands or a single strand)
    both_dir = is_both_directions(df_tf_data, df_occupancy)
    # load target sites
    ts_filename = os.path.join(PROJ_DIR, *dict_params[TS_FILE].split(DIR_SEP)).replace('"', '')
    if os.path.isfile(ts_filename):
        df_target_sites = create_target_site_df(ts_filename, df_site_affinity, both_dir, df_tf_data)
    else:
        df_target_sites = None
    # some important numbers from parameters
    total_time = float(dict_params[TOTAL_SIM_TIME])
    ensemble_size = int(dict_params[ENSEMBLE_SIZE])
    print(f'{TOTAL_SIM_TIME}: {total_time}')
    print(f'{ENSEMBLE_SIZE}: {ensemble_size}')
    print(df_tf_data[[COPY_NUMBER, ASSOC_RATE, REPRESSION_RATE, DEREPRESSION_ATTENUATION_FACTOR]])

    positions_i = np.array(df_occupancy.index)
    positions_num = len(positions_i)
    # load positions availability (btrack)
    btrack_path = os.path.join(PROJ_DIR, *dict_params[BTRACK].split(DIR_SEP)).replace('"','')
    if os.path.isfile(btrack_path):
        btrack = np.loadtxt(btrack_path, dtype=int)
    else:
        print('No btrack file, the whole chromatin is assumed to be open')
        btrack = np.ones_like(positions_i)
    
    #total_bind_time = np.sum(df_occupancy.to_numpy(), axis=0)[2:] / ensemble_size
    total_bind_time = np.sum(get_occupancy_from_df(df_occupancy, ensemble_size), axis=1)
    print(f'Binding time of each TF on each strand: {array2str(total_bind_time, 0)}')
    # if both_dir:
    #     total_bind_time = np.sum(total_bind_time.reshape(-1,2), axis=1) # sum each pair of neighbours
    #     print(f'Binding time of each TF on both strands: {array2str(total_bind_time, 0)}')
    #     print(f'timeBoundAvg: {array2str(total_bind_time / float(dict_params[TOTAL_SIM_TIME]), 9)}')
    #     total_bind_time = double_each_element(total_bind_time)
    
    # Probability of each affinity
    tf_list = np.array(df_tf_data.index, dtype=str)
    tf_size = get_tf_size(df_tf_data)

    specific_time = np.array(df_tf_data[SPEC_WAITING_TIME], dtype=float)
    assoc_rate    = np.array(df_tf_data[ASSOC_RATE], dtype=float) / len(positions_i)
    jump_prob     = np.array(df_tf_data[JUMP_PROB], dtype=float)
    unbind_prob   = np.array(df_tf_data[UNBIND_PROB], dtype=float)
    molecule_num  = np.array(df_tf_data[COPY_NUMBER], dtype=int)
    if both_dir:
        assoc_rate    = double_each_element(assoc_rate) / 2
        specific_time = double_each_element(specific_time)
        jump_prob     = double_each_element(jump_prob)
        unbind_prob   = double_each_element(unbind_prob)
        tf_list = df_site_affinity.columns[1:]
        tf_num  = len(tf_list)
        tf_size = double_each_element(tf_size)
        molecule_num = double_each_element(molecule_num)
        affinity_thresh = double_each_element(affinity_thresh)
        
    affinity = get_affinity_from_df(df_site_affinity, tf_size, btrack, affinity_thresh)

    #affinity_exponent = np.einsum('i,ij->ij', specific_time, np.exp(affinity)) / total_time
    affinity_exponent = np.exp(affinity)
    affinity_exp_sum  = np.zeros(tf_num)
    for tf_i in range(tf_num):
        #if both_dir:
            #affinity_exp_sum[tf_i] = np.sum(np.exp(affinity[tf_i]))
        #    ind = (tf_i // 2) * 2 # or tf_i - (tf_i%2)
        #    for s in range(2):
                #sites_strand = df_target_sites.loc[df_target_sites.name_strand == tf_list[ind+s]].pos.to_numpy(dtype=int)
                #affinity_exp_sum[tf_i] += np.sum(np.exp(affinity[ind+s][sites_strand]))
        #        affinity_exp_sum[tf_i] += np.sum(np.exp(affinity[ind+s]))
        #else:
        affinity_exp_sum[tf_i] = np.sum(np.exp(affinity[tf_i]))
    print(f'Affinity exponent sum for each TF: {array2str(affinity_exp_sum, 2)}')
    #affinity_prob = np.einsum('i,ij->ij', molecule_num / (unbind_prob*jump_prob/assoc_rate/specific_time + affinity_exp_sum), affinity_exponent)
    affinity_prob = np.einsum('i,ij->ij', 1/affinity_exp_sum, affinity_exponent)
    #affinity_prob = affinity_exponent / np.sum(affinity_exponent)*2
    occupancy_from_exper = get_occupancy_from_df(df_occupancy, ensemble_size)
    #occupancy_from_exper /= total_time
    occupancy_from_exper = np.einsum('i,ij->ij', 1/total_bind_time, occupancy_from_exper)
        
    return affinity_prob, occupancy_from_exper, tf_list, df_target_sites

def is_both_directions(*args):
    """
    args should contain either df_site_affinity or df_tf_data and df_occupancy
    """
    if len(args) == 2:
        df_tf_data, df_occupancy = args
        if len(df_tf_data.index) == len(df_occupancy.columns[2:]):
            # 1 direction
            return False
        elif 2*len(df_tf_data.index) == len(df_occupancy.columns[2:]):
            # 2 directions
            return True
        raise ValueError('Error in data.')
    elif len(args) == 1:
        df_site_affinity = args
        if df_site_affinity.columns.to_list()[-1].find("3'5'") == -1:
            # 1 direction
            return False
        # 2 directions
        return True
    else:
        raise ValueError('Error in args.')


def get_params_and_ids_by_vals(param_vals, res_dirs, prefixes, ids, vals):
    selected_params = param_vals
    selected_ids = np.array(ids)
    selected_dirs = np.array(res_dirs)
    selected_pref = np.array(prefixes)
    for i, val in enumerate(vals):
        if val and not np.isnan(val):
            idx = selected_params[:,i] == val
            selected_params = selected_params[idx]
            selected_ids = selected_ids[idx]
            selected_dirs = selected_dirs[idx]
            selected_pref = selected_pref[idx]
    return selected_params, selected_dirs, selected_pref, selected_ids

def get_ids(name):
    f = open(os.path.join(PROJ_DIR, ID_DIR, name), 'r')
    lines = f.readlines() # line 0 is description string, line 1 is for parameters vals and line 2 is for ids
    first_line = 0
    for line in lines:
        if line.startswith('#'):
            first_line += 1
        else:
            break    
    param_name = lines[first_line].strip().split(' = ')[0]
    param_vals = lines[first_line].strip().split(' = ')[1].split(', ')
    ids_to_get = lines[first_line + 1].strip().split(' = ')[1].split(', ')
    return param_name, param_vals, ids_to_get

def remove_results_not_having_these_ids(dir, ids):
    for name in os.listdir(dir):
        full_name = os.path.join(dir, name)
        if os.path.isdir(full_name):
            remove_results_not_having_these_ids(full_name, ids)
        elif os.path.isfile(full_name):
            id = name.split('_')[-1]
            if not (id == name):
                id = id.split('.')[0] # remove extension
                if not (id in ids):
                    os.remove(full_name)

def get_results_by_ids(dir, ids_to_get):
    # this function is based on get_results_match_params
    # probably this is not the best way to do it, but it works 
    res_dirs = []
    prefixes = []
    ids = []
    all_ids_in_dir = []
    file_groups = {}
    for name in os.listdir(dir):
        if type(name) == bytes:
            name = name.decode("utf-8")
        full_name = os.path.join(dir, name)
        if os.path.isdir(full_name):
            vals_from_subdir = get_results_by_ids(full_name, ids_to_get)
            for i, id in enumerate(vals_from_subdir[2]):
                if id in ids_to_get:
                    res_dirs.append(vals_from_subdir[0][i])
                    prefixes.append(vals_from_subdir[1][i])
                    ids.append(id)
        elif os.path.isfile(full_name):
            id = full_name.split('_')[-1]
            if not (id == full_name):
                id = id.split('.')[0]
                if not (id in all_ids_in_dir):
                    all_ids_in_dir.append(id)
                    file_groups[id] = []
                file_groups[id].append(full_name)
    for id in file_groups:
        # set flag
        res_found = False
        # find params file and get prefix and dir name
        for file in file_groups[id]:
            file_name = os.path.basename(file)
            file_dir = os.path.dirname(file)
            if PARAMS in file_name:
                file_prefix = file_name.split(PARAMS)[0]
                res_found = True
        if res_found:
            res_dirs.append(file_dir)
            prefixes.append(file_prefix)
            ids.append(id)
    return res_dirs, prefixes, ids

def get_results_match_params(dir, check_dict_params=None, check_tf_params=None):
    res_dirs = []
    prefixes = []
    ids = []
    all_ids_in_dir = []
    file_groups = {}
    for name in os.listdir(dir):
        full_name = os.path.join(dir, name)
        if os.path.isdir(full_name):
            vals_from_subdir = get_results_match_params(full_name, check_dict_params, check_tf_params)
            res_dirs.extend(vals_from_subdir[0])
            prefixes.extend(vals_from_subdir[1])
            ids.extend(vals_from_subdir[2])
        elif os.path.isfile(full_name):
            id = full_name.split('_')[-1]
            if not (id == full_name):
                id = id.split('.')[0]
                if not (id in all_ids_in_dir):
                    all_ids_in_dir.append(id)
                    file_groups[id] = []
                file_groups[id].append(full_name)

    for id in file_groups:
        # set flag
        res_found = False
        # find params file and get prefix and dir name
        for file in file_groups[id]:
            file_name = os.path.basename(file)
            file_dir = os.path.dirname(file)
            if PARAMS in file_name:
                file_prefix = file_name.split(PARAMS)[0]
                res_found = True
        # check params file if needed
        if res_found and check_dict_params:#dict_params_desired:
            dict_params_local = read_params_to_dict(get_name(file_dir, file_prefix, PARAMS, id, 'grp'))
            res_found = check_dict_params(dict_params_local)
            #for key in dict_params_desired:
            #    if not dict_params_desired[key] == dict_params_local[key]:
            #        res_found = False
            #        break
        # open tf_file and check params if needed
        if res_found and check_tf_params:#isinstance(df_tf_data_desired, pd.DataFrame):
            tf_file_name = get_name(file_dir, file_prefix, TF_SPECIES, id, 'csv')
            if not os.path.exists(tf_file_name):
                tf_file_name = get_name(file_dir, file_prefix, TF_SPECIES + '0.0s_', id, 'csv')
                if not os.path.exists(tf_file_name):
                    print('Something wrong with the TF file in {}, id {}.'.format(file_dir, id))
                    continue
            df_tf_data_local = read_csv_to_dataframe(tf_file_name)
            res_found = check_tf_params(df_tf_data_local)
            #df_tf_data_selected = df_tf_data_local[df_tf_data_desired.columns]
            #if not df_tf_data_selected.equals(df_tf_data_desired):
            #    res_found = False
        # add dir name, prefix and id to the output lists
        if res_found:
            res_dirs.append(file_dir)
            prefixes.append(file_prefix)
            ids.append(id)
    
    return res_dirs, prefixes, ids

def array2str(array, frac_digits):
    return ', '.join(['{1:.{0}f}'.format(frac_digits, a) for a in array])

def str2list(s, sep=', '):
    return s.replace('[', '').replace(']', '').split(sep)

def get_results_match_params_as_df(dir, check_dict_params=None, check_tf_params=None):
    res_dirs, prefixes, ids = get_results_match_params(dir, check_dict_params, check_tf_params)
    return pd.DataFrame(list(zip(res_dirs, prefixes, ids)), columns=['path', 'prefix', 'id'])

#
# Reading functions
#
def create_target_site_df(ts_filename, df_site_affinity, both_dir, df_tf_data=None):
    df = pd.DataFrame(columns=['repressor', 'name', 'name_strand', 'pos', 'size', 'strand', 'affinity'])
    f = open(ts_filename)
    i = 0
    for r in f:
        ts_info = r.split(':')
        ts_name_dir = ts_name = ts_info[0]
        ts_pos = int(ts_info[2].split('..')[0])
        ts_size = int(ts_info[2].split('..')[1]) + 1 - ts_pos
        ts_dir = int(ts_info[-1])
        if both_dir:
            ts_name_dir += ("5'3'" if ts_dir == 0 else "3'5'")
        if df_tf_data is not None:
            is_site_of_repr = float(df_tf_data.at[ts_name, REPRESSION_RATE]) > 0.0
        else:
            is_site_of_repr = False
        ts_affinity = float(df_site_affinity.loc[df_site_affinity.position == ts_pos][ts_name_dir])
        df.loc[i] = [is_site_of_repr, ts_name, ts_name_dir, ts_pos, ts_size, ts_dir, ts_affinity]
        i += 1
        #ts_list.append(TargetSite(is_site_of_repr, ts_info[0], ts_pos, ts_size, ts_dir, ts_affinity))
    f.close
    return df


def get_name(experiment_dir, file_prefix, name, experiment_id, ext=''):
    if ext and not ext.startswith('.'):
        ext = '.' + ext
    experiment_id = str(experiment_id)
    full_name = os.path.join(experiment_dir, file_prefix + name + experiment_id + ext)
    if not os.path.exists(full_name):
        file_list = []
        for file in os.listdir(experiment_dir):
            if file.startswith(file_prefix + name) and file.endswith(experiment_id + ext):
                file_list.append(file)
        if len(file_list) > 1:
            print('Warning: more than one {} file (experiment id {})'.format(name, experiment_id))
        if len(file_list) == 0:
            print('No {} file is found (experiment id {})'.format(name, experiment_id))
        else:
            for name in file_list:
                if not ('0.0s' in name):
                    full_name = os.path.join(experiment_dir, name)
    return full_name

def read_wig_to_dataframe(filename, index=None, skip_top_lines=1): 
    # skip_top_lines=1 is legacy
    f = open(filename, 'r')
    for _ in range(skip_top_lines):
        f.readline() # skip top lines
    str_header = f.readline() # read header
    list_header = str_header.replace('"', '').strip().split(', ')
    if index is not None: list_header = list_header[1:]
    result = pd.read_csv(f, names=list_header, index_col=index)
    f.close()
    return result

def read_csv_to_dataframe(filename, sep=',', dict_types=None, index_col=True):
    f = open(filename, 'r')
    reader = csv.reader(f, delimiter=sep, quotechar='"', skipinitialspace=True)
    if index_col:
        list_header = next(reader)[1:]
        list_data = []
        list_index = []
        for row in reader:
            list_index.append(row[0])
            list_data.append(row[1:])
        df = pd.DataFrame(list_data, columns=list_header, index=list_index)
    else:
        df = pd.read_csv(f, names=next(reader), delimiter=sep)
    f.close()
    if isinstance(dict_types, dict):
        df = df.astype(dict_types)
    return df

def read_params_to_dict(filename):
    f = open(filename, 'r')
    d = {}
    for row in f:
        row = row.strip()
        if row and not row.startswith('#'):
            l = row.replace(';', '').split('=')
            d.update({l[0].strip(): l[1].strip()})
    f.close()
    return d

def write_dict_to_params(filename, d):
    f = open(filename, 'w')
    for key in d.keys():
        line = str(key) + ' = ' + str(d[key]) + ';\n'
        f.write(line)
    f.close()

def read_ini_to_dict(filename):
    f = open(filename, 'r')
    d = {}
    name = ''
    for row in f:
        if row and row.startswith('name'):
            l = row.replace(';', '').split('=')
            name = l[1].replace('"', '').strip()
            d.update({name: 0})
        if row and row.startswith('value'):
            if not (name in d):
                raise ValueError('Error in .ini file.')
            l = row.replace(';', '').split('=')
            val = l[1].strip()
            d.update({name: val})
    f.close()
    return d