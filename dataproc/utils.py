import os
import sys
import csv
import numpy as np
import pandas as pd
import pathlib
from dataclasses import dataclass

# special symbols
DIR_SEP = '/'

# paths
# assuming the structure PROJ_DIR/dataproc/utils.py
PROJ_DIR = str(pathlib.Path(__file__).parent.parent.resolve())
#PROJ_DIR = os.path.abspath(os.pardir) 
DATAPROC_DIR = 'dataproc'
BIODATA_DIR = 'biodata'
# default directory with results'
RES_DIR  = 'results'
FIG_DIR = 'figures'
JAVA_CLASS = 'simulator.SimulatorCLI'
JAVA_CLASS_PATH = os.path.join('out', 'production', 'GRiPE')
SITE_OCCUPANCY_SCRIPT = 'site_occupancy.py'

# parts of file names
OCCUPANCY = '_occupancy_'
PARAMS = '_params_'
ENERGY = '_affinity_landscape_'
TF_SPECIES = '_TF_species_'

# file fields
NONSPEC_WAITING_TIME = 'SPECIFICWAITINGTIME'
TOTAL_SIM_TIME = 'STOP_TIME'
ENSAMBLE_SIZE = 'ENSAMBLE_SIZE'
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
COPY_NUMBER = 'COPYNUMBER'

# output file names
PLOT_PROB_SITE = 'sites_probabilities_'
PLOT_PROB_ENERGY = 'energy_probabilities_'


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


@dataclass
class TargetSite:
    repressor: bool
    name: str
    pos: int
    size: int
    dir: int
    pwm_energy: float
        
def bp_to_idx(bp):
    bps = list('acgt')
    return bps.index(bp)

def complementary(bp):
    compl_bps = list('tgca')
    return compl_bps[bp_to_idx(bp)]

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

def double_each_element(array):
    n = len(array)
    return np.array([[el,]*2 for el in array]).reshape(2*n)

# reading functions
def create_target_site_df(ts_filename, df_tf_data, df_site_energy):
    df = pd.DataFrame(columns=['repressor', 'name', 'name_strand', 'pos', 'size', 'strand', 'energy'])
    f = open(ts_filename)
    i = 0
    for r in f:
        ts_info = r.split(':')
        ts_name = ts_info[0]
        is_site_of_repr = float(df_tf_data.at[ts_name, REPRESSION_RATE]) > 0.0
        ts_pos = int(ts_info[2].split('..')[0])
        ts_size = int(ts_info[2].split('..')[1]) + 1 - ts_pos
        ts_dir = int(ts_info[-1])
        ts_name_dir = ts_name + ("5'3'" if ts_dir == 0 else "3'5'")
        ts_energy = df_site_energy.at[ts_pos, ts_name_dir]
        df.loc[i] = [is_site_of_repr, ts_name, ts_name_dir, ts_pos, ts_size, ts_dir, ts_energy]
        i += 1
        #ts_list.append(TargetSite(is_site_of_repr, ts_info[0], ts_pos, ts_size, ts_dir, ts_energy))
    f.close
    return df


def get_name(experiment_dir, file_prefix, name, experiment_id, ext=''):
    if ext and not ext.startswith('.'):
        ext = '.' + ext
    return os.path.join(experiment_dir, file_prefix + name + str(experiment_id) + ext)

def read_wig_to_dataframe(filename, index=None):
    f = open(filename, 'r')
    f.readline() # skip first line
    str_header = f.readline() # read header
    list_header = str_header.replace('"', '').strip().split(', ')
    result = pd.read_csv(f, names=list_header, index_col=index)
    f.close()
    return result

def read_wig_to_numpy(filename):
    df = read_wig_to_dataframe(filename)
    return df.to_numpy()

def read_csv_to_dataframe(filename, sep=','):
    f = open(filename, 'r')
    reader = csv.reader(f, delimiter=sep, quotechar='"', skipinitialspace=True)
    list_header = next(reader)[1:]
    list_data = []
    list_index = []
    for row in reader:
        list_index.append(row[0])
        list_data.append(row[1:])
    f.close()
    df = pd.DataFrame(list_data, columns=list_header, index=list_index)
    return df

def read_params_to_dict(filename):
    f = open(filename, 'r')
    d = {}
    for row in f:
        row = row.strip()
        if row and not row.startswith('#'):
            l = row.replace(';', '').split()
            d.update({l[0]: l[2].replace('"', '')})
    f.close()
    return d

def get_motif_size(str_pwm):
    if not str_pwm.startswith('PWM:'):
        print('Non-PWM string is not supported.')
        sys.exit(1)
    # pwm
    return len(str_pwm.split(';')[0].split(','))