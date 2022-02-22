#! ~/anaconda3/bin/python

from utils import *
import numpy as np
#import scipy as sp
import pandas as pd
import matplotlib.pyplot as plt
import argparse

def main(dir, prefix, id):
    # load occupancy (occupancy.wig)
    df_occupancy = read_wig_to_dataframe(get_name(dir, prefix, OCCUPANCY, id, 'wig'))
    # load energies (PFM-based)
    df_site_energy = read_wig_to_dataframe(get_name(dir, prefix, ENERGY, id, 'wig'))
    # load tf data
    df_tf_data = read_csv_to_dataframe(get_name(dir, prefix, TF_SPECIES, id, 'csv')) 
    # load params
    dict_params = read_params_to_dict(get_name(dir, prefix, PARAMS, id, 'grp'))
    
    # probability of each site
    tf_list = np.array(df_tf_data.index, dtype=str)
    tf_num = len(tf_list)
    tf_size = np.array([get_motif_size(df_tf_data.at[tf_list[i], MOTIF])
                        + int(df_tf_data.at[tf_list[i], SIZE_LEFT])
                        + int(df_tf_data.at[tf_list[i], SIZE_RIGHT]) for i in range(tf_num)], dtype=int)
    nonspecific_time = np.array(df_tf_data[NONSPEC_WAITING_TIME], dtype=float)
    total_time = float(dict_params[TOTAL_SIM_TIME])
    assoc_rate = np.array(df_tf_data[ASSOC_RATE], dtype=float)
    unbind_prob = np.array(df_tf_data[UNBIND_PROB], dtype=float) * np.array(df_tf_data[JUMP_PROB], dtype=float)
    energy = df_site_energy.to_numpy().transpose()[1:]
    for i in range(tf_num):
        energy[:,-tf_size[i]+1:] = +np.inf
    energy_exponent = np.exp(-energy)
    prob_from_energy = energy_exponent / (unbind_prob/assoc_rate/nonspecific_time + np.sum(energy_exponent, axis=1))
    prob_from_exper = df_occupancy.to_numpy().transpose()[2:] / total_time

    for tf_i in range(tf_num):
        sort_idx = np.argsort(energy[tf_i])
        #energy_cont = np.linspace(energy[tf_i][sort_idx][0], energy[tf_i][sort_idx][-tf_size[tf_i]], 101)
        #prob_from_energy_cont = np.exp(-energy_cont) / (unbind_prob[tf_i]/assoc_rate[tf_i]/nonspecific_time[tf_i] + np.exp(-energy_cont[0]) - np.exp(-energy_cont[-1]))
        
        plt.figure(figsize=(4,3))
        #plt.plot(energy_cont, prob_from_energy_cont, label='Boltzmann: continuous')
        plt.plot(energy[tf_i][sort_idx], prob_from_energy[tf_i][sort_idx], marker='.', label='Boltzmann: discrete')
        plt.plot(energy[tf_i][sort_idx], prob_from_exper[tf_i][sort_idx], marker='.', label='experimental')
        plt.xlabel('energy')
        plt.ylabel('probability')
        plt.legend()
        plt.title(tf_list[tf_i])
        plt.tight_layout()
        plt.savefig(PROJ_DIR + DIR_SEP + FIG_DIR + tf_list[tf_i] + '_' + PLOT_PROB_ENERGY + str(id) + '.pdf')

        plt.figure(figsize=(4,3))
        plt.plot(prob_from_energy[tf_i], marker='.', label='Boltzmann')
        plt.plot(prob_from_exper[tf_i], marker='.', label='experiment')
        plt.xlabel('site position')
        plt.ylabel('probability')
        plt.legend()
        plt.title(tf_list[tf_i])
        plt.tight_layout()
        plt.savefig(PROJ_DIR + DIR_SEP + FIG_DIR + tf_list[tf_i] + '_' + PLOT_PROB_SITE + str(id) + '.pdf')

    # plot probability for each site
    # plt.figure(figsize=(4,3))
    # plt.plot(prob_from_energy, marker='o', label='Boltzmann')
    # plt.plot(prob_from_exper, marker='o', label='experiment')
    # plt.xlabel('site position')
    # plt.ylabel('probability')
    # plt.legend()
    # plt.tight_layout()
    # plt.savefig(PROJ_DIR + DIR_SEP + FIG_DIR + PLOT_PROB_SITE + str(id) + '.pdf')

    # # plot probability for each PFM-energy
    # energy = df_site_energy['Kr'].to_numpy()[:-2]
    # sort_idx = np.argsort(energy)
    # energy_cont = np.linspace(energy[sort_idx][0], energy[sort_idx][-1], 101)
    # prob_from_energy_cont = np.exp(-energy_cont) / (unbind_prob/assoc_rate/nonspecific_time + np.exp(-energy_cont[0]) - np.exp(-energy_cont[-1]))
    # plt.figure(figsize=(4,3))
    # plt.plot(energy_cont, prob_from_energy_cont, label='Boltzmann: continuous')
    # plt.plot(energy[sort_idx], prob_from_energy[sort_idx], marker='o', label='Boltzmann: discrete')
    # plt.plot(energy[sort_idx], prob_from_exper[sort_idx], marker='o', label='experimental')
    # plt.xlabel('energy')
    # plt.ylabel('probability')
    # plt.legend()
    # plt.tight_layout()
    # plt.savefig(PROJ_DIR + DIR_SEP + FIG_DIR + PLOT_PROB_ENERGY + str(id) + '.pdf')


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--id", type=int, help="unique long number of simulation", required=True)
    parser.add_argument("-p", "--prefix", type=str, help="files prefix", required=True)
    parser.add_argument("-d", "--dir", type=str, help="files directory (relative path)", required=True)
    args = parser.parse_args()
    dir = PROJ_DIR + DIR_SEP + args.dir
    prefix = args.prefix
    id = args.id
    print('files directory:', dir)
    print('files prefix:   ', prefix)
    print('files id:       ', id)
    main(dir, prefix, id)