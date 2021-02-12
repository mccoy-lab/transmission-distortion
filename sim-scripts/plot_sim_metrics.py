#!/usr/bin/env python3

import numpy as np
import os.path
from scipy.stats import poisson
import matplotlib.pyplot as plt
import seaborn as sns
import argparse as ap

parser = ap.ArgumentParser(description="plot the simulation results")
parser.add_argument('--ngams', action='store', nargs='+', dest='ngams', type=int, default=[3, 15, 50, 150, 500, 1000, 2500, 5000, 9600, 12800, 16000])
parser.add_argument('--nsnps', action='store', nargs='+', dest='nsnps', type=int, default=[5000, 10000, 30000, 50000, 100000])
parser.add_argument('--coverages', action='store', nargs='+', dest='coverages', type=float, default=[0.0001, 0.001, 0.01, 0.1, 0.357, 0.511, 0.693, 0.916, 1.204, 1.609, 2.303])
parser.add_argument('--sampleName', action='store', dest='sampleName', type=str, default='sim2')
parser.add_argument('--chrom', action='store', dest='chrom', type=str, default='chrT')
parser.add_argument('--seqError', action='store', dest='seqError', type=float, default=0.05)
parser.add_argument('--windowLength', action='store', dest='windowLength', type=int, default=3000)
parser.add_argument('--randsd', action='store', nargs='+', dest='randsd', type=int, default=[42])
parser.add_argument('--rlam', action='store', dest='rlam', type=int, default=1)
parser.add_argument('--ase', action='store', dest='ase', type=str, default='TRUE')
parser.add_argument('--adnm', action='store', dest='adnm', type=str, default='TRUE')
parser.add_argument('--sea', action='store', dest='sea', type=float, default=0.05)
parser.add_argument('--dnl', action='store', dest='dnl', type=int, default=5)
parser.add_argument('--dna', action='store', dest='dna', type=float, default=7.5)
parser.add_argument('--dnb', action='store', dest='dnb', type=int, default=10)
parser.add_argument('--baseDir', action='store', dest='baseDir', type=str, default='/home-3/kweave23@jhu.edu/work-rmccoy22/kweave23/sc_transmission_distortion/run_sim2_20210208/')
parser.add_argument('--metricOI', action='store', dest='metricOI', type=str, default='gam_hap_rec_cor_acc', help='one of {gam_hap_rec_cor_acc, gam_hap_rec_raw_acc, parent_hap_rec_acc, lib_precision, cons_precision, lib_recall, cons_recall, lib_accuracy, cons_accuracy, lib_f1, cons_f1, lib_specificity, cons_specificity, lib_fdr, cons_fdr, lib_fpr, cons_fpr}')
args = parser.parse_args()

def file_based_on_metricOI(metricOI):
    dict_for_file = {"gam_hap_rec_cor_acc":"sim_hap_reconstruction_acc.txt",
                     "gam_hap_rec_raw_acc":"sim_hap_reconstruction_acc.txt",
                     "parent_hap_rec_acc":"sim_hap_reconstruction_acc.txt",
                     "lib_precision":"sim_lib_recomb.txt",
                     "cons_precision":"sim_cons_recomb.txt",
                     "lib_recall":"sim_lib_recomb.txt",
                     "cons_recall":"sim_cons_recomb.txt",
                     "lib_accuracy":"sim_lib_recomb.txt",
                     "cons_accuracy":"sim_cons_recomb.txt",
                     "lib_f1":"sim_lib_recomb.txt",
                     "cons_f1":"sim_cons_recomb.txt",
                     "lib_specificity":"sim_lib_recomb.txt",
                     "cons_specificity":"sim_cons_recomb.txt",
                     "lib_fdr":"sim_lib_recomb.txt",
                     "cons_fdr":"sim_cons_recomb.txt",
                     "lib_fpr":"sim_lib_recomb.txt",
                     "cons_fpr":"sim_cons_recomb.txt"
                    }
    return dict_for_file[metricOI]

def file_line_index_based_on_metricOI(metricOI):
    dict_for_index = {"gam_hap_rec_cor_acc": -1,
                     "gam_hap_rec_raw_acc": -2,
                     "parent_hap_rec_acc":0,
                     "lib_precision":1,
                     "cons_precision":1,
                     "lib_recall":2,
                     "cons_recall":2,
                     "lib_accuracy":3,
                     "cons_accuracy":3,
                     "lib_f1":4,
                     "cons_f1":4,
                     "lib_specificity":5,
                     "cons_specificity":5,
                     "lib_fdr":6,
                     "cons_fdr":6,
                     "lib_fpr":7,
                     "cons_fpr":7
                    }
    return dict_for_index[metricOI]

def add_to_data_storage(the_storage_arr, need_to_avg, fullDir, metricOI, i, j, k, l):
    outputFile = file_based_on_metricOI(metricOI)
    if os.path.exists("{}{}".format(fullDir, outputFile)):
        metric_file_info = open("{}{}".format(fullDir, outputFile)).readlines()
        if len(metric_file_info) > 0:
            metric = float(metric_file_info[file_line_index_based_on_metricOI(metricOI)].split(":")[-1])
            if need_to_avg:
                the_storage_arr[i, j, k, l] = metric
            else:
                the_storage_arr[i, j, k] = metric
            return (the_storage_arr, False)
        else:
            return (the_storage_arr, False)
    else:

        return (the_storage_arr, True)

def find_missing_genotype_rate(coverages):
    mgr = poisson.pmf(0, coverages)*100
    return np.around(mgr, 2)

need_to_avg = False
if len(args.randsd) > 1:
    need_to_avg = True
    data_storage_arr = np.full((len(args.ngams), len(args.coverages), len(args.nsnps), len(args.randsd)), np.nan)
    for l, randsd in enumerate(args.randsd):
        for i, ngam in enumerate(args.ngams):
            for k, nsnp in enumerate(args.nsnps):
                for j, coverage in enumerate(args.coverages):
                    fullDir = "{}{}_{}_{}_{}_{}_{}_{}_{}_{}_{}_{}_{}_{}_{}_{}/".format(args.baseDir, args.sampleName, args.chrom, randsd, args.seqError, args.ase, args.sea, args.windowLength, ngam, nsnp, coverage, args.rlam, args.adnm, args.dnl, args.dna, args.dnb)
                    data_storage_arr, bool_val = add_to_data_storage(data_storage_arr, need_to_avg, fullDir, args.metricOI, i, j, k, l)
else:
    num_no_ospath = 0
    data_storage_arr = np.full((len(args.ngams), len(args.coverages), len(args.nsnps)), np.nan)
    for i, ngam in enumerate(args.ngams):
        for k, nsnp in enumerate(args.nsnps):
            for j, coverage in enumerate(args.coverages):
                fullDir = "{}{}_{}_{}_{}_{}_{}_{}_{}_{}_{}_{}_{}_{}_{}_{}/".format(args.baseDir, args.sampleName, args.chrom, args.randsd[0], args.seqError, args.ase, args.sea, args.windowLength, ngam, nsnp, coverage, args.rlam, args.adnm, args.dnl, args.dna, args.dnb)
                data_storage_arr, bool_val = add_to_data_storage(data_storage_arr, need_to_avg, fullDir, args.metricOI, i, j, k, 0)
                if bool_val:
                    num_no_ospath += 1

if need_to_avg:
    data_storage_arr = np.mean(data_storage_arr, axis=3)

print(num_no_ospath, flush=True)
f=open('sim2_storage_arr.npz', 'wb')
np.savez(f, sim2_results=data_storage_arr)
f.close()

fig, axes = plt.subplots(nrows=1, ncols=len(args.nsnps), figsize=(20,10), sharey=True)
for k in np.arange(len(args.nsnps)):
    g = sns.heatmap(data_storage_arr[:,:,k], ax=axes[k], cmap="plasma", linecolor='black', linewidths=0.1, yticklabels=args.ngams, xticklabels=find_missing_genotype_rate(args.coverages))
    axes[k].set_title("{} SNPs".format(args.nsnps[k]))
    if k == 0:
        axes[k].set_ylabel("Number of Gametes")
        g.set_yticklabels(g.get_yticklabels(), rotation = 90)
    axes[k].set_xlabel("Missing Genotype Rate")
    axes[k].tick_params(axis='x', labelrotation=45)
    g.set_facecolor('gray')

plt.tight_layout()
fig.savefig('simulation_varying_coverage_sperm_{}.png'.format(args.metricOI))
plt.close(fig)
