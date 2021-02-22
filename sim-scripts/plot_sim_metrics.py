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
parser.add_argument('--dnmLoc', action='store', dest='dnmLoc', type=str, default='sim_loc_dnms.txt')
parser.add_argument('--pvalsPred', action='store', dest='pvalsPred', type=str, default='sim2_chrT_pval_sim.csv')
parser.add_argument('--pvalsTrue', action='store', dest='pvalsTrue', type=str, default='sim2_chrT_pval_sim_truth.csv')
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
            return (the_storage_arr, False, False)
        else:
            return (the_storage_arr, False, True)
    else:

        return (the_storage_arr, True, False)

def get_from_pval_file(fullDir, outputFile_csv):
    df = pd.read_csv("{}{}".format(fullDir, outputFile_csv))
    h1_s, h2_s = df["h1_count"], df["h2_count"]
    x_s, y_s = df["genomic_position"], -np.log10(df["pval"])
    color_vals = np.maximum(h1_s, h2_s)/(h1_s+h2_s)*100
    return x_s, y_s, color_vals

def get_from_dnm_file(fullDir, outputFile_dnm):
    dnm_loc_info = open("{}{}".format(fullDir, outputFile_dnm)).readlines()
    dnm_locs = np.array([int(dnm_loc_spec.split(":")[-1]) for dnm_loc_spec in dnm_loc_info[0].split(';')])
    to_keep_ind = np.array([True if dnm_bool_line.split(":")[-1].strip() == "FALSE" else False for dnm_bool_line in dnm_loc_info[1:]])
    dnm_locs = dnm_locs[to_keep_ind]
    return dnm_locs

def find_missing_genotype_rate(coverages):
    mgr = poisson.pmf(0, coverages)*100
    return np.around(mgr, 2)

mgrs = find_missing_genotype_rate(args.coverages)
cov_ind = np.argsort(mgrs)
mgrs = mgrs[cov_ind]
coverages = np.array(args.coverages)[cov_ind]
ngam_ind = np.argsort(args.ngams)[::-1]
ngams = np.array(args.ngams)[ngam_ind]
nsnp_ind = np.argsort(args.nsnps)
nsnps = np.array(args.nsnps)[nsnp_ind]

need_to_avg = False
if len(args.randsd) > 1:
    need_to_avg = True
    data_storage_arr = np.full((len(ngams), len(coverages), len(nsnps), len(args.randsd)), np.nan)
    for l, randsd in enumerate(args.randsd):
        for i, ngam in enumerate(ngams):
            for k, nsnp in enumerate(nsnps):
                for j, coverage in enumerate(coverages):
                    fullDir = "{}{}_{}_{}_{}_{}_{}_{}_{}_{}_{}_{}_{}_{}_{}_{}/".format(args.baseDir, args.sampleName, args.chrom, randsd, args.seqError, args.ase, args.sea, args.windowLength, ngam, nsnp, coverage, args.rlam, args.adnm, args.dnl, args.dna, args.dnb)
                    data_storage_arr, bool_val, bool_val2 = add_to_data_storage(data_storage_arr, need_to_avg, fullDir, args.metricOI, i, j, k, l)
else:
    figp, axesp = plt.subplots(nrows=len(ngams), ncols=len(coverages), sharey=True)
    cm = plt.cm.get_cmap("plasma")
    random_snp = np.random.choice(np.arange(len(nsnps)), 1)
    num_no_ospath = 0
    num_fail = 0
    data_storage_arr = np.full((len(ngams), len(coverages), len(nsnps)), np.nan)
    for i, ngam in enumerate(ngams):
        for k, nsnp in enumerate(nsnps):
            for j, coverage in enumerate(coverages):
                fullDir = "{}{}_{}_{}_{}_{}_{}_{}_{}_{}_{}_{}_{}_{}_{}_{}/".format(args.baseDir, args.sampleName, args.chrom, args.randsd[0], args.seqError, args.ase, args.sea, args.windowLength, ngam, nsnp, coverage, args.rlam, args.adnm, args.dnl, args.dna, args.dnb)
                data_storage_arr, bool_val, bool_val2 = add_to_data_storage(data_storage_arr, need_to_avg, fullDir, args.metricOI, i, j, k, 0)
                if bool_val:
                    num_no_ospath += 1
                if bool_val2:
                    num_fail += 1
                if k == random_snp:
                    x_s_p, y_s_p, colors_p = get_from_pval_file(fullDir, args.pvalsPred)
                    x_s_t, y_s_t, colors_t = get_from_pval_file(fullDir, args.pvalsTrue)
                    posp = axesp[i,j].scatter(x_s_p, y_s_p, c=colors_p, cmap=cm)
                    clbp = plt.colorbar(posp, ax=axesp[i,j])
                    clbp.set_label('td ratio', labelpad=-40, y=1.05, rotation=0)
                    if args.adnm == "TRUE":
                        dnm_locs = get_from_dnm_file(fullDir, args.dnmLoc)
                        for i, ind in enumerate(dnm_locs):
                            axesp[i,j].axvline(np.where(x_s_p == ind)[0], linestyle='--', color='black')
                    axeps[i,j].set_xlabel('SNP index')
                    axesp[i,j].set_ylabel('-log10(pval)')
    plt.tight_layout()
    fig.savefig('panel_manhattan_plot.png')
    plt.close(fig)


if need_to_avg:
    data_storage_arr = np.mean(data_storage_arr, axis=3)

print(num_no_ospath, num_fail, flush=True)
f=open('sim2_storage_arr.npz', 'wb')
np.savez(f, sim2_results=data_storage_arr)
f.close()

fig, axes = plt.subplots(nrows=1, ncols=len(nsnps), figsize=(20,10), sharey=True)
for k in np.arange(len(nsnps)):
    g = sns.heatmap(data_storage_arr[:,:,k], ax=axes[k], cmap="plasma", linecolor='black', linewidths=0.1, yticklabels=ngams, xticklabels=mgrs)
    axes[k].set_title("{} SNPs".format(nsnps[k]))
    if k == 0:
        axes[k].set_ylabel("Number of Gametes")
        g.set_yticklabels(g.get_yticklabels(), rotation = 90)
    axes[k].set_xlabel("Missing Genotype Rate")
    axes[k].tick_params(axis='x', labelrotation=45)
    g.set_facecolor('gray')

plt.tight_layout()
fig.savefig('simulation_varying_coverage_sperm_{}.png'.format(args.metricOI))
plt.close(fig)
