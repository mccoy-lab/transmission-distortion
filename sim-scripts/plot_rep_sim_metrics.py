#!/usr/bin/env python3

import numpy as np
import argparse as ap
import os.path

parser = ap.ArgumentParser(description="plot the simulation results")
parser.add_argument('--ngam', action='store', dest='ngam', type=int, default=1000)
parser.add_argument('--nsnp', action='store', dest='nsnp', type=int, default=30000)
parser.add_argument('--coverage', action='store', dest='coverage', type=float, default=0.01)
parser.add_argument('--windowLength', action='store', dest='windowLength', type=int, default=3000)
parser.add_argument('--randsds', action='store', nargs='+', dest='randsds', type=int, default=[27, 42, 386, 651, 1059, 2556, 2563, 2862, 3417, 4900])
parser.add_argument('--metricOIs', action='store', dest='metricOIs', nargs='+', type=str, default=['gam_hap_rec_raw_acc','gam_hap_rec_cor_acc', 'parent_hap_rec_acc', 'cons_f1', 'lib_f1', 'cons_accuracy', 'lib_accuracy', 'cons_specificity', 'lib_specificity', 'nrecomb', 'nsnps'], help='one of {gam_hap_rec_cor_acc, gam_hap_rec_raw_acc, parent_hap_rec_acc, lib_precision, cons_precision, lib_recall, cons_recall, lib_accuracy, cons_accuracy, lib_f1, cons_f1, lib_specificity, cons_specificity, lib_fdr, cons_fdr, lib_fpr, cons_fpr, nrecomb, nsnps}')
parser.add_argument('--dirBase', action='store', dest='dirBase', type=str, default='/home/kweave23/gamete_data/run_sim3_20210218/')
args = parser.parse_args()

def file_based_on_metricOI(metricOI):
    dict_for_file = {"gam_hap_rec_cor_acc":"sim3_hap_rec_{}.txt",
                     "gam_hap_rec_raw_acc":"sim3_hap_rec_{}.txt",
                     "parent_hap_rec_acc":"sim3_hap_rec_{}.txt",
                     "lib_precision":"sim3_lib_recomb_{}.txt",
                     "cons_precision":"sim3_cons_recomb_{}.txt",
                     "lib_recall":"sim3_lib_recomb_{}.txt",
                     "cons_recall":"sim3_cons_recomb_{}.txt",
                     "lib_accuracy":"sim3_lib_recomb_{}.txt",
                     "cons_accuracy":"sim3_cons_recomb_{}.txt",
                     "lib_f1":"sim3_lib_recomb_{}.txt",
                     "cons_f1":"sim3_cons_recomb_{}.txt",
                     "lib_specificity":"sim3_lib_recomb_{}.txt",
                     "cons_specificity":"sim3_cons_recomb_{}.txt",
                     "lib_fdr":"sim3_lib_recomb_{}.txt",
                     "cons_fdr":"sim3_cons_recomb_{}.txt",
                     "lib_fpr":"sim3_lib_recomb_{}.txt",
                     "cons_fpr":"sim3_cons_recomb_{}.txt",
                     "nrecomb": "sim3_nrecomb_{}.txt" ,
                     "nsnps": "sim3_nsnps_{}.txt"
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
                     "cons_fpr":7,
                     "nrecomb": 0,
                     "nsnps": 0
                    }
    return dict_for_index[metricOI]

def find_summary_stat(metricOI, values):
    if metricOI in ["gam_hap_rec_cor_acc", "gam_hap_rec_raw_acc", "parent_hap_rec_acc", "nrecomb", "nsnps", "lib_precision", "cons_precision", "lib_recall", "cons_recall", "lib_accuracy", "cons_accuracy", "lib_specificity", "cons_specificity", "lib_fdr", "cons_fdr", "lib_fpr", "cons_fpr"]:
        return "mean", np.nanmean(values)
    elif metricOI in ["lib_f1", "cons_f1"]:
        return "median", np.nanmedian(values)

for metricOI in args.metricOIs:
    values_to_store = []
    for randsd in args.randsds:
        file_base = 'wl{}_gam{}_snp{}_cov{}/'.format(args.windowLength, args.ngam, args.nsnp, args.coverage)
        file = args.dirBase + file_base + file_based_on_metricOI(metricOI).format(randsd)
        if os.path.exists(file) and len(open(file).readlines()) > 0:
            values_to_store.append(float(open(file).readlines()[file_line_index_based_on_metricOI(metricOI)].split(':')[-1]))
    stat_type, metric_val = find_summary_stat(metricOI, values_to_store)
    print(metricOI, stat_type, metric_val)
