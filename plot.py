#!/usr/bin/env python

import pylab as pl
import os

import config

def plot_in_one():
    f1 = pl.figure()
    f2 = pl.figure()
    ax1 = f1.add_subplot(111)
    ax1.plot(range(0, 10))

    ax2 = f2.add_subplot(111)
    ax2.plot(range(10, 20))

    pl.show()


def iteration_recallprecision_curve():
    name = "graphlet96"
    folder = config.folder + "results/%s/filtered_iterated_wang/" % name
    max_go_num = 9
    iter_num = 10
    pl.clf()

    recall_precision_list = []
    for i in range(1, max_go_num+1):
        recall_precision = []
        recall = []
        precision = []
        fpath = os.path.join(folder, "%d.csv"%i)
        ifile = open(fpath, "r")
        for line in ifile:
            line = line.strip()
            if len(line)>0:
                arr = line.split(",")
                recall.append(float(arr[0]))
                precision.append(float(arr[1]))
        ifile.close()
        for iter in range(0, iter_num):
            recall_avg = 0.0
            precision_avg = 0.0
            for cv in range(0, 10):
                recall_avg += recall[cv*10+iter]
                precision_avg += precision[cv*10+iter]
            recall_avg /= 10
            precision_avg /= 10
            recall_precision.append((recall_avg, precision_avg))
        recall_precision_list.append(recall_precision)

    for iter in range(0, iter_num):
        recall = []
        precision = []
        for i in range(0, max_go_num): # GONUM
            recall.append(recall_precision_list[i][iter][0])
            precision.append(recall_precision_list[i][iter][1])

        pl.plot(recall, precision, '--+', label='Iteration %d'%iter, lw=1)

    pl.xlabel('Recall')
    pl.ylabel('Precision')
    pl.xlim([0.3, 0.48])
    pl.ylim([0.45, 0.55])
    pl.title('Precision-Recall')
    pl.legend(loc="lower right")
    pl.savefig("%s.pdf" % name)
    pl.show()


def precisionrecall_curve(recall, precision):
    pl.clf()
    pl.plot(recall, precision, '--+', label='Precision-Recall Curve')

    pl.xlabel('Recall')
    pl.ylabel('Precision')
    pl.xlim([0.2, 1.0])
    pl.ylim([0.2, 1.0])
    pl.title('Precision-Recall')
    pl.legend(loc="lower left")
    pl.show()

def precisionrecall_curve_multiple():
    pl.clf()

    pl.plot(recall1, precision1, '--o', label='iterated_weighted_%s'%config.sim_metric, lw=1)
    #pl.plot(recall2, precision2, '--+', label='weighted_%s'%config.sim_metric, lw=1)
    pl.plot(recall3, precision3, '--v', label='iterated_sim_%s'%config.sim_metric, lw=1)
    pl.plot(recall4, precision4, '--^', label='NC', lw=1)
    
    pl.xlabel('Recall')
    pl.ylabel('Precision')
    pl.xlim([0.3, 0.9])
    pl.ylim([0.3, 0.9])
    #pl.title('Precision-Recall plot')
    pl.legend(loc="lower left")

    pl.savefig("plots/%s.pdf"%config.sim_metric)

    #pl.show()


def get_precision_recall(folder):
    recall = []
    precision = []

    for i in range(1, 21):
        fpath = os.path.join(folder + "%d.csv"%i)
        ifile = open(fpath, "r")
        recall_avg = 0.0
        precision_avg = 0.0
        i = 0
        for line in ifile:
            line = line.strip()
            if len(line)>0:
                i+=1
                arr = line.split(",")
                recall_avg += float(arr[0])
                precision_avg += float(arr[1])
        ifile.close()
        recall_avg /= i
        precision_avg /= i
        recall.append(recall_avg)
        precision.append(precision_avg)
    return (recall, precision)


if __name__ == "__main__":
    result_folder = config.folder + "results/old/"
    folder1 = result_folder + config.filtered + "iterated_weighted_" + config.sim_metric + "/"
    folder2 = result_folder + config.filtered + "weighted_" + config.sim_metric + "/"
    folder3 = result_folder + config.filtered + "iterated_" + config.sim_metric + "/"
    folder4 = result_folder + config.filtered + "mv_wang/"

    
    (recall1, precision1) = get_precision_recall(folder1)
    (recall2, precision2) = get_precision_recall(folder2)
    (recall3, precision3) = get_precision_recall(folder3)
    (recall4, precision4) = get_precision_recall(folder4)

    #precisionrecall_curve_multiple()

    #plot_in_one()

    iteration_recallprecision_curve()

