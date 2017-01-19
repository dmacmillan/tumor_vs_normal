__version__ = '1.0.0'

import argparse
import os
import sys
import logging
import time
import random
import colorsys
import subprocess
import cPickle as pickle
import pysam
import firebrowse
from Matrix import *

# sample types from firehose
# TP == only primary tumor samples
# TM == only metastatic tumor samples
# TR == only tumor recurrance samples
# NT == only tissue normal samples
# NB == only blood normal samples

def createExpressionMapping(expr):
    result = {}
    m = Matrix.readMatrix(expr)
    for i in xrange(m.m):
        gene = m.rownames[i]
        if gene not in result:
            result[gene] = {}
        for j in xrange(m.n):
            sample = m.colnames[j]
            result[gene][sample] = m.mtx[i][j]
    return result

def buildExpression(genes, mycohort):
    result = {}
    for mygene in genes:
        try:
            values = firebrowse.Samples().mRNASeq(gene=mygene, cohort=mycohort, format=firebrowse.CODEC_DJSON)
        except ValueError:
            continue
        for sample in values['mRNASeq']:
            stype = sample['sample_type']
            name = sample['tcga_participant_barcode']
            expr = sample['expression_log2']
            if not expr:
                continue
            if mygene not in result:
                result[mygene] = {name: {stype: [expr]}}
            elif name not in result[mygene]:
                result[mygene][name] = {stype: [expr]}
            elif stype not in result[mygene][name]:
                result[mygene][name][stype] = [expr]
            else:
                result[mygene][name][stype].append(expr)
    return result

def getFirebrowseExpression(mycohort, mygene, samples):
    values = firebrowse.Samples().mRNASeq(gene=mygene, cohort=mycohort, format=firebrowse.CODEC_DJSON)
    values = [x['expression_log2'] for x in values['mRNASeq'] if x['tcga_participant_barcode'] in samples]
    return values

def genRHeatmap(gene, path, clust='TRUE'):
    out = os.path.join(path, 'plot_{}.r'.format(gene))
    with open(out, 'w') as f:
        f.write(
            'setwd("{}")\n' \
            'library(pheatmap)\n' \
            'library(RColorBrewer)\n' \
            'mtx = read.table("{}")\n' \
            'mybreaks = seq(min(mtx), max(mtx), by = 1)\n' \
            'colors = rev(brewer.pal(max(min(length(mybreaks),11),3), "RdYlBu"))\n' \
            'png("{}.png", width=1000, height=1000, res=150)\n' \
            'pheatmap(mtx, main="{}", legend_breaks = mybreaks, cluster_rows = {}, cluster_cols = {})\n' \
            'dev.off()'.format(path, gene, gene, gene, clust, clust)
        )
    return out

def sprint(text):
    sys.stdout.write(text)
    sys.stdout.flush()

def centroid(_list):
    return float(sum(_list))/len(_list)

# List of lists must be sorted in ascending order
def getListIndexMinMetric(_list, linkage='centroid', metric='euclidian'):
    if metric == 'euclidian':
        index = None
        _min = float('inf')
        for i in xrange(len(_list) - 1):
            if linkage == 'centroid':
                dist = abs(centroid([x.end for x in _list[i+1]]) - centroid([x.end for x in _list[i]]))
            elif linkage == 'single':
                dist = float('inf')
                for p in xrange(len(_list[i+1])):
                    for q in xrange(len(_list[i])):
                        temp = abs(_list[i+1][p].end - _list[i][q].end)
                        if temp < dist:
                            dist = temp
            if dist < _min:
                index = i
                _min = dist
        return index, _min

# Given a sorted list, a linkage criteria, 
# and a threshold, return a
# list of lists where each internal list
# represents a cluster
def iterAHC(_list, linkage='centroid', threshold=20, metric='euclidian'):
    clusters = sorted([x for x in _list], key = lambda x: x.end)
    clusters = [[x] for x in clusters]
    index, _min = getListIndexMinMetric(clusters)
    while (_min <= threshold) and (len(clusters) > 1):
        clusters[index] += clusters[index+1]
        del(clusters[index+1])
        index, _min = getListIndexMinMetric(clusters, linkage=linkage, metric=metric)
    return clusters

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Given KLEAT files identify tissue-specific cleavage sites as well as differentially utilized sites')
    parser.add_argument('kleat_bed_tumor', help='A tabix-indexed BED file containing all cancerous KLEAT calls you wish to analyze for a given tissue')
    parser.add_argument('kleat_bed_normal', help='A tabix-indexed BED file containing all normal KLEAT calls you wish to analyze from all tissues')
    parser.add_argument("cohort", choices=['BRCA', 'LUSC', 'LUAD'], help='Set the cohort, this is used for expression fetching')
    parser.add_argument('-me', '--min_expression', type=int, default=1, help='Reject cleavage site call for gene if the expression level origin sample for that gene is less than this threshold. Default == 1')
    parser.add_argument('-ms', '--max_sites', type=int, default=2, help='Maximum number of cleavage sites (clusters) per 3\'UTR in either tumor or normal. Default == 2')
    parser.add_argument('-mc', '--min_sites_cluster', type=int, default=2, help='The minimum number of sites per cluster. Default == 2')
    parser.add_argument('-u', '--utr3s', help='Bgzipped tabix-indexed BED file containing all viable 3\'UTRs to use')
    parser.add_argument('-n', '--normalize', action='store_true', help='Normalize frequencies by number of samples in tissue')
    parser.add_argument("-l", "--log", dest="logLevel", default='WARNING', choices=['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL'], help='Set the logging level. Default = "WARNING"')
    parser.add_argument('-o', '--outdir', default=os.getcwd(), help='Path to output to. Default is {}'.format(os.getcwd()))
    
    args = parser.parse_args()

    #utr3s = pysam.tabix_iterator(open(args.utr3s), parser=pysam.asGTF())
    utr3s = pysam.tabix_iterator(open(args.utr3s), parser=pysam.asBed())
    genes = set()
    for utr3 in utr3s:
        genes.add(utr3.name)
    utr3s = pysam.tabix_iterator(open(args.utr3s), parser=pysam.asBed())

    if args.normalize:
        expr_path = os.path.join(args.outdir, '.expression')
        if not os.path.isfile(expr_path):
            expression = buildExpression(genes, args.cohort)
            pickle.dump(expression, open(expr_path, 'wb'))
        else:
            sprint('Loading expression data ... ')
            expression = pickle.load(open(expr_path, 'rb'))
            print 'DONE'

    kleats_tumor = pysam.TabixFile(args.kleat_bed_tumor, parser=pysam.asBed())
    kleats_normal = pysam.TabixFile(args.kleat_bed_normal, parser=pysam.asBed())

    count_normal_samples = pysam.tabix_iterator(open(args.kleat_bed_normal), parser=pysam.asBed())
    normal_samples = set()
    for sample in count_normal_samples:
        normal_samples.add(sample.name)

    nsamples = len(normal_samples)

    if not os.path.isdir(args.outdir):
        try:
            os.makedirs(args.outdir)
        except OSError:
            pass

    # Logging
    logging.basicConfig(filename=os.path.join(args.outdir, 'log_{}'.format(os.path.basename(__file__))), level=getattr(logging, args.logLevel), filemode='w')

    mapping = {}
    for utr3 in utr3s:
        skip = False
        logging.debug('utr3: {}'.format(utr3))
        gene = utr3.name
        #gene = utr3.asDict()['gene_name']
        logging.debug('gene: {}'.format(gene))
        lutr3 = utr3.end - utr3.start
        try:
            tumor_calls = [x for x in kleats_tumor.fetch(utr3.contig, utr3.start-20, utr3.end+20)]
            logging.debug('tumor_calls: {}'.format([x.end for x in tumor_calls]))
        except ValueError:
            tumor_calls = None
        try:
            normal_calls = [x for x in kleats_normal.fetch(utr3.contig, utr3.start-20, utr3.end+20)]
            logging.debug('normal_calls: {}'.format([x.end for x in normal_calls]))
        except ValueError:
            normal_calls = None
        if not normal_calls and not tumor_calls:
            logging.debug('No calls for normal or tumor')
            continue
        for call in normal_calls:
            call.name = '{}_normal'.format(call.name)
        for call in tumor_calls:
            call.name = '{}_tumor'.format(call.name)
        tkeys = ('tumor', 'normal')
        clusters = iterAHC([x for x in tumor_calls + normal_calls])
#        with open('clusters.single.100.bedgraph', 'w') as o:
#            for c in clusters:
#                value = int(centroid([x.end for x in c]))
#                o.write('{}\t{}\t{}\t{}\n'.format(utr3.contig, value-1, value, len(c)))
#        continue
        if any([len(x) < args.min_sites_cluster for x in clusters]):
            logging.debug('Not enough sites in cluster, skipping')
            continue
        nclusters = len(clusters)
        logging.debug('nclusters = {}'.format(nclusters))
        if nclusters > args.max_sites:
            logging.debug('Number of clusters exceeds maximum allowed, skipping')
            continue
        mtx = Matrix(2, nclusters)
        centroids = []
        for cluster in clusters:
            values = [x.end for x in cluster]
            centroids.append(centroid(values))
        logging.debug('centroids: {}'.format(centroids))
        for i, _type in enumerate(tkeys):
            _type = tkeys[i]
            logging.debug('type: {}'.format(_type))
            for j in xrange(mtx.n):
                values = [x.end for x in clusters[j] if x.name.split('_')[1] == _type]
                logging.debug('values: {}'.format(values))
                if not values:
                    continue
                if args.normalize:
                    if _type == 'tumor':
                        ttype = 'TP'
                    elif _type == 'normal':
                        ttype = 'NT'
                    try:
                        logging.debug('expression: {}'.format([expression[gene][x.name.split('_')[0]] for x in clusters[j]]))
                        expr = [expression[gene][x.name.split('_')[0]][ttype] for x in clusters[j] if x.name.split('_')[1] == _type]
                        expr = [y for x in expr for y in x]
                    except KeyError:
                        expr = None
                    if not expr:
                        skip = True
                        break
                    expression_mean = sum(expr)/len(expr)
                    logging.debug('expression mean: {}'.format(expression_mean))
                cs = centroid(values)
                score = len(values)
                #score /= nsamples
                logging.debug('nsamples: {}'.format(nsamples))
                logging.debug('pre-normal: {}'.format(score))
                if args.normalize:
                    score /= float(expression_mean)
                mtx.mtx[i][j] = score
                logging.debug('post-normal: {}'.format(score))
        if skip:
            continue
        mtx.rownames = tkeys
        if utr3.strand == '+':
            mtx.colnames = ['c_{}_{}'.format(int(x), int(100*abs(x - utr3.start)/(lutr3))) for x in centroids]
        else:
            mtx.colnames = ['c_{}_{}'.format(int(x), int(100*abs(utr3.end - x)/(lutr3))) for x in centroids]
        if (mtx.n == 1) or (mtx.m < 2):
            continue
        gname = '{}_{}_{}_{}'.format(gene, utr3.start, utr3.end, utr3.strand)
        with open(os.path.join(args.outdir, gname), 'w') as f:
            f.write(str(mtx))
