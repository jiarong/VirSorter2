#!/usr/bin/env python

import sys
import os
import pandas as pd
import numpy as np
import warnings
import joblib

script_dir = os.path.dirname(os.path.abspath(__file__))
snakefile_dir = os.path.dirname(script_dir)
pkg_dir = os.path.dirname(snakefile_dir)
sys.path.append(pkg_dir)
from virsorter.config import get_default_config

DEFAULT_CONFIG = get_default_config()

TOTAL_FEATURE_LIST = DEFAULT_CONFIG['TOTAL_FEATURE_LIST']
SELECT_FEATURE_LIST = DEFAULT_CONFIG['SELECT_FEATURE_LIST']
CLASSIFY_THREADS = DEFAULT_CONFIG['CLASSIFY_THREADS'] 

np.random.seed(seed=99)

cpu = CLASSIFY_THREADS  ### add this to CL interface

def main():
    ''' classify seqs as viral or nonviral based on features. 

    Example:
        python classify.py <all.pdg.ftr> <model> label <out.tsv>

        <all.pdg.ftr>: feature file
        <model>: RandomForest model trained for a viral group
        label: viral group name
        <out.tsv>: tab delimited file with first col as seqname, second
            col as probability of being viral

    '''
    if len(sys.argv) != 5:
        mes = '*** Usage: python {} <all.pdg.ftr> <model> label <out.tsv>\n'
        sys.stderr.write(mes.format(os.path.basename(sys.argv[0])))
        sys.exit(1)

    ftr_f = sys.argv[1]
    model_f = sys.argv[2]
    label = sys.argv[3]
    outfile = sys.argv[4]

    X = pd.read_csv(ftr_f, sep='\t', header=0, index_col=0)
    if len(X) == 0:
        cols = ['seqname', label]
        with open(outfile, 'w') as fw:
            fw.write('{}\n'.format('\t'.join(cols)))
            sys.exit(0)

    X = X.loc[:, SELECT_FEATURE_LIST]
    with warnings.catch_warnings():
        warnings.filterwarnings('ignore', category=DeprecationWarning)
        warnings.filterwarnings('ignore', category=FutureWarning)
        model = joblib.load(model_f)
        try:
            # pipe on grid search
            model.named_steps.gs.best_estimator_.set_params(n_jobs=cpu)
        except AttributeError as e:
            # grid search on pipe
            # steps = [('scaler', MinMaxScaler()), ('rf', clf_to_train)]
            model.named_steps.rf.set_params(n_jobs=cpu)
        pred_prob = model.predict_proba(X)
        labels = model.classes_
        df = pd.DataFrame(pred_prob, columns= labels, index=X.index)
        df = df.rename({0: 'nonviral', 1: label}, axis=1, errors='raise')
        df = df.drop('nonviral', axis=1)
        # print result 2 column table: seqname and predicted prob as group
        df.to_csv(outfile, sep='\t', index_label='seqname')

if __name__ == '__main__':
    main()
