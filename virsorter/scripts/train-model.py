import sys
import os
import multiprocessing
import click
import pandas as pd
import joblib

script_dir = os.path.dirname(os.path.abspath(__file__))
snakefile_dir = os.path.dirname(script_dir)
pkg_dir = os.path.dirname(snakefile_dir)
sys.path.append(pkg_dir)
from virsorter.config import set_logger, get_default_config

DEFAULT_CONFIG = get_default_config()

from sklearn.preprocessing import MinMaxScaler
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import GridSearchCV
from sklearn.pipeline import Pipeline

set_logger()

TOTAL_FEATURE_LIST = list(DEFAULT_CONFIG['TOTAL_FEATURE_LIST'])
SELECT_FEATURE_LIST = list(DEFAULT_CONFIG['SELECT_FEATURE_LIST'])
RANDOM_STATE = DEFAULT_CONFIG['RANDOM_STATE']
N_JOBS=1

CONTEXT_SETTINGS = {'help_option_names': ['-h', '--help']}
@click.command(context_settings=CONTEXT_SETTINGS)
@click.argument('viral-ftr', type=click.Path(resolve_path=True))
@click.argument('nonviral-ftr', type=click.Path(resolve_path=True))
@click.option(
    '--balanced', 
    is_flag=True, 
    default=False,
    show_default=True, 
    type=bool,
    help=('random undersample the larger to the size of the smaller '
        'feature file')
)
@click.option(
    '-j',
    '--jobs',
    default=multiprocessing.cpu_count(),
    type=int,
    show_default=True,
    help='number of threads for training classifier',
)
def main(viral_ftr, nonviral_ftr, balanced, jobs):
    df_v = pd.read_csv(viral_ftr, sep='\t', header=0,)
    df_nv = pd.read_csv(nonviral_ftr, sep='\t', header=0,)



    if balanced:
        min_size = min(len(df_v), len(df_nv))
        if len(df_v) != min_size:
            df_v = df_v.sample(min_size, replace=False)
        if len(df_nv) != min_size:
            df_nv = df_nv.sample(min_size, 
                    replace=False, random_state=RANDOM_STATE)
    df_v['class'] = 1
    df_nv['class'] = 0



    df = pd.concat([df_v, df_nv])
    y = df['class']
    X = df.loc[:, SELECT_FEATURE_LIST]

    clf_to_train = RandomForestClassifier(random_state=RANDOM_STATE, 
            n_jobs=1)

    ## option 1: pipe on gridsearch 
    #params_to_tune = [{'n_estimators':[20, 50, 100, 150, 200], 'criterion':['gini','entropy']}]
    #gs = GridSearchCV(clf_to_train, params_to_tune, 
    #        scoring='f1', n_jobs=jobs)
    #steps = [('scaler', MinMaxScaler()), ('gs', gs)]
    #pipe = Pipeline(steps)
    #pipe.fit(X, y)
    #gs_clf = pipe.named_steps.gs
    #rf_clf = gs_clf.best_estimator_
    #best_params = gs_clf.best_params_
    #model = pipe

    #optino 2: gridsearch on pipe
    steps = [('scaler', MinMaxScaler()), ('rf', clf_to_train)]
    pipe = Pipeline(steps)
    params_to_tune = [{'rf__n_estimators':[20, 50, 100, 150, 200], 'rf__criterion':['gini','entropy']}]
    gs = GridSearchCV(pipe, params_to_tune, scoring='f1', n_jobs=jobs)
    gs_clf = gs.fit(X, y)
    best_params = gs_clf.best_params_
    best_pipe = gs_clf.best_estimator_
    rf_clf = best_pipe.named_steps.rf
    model = best_pipe

    #print(best_params)

    # save pipe
    joblib.dump(model, 'model')

    ### Note that not all clf have feature_importances_ attr
    feature_importances = pd.Series(rf_clf.feature_importances_, 
            index=SELECT_FEATURE_LIST)
    feature_importances = feature_importances.sort_values(ascending=False)
    #print(feature_importances)
    feature_importances.to_csv('feature-importances.tsv', 
            sep='\t', index=True, index_label='feature', 
            header=['importance'])


if __name__ == '__main__':
    main()
