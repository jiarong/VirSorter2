import sys
import os
import multiprocessing
import logging

from ruamel.yaml import YAML

from snakemake import utils
from snakemake.io import load_configfile

user_config_dir = os.path.join(os.path.expanduser('~'), '.virsorter')
user_template = os.path.join(user_config_dir, 'template-config.yaml')
src_config_dir = os.path.dirname(os.path.abspath(__file__))
src_template = os.path.join(src_config_dir, 'template-config.yaml')

if os.path.isfile(user_template):
    # check .virsorter in user home direcory first
    template = user_template
    #os.makedirs(user_config_dir, exist_ok=True)
else:
    assert os.path.isfile(src_template)
    template = src_template

TEMPLATE = template
DEFAULT_CONFIG = YAML().load(open(TEMPLATE))

def set_logger():
    logging.basicConfig(
        level=logging.INFO,
        datefmt="%Y-%m-%d %H:%M",
        format="[%(asctime)s %(levelname)s] %(message)s",
    )

### functions needes:
# make_config
# load_config
# validate_config
# No make_sample_table function is needed for virsorter

def make_config(db_dir, seqfile, config_f, include_groups, tmpdir, provirus=True, hallmark_required_on_short=False, threads=None, max_orf_per_seq=20):
    '''
    read config params from template-config.yaml
    then update the some params provided by command line

    Args:
        dbdir (str): location of database
        threads (int): number of threads per node to utilize
        config (str): path of output yaml file
    '''

    yaml = YAML()
    #yaml.version = (1, 1) # default is (1, 2)
    #yaml.default_flow_style = False # only needed for pre 1.2 version
    with open(TEMPLATE) as fp:
        config = yaml.load(fp)

    if db_dir != None:
        config['DBDIR'] = db_dir
    config['SEQFILE'] = seqfile
    config['PROVIRUS'] = provirus
    config['HALLMARK_REQUIRED_ON_SHORT'] = hallmark_required_on_short
    config['MAX_ORF_PER_SEQ'] = max_orf_per_seq
    config['TMPDIR'] = tmpdir

    config['THREADS'] = multiprocessing.cpu_count() if not threads else threads

    groups = [i.strip() for i in include_groups.split(',')]
    config['GROUPS'] = groups

    with open(config_f, 'w') as fw:
        yaml.dump(config, fw)
    mes = 'conig file written to {}\n'
    logging.info(mes.format(config_f))

def validate_config(config_f, workflow):
    config = load_configfile(config_f)
    # validate_sample_defs(config, workflow)
    # low priority; could later add more validation steps
