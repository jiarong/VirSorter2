import sys
import os
import multiprocessing
import logging

from ruamel.yaml import YAML

#from snakemake.io import load_configfile

USER_CONFIG_DIR = os.path.join(os.path.expanduser('~'), '.virsorter')
SRC_CONFIG_DIR = os.path.dirname(os.path.abspath(__file__))

user_template = os.path.join(USER_CONFIG_DIR, 'template-config.yaml')
src_template = os.path.join(SRC_CONFIG_DIR, 'template-config.yaml')

if os.path.isfile(user_template):
    # check .virsorter in user home direcory first
    template = user_template
    #os.makedirs(user_config_dir, exist_ok=True)
else:
    template = src_template

TEMPLATE = template

def get_default_config():
    assert os.path.isfile(TEMPLATE)
    return YAML().load(open(TEMPLATE))

def set_logger():
    logging.basicConfig(
        level=logging.INFO,
        datefmt="%Y-%m-%d %H:%M",
        format="[%(asctime)s %(levelname)s] %(message)s",
    )

def init_config_template(src_config_dir, user_config_dir, db_dir):
    src_template_ori = os.path.join(src_config_dir, 
            'template-config-original.yaml')
    src_template = os.path.join(src_config_dir, 'template-config.yaml')
    user_template = os.path.join(user_config_dir, 'template-config.yaml')
    if os.access(src_template_ori, os.W_OK):
        # check .virsorter in user home direcory first
        template = src_template
        #os.makedirs(user_config_dir, exist_ok=True)
    else:
        os.makedirs(user_config_dir, exist_ok=True)
        mes = ('Attention: can not write template-config.yaml '
                'in source directory:\n'
                f'{src_config_dir}\n'
                'makeing a copy to user home direcotry:\n'
                f'{user_template}\n')

        logging.info(mes)
        template = user_template
        logging.info('Using {template} as config template')

    yaml = YAML()
    with open(src_template_ori) as fp:
        config = yaml.load(fp)
        config['DBDIR'] = db_dir
        logging.info(f'saving {db_dir} as DBDIR to config file {template}')

    with open(template, 'w') as fw:
        yaml.dump(config, fw)

    return config

### functions needes:
# make_config
# load_config
# validate_config
# No make_sample_table function is needed for virsorter

def make_config(db_dir, seqfile, config_f, include_groups, tmpdir,
        min_score=0.5, min_length=0, provirus=True, hallmark_required=False,
        hallmark_required_on_short=False, viral_gene_required=False,
        prep_for_dramv=False, threads=None, max_orf_per_seq=20, label='',):
    '''
    read config params from template-config.yaml
    then update the some params provided by command line

    Args:
        db_dir (str): location of database
        threads (int): number of threads per node to utilize
        config (str): path of output yaml file
    '''

    yaml = YAML()
    #yaml.version = (1, 1) # default is (1, 2)
    #yaml.default_flow_style = False # only needed for pre 1.2 version
    try:
        with open(TEMPLATE) as fp:
            config = yaml.load(fp)
            logging.info(f'Using {TEMPLATE} as config template')
            
        if TEMPLATE.startswith(USER_CONFIG_DIR):
            src_template_ori = os.path.join(SRC_CONFIG_DIR, 
                    'template-config-original.yaml')
            with open(src_template_ori) as fp:
                config_src = yaml.load(fp)
                st = set(config_src) - set(config)
                if len(st) != 0:
                    config.update(dict((i, config_src[i]) for i in st))
                    with open(TEMPLATE, 'w') as fw:
                        yaml.dump(config, fw)
                    logging.info(
                            f'{TEMPLATE} as config template does not '
                            'have all variables needed; '
                            'updating with those in {src_template_ori}')

    except FileNotFoundError as e:
        if db_dir != None:
            mes = ('"template-config.yaml" has not been initialized; '
                    'initializing..')
            logging.info(mes)
            config = init_config_template(SRC_CONFIG_DIR, 
                                             USER_CONFIG_DIR, db_dir)
        else:
            mes = ('--db-dir must be provided since "template-config.yaml" '
                    'has not been initialized')
            logging.critical(mes)
            sys.exit(1)

    if db_dir != None:
        config['DBDIR'] = db_dir
    db_dir = config['DBDIR']
    config['SEQFILE'] = seqfile
    config['PROVIRUS'] = provirus
    config['HALLMARK_REQUIRED'] = hallmark_required
    config['HALLMARK_REQUIRED_ON_SHORT'] = hallmark_required_on_short
    config['VIRAL_GENE_REQUIRED'] = viral_gene_required
    config['MAX_ORF_PER_SEQ'] = max_orf_per_seq
    config['TMPDIR'] = tmpdir
    config['PROBA_CUTOFF'] = min_score
    config['MIN_LENGTH'] = min_length
    config['PREP_FOR_DRAMV'] = prep_for_dramv
    config['LABEL'] = label

    config['THREADS'] = multiprocessing.cpu_count() if not threads else threads

    groups = [i.strip() for i in include_groups.split(',')]
    groups_avail = os.listdir(f'{db_dir}/group')
    groups_unavail = set(groups).difference(set(groups_avail))
    if len(groups_unavail) != 0:
        mes = (
                'Following viral groups are not available: {}\n'
                'Make sure viral group names match with direcotry '
                'under {}/group\n'
        )
        logging.critical(mes.format(', '.join(groups_unavail), db_dir))
    config['GROUPS'] = groups

    with open(config_f, 'w') as fw:
        yaml.dump(config, fw)
    mes = 'conig file written to {}\n'
    logging.info(mes.format(config_f))

def validate_config(config_f, workflow):
    pass
    #config = load_configfile(config_f)
    # validate_sample_defs(config, workflow)
    # low priority; could later add more validation steps
