import sys
import os
import logging
import multiprocessing
import subprocess
import glob
import shutil
import click

from snakemake import load_configfile
from ruamel.yaml import YAML
from virsorter import __version__
from virsorter.config import get_default_config, set_logger, make_config

set_logger()
logging.info(f'VirSorter {__version__}')

def log_exception(msg):
    logging.critical(msg)
    logging.info('Documentation is available at: '
                    'https://github.com/jiarong/VirSorter2')
    logging.info('Issues can be raised at: '
                    'https://github.com/jiarong/VirSorter2/issues')
    sys.exit(1)

@click.group(context_settings=dict(help_option_names=["-h", "--help"]))
@click.version_option(__version__)
@click.pass_context
def cli(obj):
    """
    virsorter - workflow for identifying viral sequences
    """

def get_snakefile(f="Snakefile"):
    sf = os.path.join(os.path.dirname(os.path.abspath(__file__)), f)
    if not os.path.exists(sf):
        sys.exit("Unable to locate the Snakemake workflow file; tried %s" % sf)
    return sf

## run command

@cli.command(
    'run',
    context_settings=dict(ignore_unknown_options=True),
    short_help='run virsorter main workflow'
)
@click.argument(
    'workflow',
    default='all',
    type=click.Choice(['all', 'classify']),
#    show_default=True,
#    help='Execute only subworkflow.',
)
@click.option('-w',
    '--working-dir',
    type=click.Path(dir_okay=True,writable=True,resolve_path=True),
    help='output directory',
    default='.'
)
@click.option('-d',
    '--db-dir',
    required=False,
    type=click.Path(dir_okay=True,writable=False,resolve_path=True),
    help='database directory, default to the --db-dir set during installation',
)
@click.option('-i',
    '--seqfile',
    required=True,
    type=click.Path(resolve_path=True),
    help=('sequence file in fa or fq format (could be compressed '
            'by gzip or bz2)')
)
@click.option(
    '-l',
    '--label',
    default='',
    help=('add label as prefix to output files; this is useful when re-running '
            'classify with different filtering options'),
)
@click.option(
    '--include-groups',
    default='dsDNAphage,NCLDV,RNA,ssDNA,lavidaviridae',
    type=str,
    show_default=True,
    help=('classifiers of viral groups to included '
        '(comma separated and no space in between); available options are: '
        'dsDNAphage,NCLDV,RNA,ssDNA,lavidaviridae'),
)
@click.option(
    '-j',
    '--jobs',
    default=multiprocessing.cpu_count(),
    type=int,
    show_default=True,
    help='max # of jobs allowed in parallel.',
)
@click.option(
    '--min-score',
    default=0.5,
    type=float,
    show_default=True,
    help='minimal score to be identified as viral',
)
@click.option(
    '--hallmark-required',
    default=False,
    is_flag=True,
    show_default=True,
    help='require hallmark gene on all seqs',
)
@click.option(
    '--hallmark-required-on-short',
    default=False,
    is_flag=True,
    show_default=True,
    help=('require hallmark gene on short seqs (length cutoff '
            'as "short" were set by "MIN_SIZE_ALLOWED_WO_HALLMARK_GENE" in '
            'template-config.yaml file, default 3kbp); this can reduce '
            'false positives at reasonable cost of sensitivity'),
)
@click.option(
    '--viral-gene-required',
    default=False,
    is_flag=True,
    show_default=True,
    help=('requires viral genes annotated, removing putative viral seqs '
            'with no genes annotated; this can reduce false positives at '
            'reasonable cost of sensitivity'),
)
@click.option(
    '--provirus-off',
    default=False,
    is_flag=True,
    show_default=True,
    help=('To turn off extracting provirus after classifying full contigs; '
            'Togetehr with lower --max-orf-per-seq, can speed up a run '
            'significantly'),
)
@click.option(
    '--max-orf-per-seq',
    default=-1,
    type=int,
    show_default=True,
    help=('Max # of orf used for computing taxonomic feature; this option '
            'can only be used in --provirus-off mode; if # of orf in a seq '
            'exceeds the max limit, it is sub-sampled to this # '
            'to reduce computation')
)
@click.option(
    '--min-length',
    default=0,
    type=int,
    show_default=True,
    help=('minimal seq length required; all seqs shorter than this will '
            'be removed'),
)
@click.option(
    '--prep-for-dramv',
    default=False,
    is_flag=True,
    show_default=True,
    help=('To turn on generating viral seqfile and viral-affi-contigs.tab '
            'for DRAMv'),
)
@click.option(
    '--tmpdir',
    default='iter-0',
    help='Directory name for intermediate files',
)
@click.option(
    '--verbose',
    is_flag=True,
    default=False,
    show_default=True,
    help='shows detailed rules output',
)
@click.option(
    '--profile',
    default=None,
    help='snakemake profile e.g. for cluster execution.',
)
@click.option(
    '-n',
    '--dryrun',
    is_flag=True,
    default=False,
    show_default=True,
    help='Check rules to run and files to produce',
)
@click.option(
    '--use-conda-off',
    is_flag=True,
    default=False,
    show_default=True,
    help=('Stop using the conda envs (vs2.yaml) that come with this '
            'package and use what are installed in current system; '
            'Only useful when you want to install dependencies on your own '
            'with your own prefer versions'),
)
@click.option(
    '--rm-tmpdir',
    is_flag=True,
    default=False,
    show_default=True,
    help='Remove intermediate file directory (--tmpdir)'
)
@click.argument(
    'snakemake_args', 
    nargs=-1, 
    type=click.UNPROCESSED, 
)
def run_workflow(workflow, working_dir, db_dir, seqfile, include_groups,
        jobs,  min_score, hallmark_required, hallmark_required_on_short,
        viral_gene_required, provirus_off, max_orf_per_seq, min_length,
        prep_for_dramv, tmpdir, rm_tmpdir, verbose, profile, dryrun,
        use_conda_off, snakemake_args, label):
    ''' Runs the virsorter main function to classify viral sequences

    This includes 3 steps: 1) preprocess, 2) feature extraction, and 3)
    classify. By default ("all") all steps are executed. The "classify"
    only run the 3) classify step without previous steps that are
    computationally heavy, good for rerunning with different filtering
    options (--min-score, --hallmark-required,
    --hallmark-required-on-short, --viral-gene-required). Most snakemake
    arguments can be appended to the command for more info see
    'snakemake --help'.
    '''

    # hard coded, need to change all "iter-0" to Tmpdir in smk
    tmpdir = 'iter-0'

    os.makedirs(working_dir, exist_ok=True)
    config_f = os.path.join(working_dir,'config.yaml')

    if min_score > 1 or min_score < 0:
        logging.critical('--min-score needs to be between 0 and 1')
        sys.exit(1)
    if min_length < 0:
        logging.critical('--min-length needs to be >= 0')
        sys.exit(1)
    if jobs < 0:
        logging.critical('--jobs needs to be >= 0')
        sys.exit(1)

    if provirus_off:
        provirus = False
        if max_orf_per_seq != -1 and prep_for_dramv:
            mes = ('--max-orf-per-seq CAN NOT be used with '
                    '--prep-for-dramv; '
                    'outputs with ORFs subsampled are NOT '
                    'compatible with DRAMv')
    else:
        provirus = True
        max_orf_per_seq = -1

    if workflow == 'classify':
        if not os.path.exists(config_f):
            mes = 'No config.yaml dectected from previous run'
            logging.critical(mes)
            sys.exit(1)

        config = load_configfile(config_f)
        min_length_prev = config['MIN_LENGTH']
        if min_length != min_length_prev:
            mes = (
                '--min-length has changed from '
                f'{min_length_prev} to {min_length}; '
                'but --min-length has not effect on classify step; '
                'The whole pipeline has to be rerun if --min-length changes'
            )
            logging.critical(mes)
            sys.exit(1)

        if provirus != config['PROVIRUS']:
            mes = (
                '--provirus-off setting change found; '
                'The whole pipeline has to be rerun if --provirus-off changes'
            )
            logging.critical(mes)
            sys.exit(1)

        target_f = '{working_dir}/{tmpdir}/reclassify.trigger'.format(
                working_dir=working_dir,
                tmpdir=tmpdir,
        )
        try:
            subprocess.run(['touch', target_f], check=True)
        except subprocess.CalledProcessError as e:
            # removes the traceback
            #logging.critical(e)
            sys.exit(1)

        os.rename(config_f, os.path.join(working_dir,'config.yaml.bak'))

    make_config(
            db_dir=db_dir, seqfile=seqfile, include_groups=include_groups,
            threads=jobs, config_f=config_f, provirus=provirus,
            hallmark_required=hallmark_required,
            hallmark_required_on_short=hallmark_required_on_short,
            viral_gene_required=viral_gene_required,
            prep_for_dramv=prep_for_dramv,
            max_orf_per_seq=max_orf_per_seq, 
            tmpdir=tmpdir, min_length=min_length, min_score=min_score, 
            label=label,
    )
    config = load_configfile(config_f)

    if db_dir == None:
        db_dir = config['DBDIR']

    cmd = (
        'snakemake --snakefile {snakefile} --directory {working_dir} '
        '--jobs {jobs} '
        '--configfile {config_file} {conda_prefix} '
        '--rerun-incomplete {use_conda_off} --nolock --latency-wait 600'
        ' {profile} {dryrun} {verbose} '
        ' {target_rule} '
        ' {args} '
    ).format(
        snakefile=get_snakefile(),
        working_dir=working_dir,
        jobs=jobs,
        config_file=config_f,
        profile='' if (profile is None) else '--profile {}'.format(profile),
        dryrun='--dryrun' if dryrun else '',
        use_conda_off='' if use_conda_off else '--use-conda',
        verbose='' if verbose else '--quiet',
        args=' '.join(snakemake_args),
        target_rule='-R {}'.format(workflow) if workflow!='all' else workflow,
        conda_prefix='' if use_conda_off else '--conda-prefix {}'.format(
            os.path.join(db_dir,'conda_envs')
        )
    )
    logging.info('Executing: %s' % cmd)
    try:
        subprocess.run(cmd, check=True, shell=True)
    except subprocess.CalledProcessError as e:
        # removes the traceback
        #logging.critical(e)
        sys.exit(1)

    if rm_tmpdir:
        tmpdir_path = os.path.join(working_dir, tmpdir)
        shutil.rmtree(tmpdir_path, ignore_errors=True)

# initialize
@cli.command(
    'setup',
    context_settings=dict(ignore_unknown_options=True),
    short_help='download reference files (~10GB) and install dependencies',
)
@click.option('-d',
    '--db-dir',
    help='diretory path for databases',
    type=click.Path(dir_okay=True,writable=True,resolve_path=True),
    required=True
)
@click.option(
    '-j',
    '--jobs',
    default=multiprocessing.cpu_count(),
    type=int,
    show_default=True,
    help='number of simultaneous downloads',
)
@click.option(
    '-s',
    '--skip-deps-install',
    is_flag=True,
    show_default=True,
    help=('skip dependency installation (if you want to install on your '
            'own as in development version)'),
)
@click.argument('snakemake_args', nargs=-1, type=click.UNPROCESSED)
def run_setup(db_dir,jobs, skip_deps_install, snakemake_args):
    '''Setup databases and install dependencies.
    
    Executes a snakemake workflow to download reference database files
    and validate based on their MD5 checksum, and install dependencies
    '''
    cmd = (
        'snakemake --snakefile {snakefile} '
        '--directory {db_dir} --quiet '
        '--config Skip_deps_install={skip_deps_install} '
        '--jobs {jobs} --rerun-incomplete --latency-wait 600 '
        '--nolock  --use-conda --conda-prefix {conda_prefix} '
        '{args}'
    )
    cmd_str = cmd.format(
        snakefile=get_snakefile('rules/setup.smk'),
        db_dir=db_dir,
        skip_deps_install=skip_deps_install,
        jobs=jobs,
        conda_prefix=os.path.join(db_dir,'conda_envs'),
        args=' '.join(snakemake_args),
    )

    logging.info('Setting up VirSorter2 database; this might take ~10 mins '
                    'and only needs to be done once.')
    #logging.info('Executing: %s' % cmd_str)
    # try first with zenodo
    try:
        subprocess.run(
                cmd_str, check=True, shell=True,
                #stderr=subprocess.PIPE, 
                #stdout=subprocess.PIPE,
        )
    except subprocess.CalledProcessError as e:
        logging.info('First attempt failed; trying the second time.')
        # try again with osf
        try: 
            cmd_str = cmd.format(
                snakefile=get_snakefile('rules/setup-retry.smk'),
                db_dir=db_dir,
                skip_deps_install=skip_deps_install,
                jobs=jobs,
                conda_prefix=os.path.join(db_dir,'conda_envs'),
                args=' '.join(snakemake_args),
            )
            subprocess.run(cmd_str, check=True, shell=True)
        except subprocess.CalledProcessError as e:
            # remove the traceback
            logging.critical(e)
            sys.exit(1)

# train feature
@cli.command(
    'train-feature',
    context_settings=dict(ignore_unknown_options=True),
    short_help='subcommand for training feature of customized classifier',
)
@click.option('-w',
    '--working-dir',
    help='output directory',
    type=click.Path(dir_okay=True,writable=True,resolve_path=True),
    required=True
)
@click.option(
    '--seqfile',
    help=('genome sequence file for training; for file pattern globbing, '
            'put quotes around the pattern, eg. "fasta-dir/*.fa"'),
    type=str,
    required=True,
    multiple=True,
)
@click.option(
    '--hmm',
    help=('customized viral HMMs for training; default to the one used in '
            'VirSorter2'),
    type=click.Path(resolve_path=True),
)
@click.option(
    '--hallmark',
    type=click.Path(resolve_path=True),
    help=('hallmark gene hmm list from -hmm for training (a tab separated '
            'file with three columns: 1. hmm name, 2. gene name of hmm, '
            '3. hmm bit score cutoff); default to the one used for '
            'dsDNAphage in VirSorter2'),
)
@click.option(
    '--prodigal-train',
    help=('customized training db from prodigal; default to the one used '
            'in prodigal --meta mode'),
    type=click.Path(resolve_path=True),
)
@click.option(
    '--frags-per-genome',
    default=5,
    type=int,
    show_default=True,
    help='number of random DNA fragments collected from each genome',
)
@click.option(
    '-j',
    '--jobs',
    default=multiprocessing.cpu_count(),
    type=int,
    show_default=True,
    help='max # of jobs in parallel',
)
@click.option(
    '--min-length',
    default=1000,
    type=int,
    show_default=True,
    help='minimum size of random DNA fragment for training',
)
@click.option(
    '--max-orf-per-seq',
    default=20,
    type=int,
    show_default=True,
    help=('Max # of orf used for computing taxonomic features; if # of orf '
            'in a seq exceeds the max limit, it is sub-sampled to this # '
            'to reduce computation; to turn off this, set it to -1'),
)
@click.option(
    '--genome-as-bin',
    default=False,
    is_flag=True,
    show_default=True,
    help=('if applied, each file (genome bin) is a genome in --seqfile, '
            'else each sequence is a genome'),
)
@click.option(
    '--use-conda-off',
    is_flag=True,
    default=False,
    show_default=True,
    help=('Stop using the conda envs (vs2.yaml) that come with this package '
            'and use what are installed in current system; Only useful when '
            'you want to install dependencies on your own with your own '
            'preferred versions'),
)
@click.argument('snakemake_args', nargs=-1, type=click.UNPROCESSED)
def train_feature(working_dir, seqfile, hmm, hallmark, prodigal_train, 
                    frags_per_genome, min_length, max_orf_per_seq, 
                    genome_as_bin, jobs, use_conda_off, snakemake_args):
    '''Training features for customized classifier.
    
    Executes a snakemake workflow to do the following:
    1) prepare random DNA fragments from viral and nonviral genome data 
    2) extract feature from random DNA fragments to make ftrfile
    '''

    DEFAULT_CONFIG = get_default_config()

    cwd = os.getcwd()
    lis = []
    pat_lis = []
    for pat in seqfile:
        # only works in linux
        if pat.startswith('/'):
            new_pat = pat
        else:
            new_pat = '{}/{}'.format(cwd, pat)
        fs = glob.glob(pat)
        lis.extend(fs)
        pat_lis.append(new_pat)
    
    if len(lis) == 0:
        mes = 'No files match {}'.format(viral_seqfile)
        logging.critical(mes)
        sys.exit(1)
    else:
        mes = '{} seqfiles are used for training features'.format(len(lis))
        logging.info(mes)

    if hmm == None:
        hmm = 'NA'
    if hallmark == None:
        hallmark = 'NA'

    if prodigal_train == None:
        prodigal_train = 'NA'

    cmd = (
        'snakemake --snakefile {snakefile} '
        '--directory {working_dir} '
        '--config Viral_seqfile="{seqfile}" '
            'Hmm={hmm} '
            'Hallmark={hallmark} '
            'Rbs={prodigal_train} '
            'Min_length={min_length} '
            'Max_orf_per_seq={max_orf_per_seq} '
            'Viral_genome_as_bin={genome_as_bin} '
            'Fragments_per_genome={frags_per_genome} '
        '--jobs {jobs} --rerun-incomplete --latency-wait 600 '
        '--nolock  {use_conda_off} --quiet {conda_prefix} '
        '{add_args} {args}'
    ).format(
        snakefile=get_snakefile('rules/train-feature.smk'),
        working_dir=working_dir,
        seqfile=' '.join(pat_lis),
        hmm=hmm,
        hallmark=hallmark,
        prodigal_train=prodigal_train,
        min_length=min_length,
        max_orf_per_seq=max_orf_per_seq,
        genome_as_bin=genome_as_bin,
        frags_per_genome=frags_per_genome, 
        jobs=jobs,
        use_conda_off='' if use_conda_off else '--use-conda',
        conda_prefix='' if use_conda_off else '--conda-prefix {}'.format(
            os.path.join(DEFAULT_CONFIG['DBDIR'],'conda_envs')
        ),
        add_args=('' if snakemake_args and snakemake_args[0].startswith('-') 
                    else '--'),
        args=' '.join(snakemake_args),
    )
    logging.info('Executing: %s' % cmd)
    try:
        subprocess.run(cmd, check=True, shell=True)
    except subprocess.CalledProcessError as e:
        # removes the traceback
        #logging.critical(e)
        exit(1)


# train model
@cli.command(
    'train-model',
    context_settings=dict(ignore_unknown_options=True),
    short_help='subcommand for training customized classifier model',
)
@click.option('-w',
    '--working-dir',
    help='output directory',
    type=click.Path(dir_okay=True,writable=True,resolve_path=True),
    required=True
)
@click.option(
    '--viral-ftrfile',
    help='viral genome feature file for training',
    type=click.Path(resolve_path=True),
    required=True,
)
@click.option(
    '--nonviral-ftrfile',
    help='nonviral genome feature file for training',
    type=click.Path(resolve_path=True),
    required=True,
)
@click.option(
    '-j',
    '--jobs',
    default=multiprocessing.cpu_count(),
    type=int,
    show_default=True,
    help='number of threads for classier',
)
@click.option(
    '--balanced',
    is_flag=True,
    default=False,
    show_default=True,
    type=bool,
    help='random undersample the larger to the size of the smaller feature file'
)
@click.option(
    '--use-conda-off',
    is_flag=True,
    default=False,
    show_default=True,
    help=('Stop using the conda envs (vs2.yaml) that come with this package '
            'and use what are installed in current system; Only useful when '
            'you want to install dependencies on your own with your own '
            'prefer versions'),
)
@click.argument('snakemake_args', nargs=-1, type=click.UNPROCESSED)
def train_model(working_dir, viral_ftrfile, nonviral_ftrfile, balanced, 
                jobs, use_conda_off, snakemake_args):
    '''Training customized classifier model.
    '''

    DEFAULT_CONFIG = get_default_config()

    if balanced == None:
        balanced = False
    cmd = (
        'snakemake --snakefile {snakefile} '
        '--directory {working_dir} '
        '--config '
            'Viral_ftrfile={viral_ftrfile} '
            'Nonviral_ftrfile={nonviral_ftrfile} '
            'Balanced={balanced} '
            'Jobs={jobs} '
        '--jobs {jobs} --rerun-incomplete --latency-wait 600 '
        '--nolock {use_conda_off} --quiet {conda_prefix} '
        '{add_args} {args}'
    ).format(
        snakefile=get_snakefile('rules/train-model.smk'),
        working_dir=working_dir,
        viral_ftrfile=viral_ftrfile,
        nonviral_ftrfile=nonviral_ftrfile,
        balanced=balanced,
        jobs=jobs,
        use_conda_off='' if use_conda_off else '--use-conda',
        conda_prefix='' if use_conda_off else '--conda-prefix {}'.format(
            os.path.join(DEFAULT_CONFIG['DBDIR'],'conda_envs')
        ),
        add_args=('' if snakemake_args and snakemake_args[0].startswith('-') 
                    else '--'),
        args=' '.join(snakemake_args),
    )
    logging.info('Executing: %s' % cmd)
    try:
        subprocess.run(cmd, check=True, shell=True)
    except subprocess.CalledProcessError as e:
        # removes the traceback
        #logging.critical(e)
        exit(1)


# config management
@cli.command(
    'config',
    context_settings=dict(ignore_unknown_options=True),
    #no_args_is_help=True,
    short_help='subcommand for configuration management',
)
@click.option(
    '--show',
    is_flag=True,
    default=False,
    show_default=True,
    help='show all configuration values',
)
@click.option(
    '--show-source',
    is_flag=True,
    default=False,
    show_default=True,
    help='show path of the configuration file',
)
@click.option(
    '--init-source',
    is_flag=True,
    default=False,
    show_default=True,
    help='initialize configuration file',
)
@click.option(
    '--db-dir',
    type=click.Path(),
    help='directory for databases; required for --init-source',
)
@click.option(
    '--set',
    help=('set KEY to VAL with the format: KEY=VAL; for nested dict in YAML '
            'use KEY1.KEY2=VAL (e.g. virsorter config --set '
            'GROUP_INFO.RNA.MIN_GENOME_SIZE=2000)'),
)
@click.option(
    '--get',
    help=('the value of a KEY (e.g. virsorter config '
            '--get GROUP_INFO.RNA.MIN_GENOME_SIZE'),
)
def config(show, show_source, init_source, db_dir, set, get):
    '''CLI for managing configurations.

    There are many configurations kept in "template-config.yaml" in source 
    code directory or "~/.virsorter" (when source code directory is not 
    writable for user). This file can located with 
    `virsorter config --show-source`. You can set the configurations with 
    `virsorter config --set KEY=VAL`. Alternative, you can edit in the 
    configuration file ("template-config.yaml") directly.
    '''

    from virsorter.config import (
            TEMPLATE, SRC_CONFIG_DIR, 
            USER_CONFIG_DIR, init_config_template
    )

    if init_source:
        if db_dir == None:
            mes = '--db-dir is required for --init-source'
            logging.critical(mes)
            sys.exit(1)
        else:
            if not os.path.isdir(db_dir):
                mes = (f'--db-dir {db_dir} does NOT exist yet; Make sure it '
                        'is created later\n')
                logging.warning(mes)

            init_config_template(SRC_CONFIG_DIR, USER_CONFIG_DIR, db_dir)
            sys.exit(0)

    if not os.path.isfile(TEMPLATE):
        mes = ('config file "template-config.yaml" has not been '
                'initialized yet; Please use '
                '`virsorter config --init-source --db-dir PATH` to initialize')
        logging.critical(mes)
        sys.exit(1)

    config = get_default_config()

    if show:
        YAML().dump(config, sys.stdout)
        sys.exit(0)

    if show_source:
        mes = f'config file path: {TEMPLATE}\n'
        sys.stdout.write(mes)
        sys.exit(0)

    if get != None:
        s = get
        lis = [var.strip() for var in s.split(',')]
        for var in lis:
            temp = config
            for i in var.split('.'):
                i = i.strip()
                try:
                    temp = temp[i]
                except KeyError as e:
                    mes = f'{i} is not a key in config file ({TEMPLATE})'
                    logging.critical(mes)
                    sys.exit(1)

            mes = f'{var}: {temp}\n'
            sys.stdout.write(mes)

        sys.exit(0)

    if set != None:
        s = set
        lis = [item.strip() for item in s.split(',')]
        for item in lis:
            temp = config
            var, val = item.split('=')
            var = var.strip()
            val = val.strip()
            keys = [key.strip() for key in var.split('.')]
            for i in range(len(keys)):
                if i == (len(keys) - 1):
                    # stop at 2nd last key
                    break
                key = keys[i]
                try:
                    temp = temp[key]
                except KeyError as e:
                    mes = f'{key} is not a key in config file ({TEMPLATE})'
                    logging.critical(mes)
                    sys.exit(1)

            last_key = keys[-1]

            try:
                old_val = temp[last_key]
                if os.path.exists(old_val):
                    val = os.path.abspath(val)
                temp[last_key] = val
            except KeyError as e:
                mes = f'{last_key} is not a key in config file ({TEMPLATE})'
                logging.critical(mes)
                sys.exit(1)

            mes = f'{var}: {old_val} ==> {val}\n'
            sys.stdout.write(mes)
            with open(TEMPLATE, 'w') as fw:
                YAML().dump(config, fw)

        sys.exit(0)
