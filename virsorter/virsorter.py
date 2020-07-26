import sys
import os
import logging
import multiprocessing
import subprocess
import glob
import shutil
import click

from snakemake import load_configfile
from virsorter import __version__
from virsorter.config import get_default_config, set_logger, make_config

set_logger()

def log_exception(msg):
    logging.critical(msg)
    logging.info("Documentation is available at: https://github.com/jiarong/VirSorter2")
    logging.info("Issues can be raised at: https://github.com/jiarong/VirSorter2/issues")
    sys.exit(1)

@click.group(context_settings=dict(help_option_names=["-h", "--help"]))
@click.version_option(__version__)
@click.pass_context
def cli(obj):
    """
    virsorter - workflow for identifying viral sequences
    """

#cli.command('train') # work on this later

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
    help='sequence file in fa or fq format (could be compressed by gzip or bz2)'
)
@click.option(
    '--include-groups',
    default='dsDNAphage,NCLDV,RNA,ssDNA,lavidaviridae',
    type=str,
    show_default=True,
    help='classifiers of viral groups to included (comma separated and no space in between)'
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
    '--hallmark-required-on-short',
    default=False,
    is_flag=True,
    show_default=True,
    help='require hallmark gene short seqs (length cutoff as "short" were set in template-config.yaml file); this can reduce false positives at reasonable cost of sensitivity',
)
@click.option(
    '--provirus-off',
    default=False,
    is_flag=True,
    show_default=True,
    help='To turn off extracting provirus after classifying full contig seqs; Togetehr with lower --max-orf-per-seq, can speed up a run significantly',
)
@click.option(
    '--max-orf-per-seq',
    default=20,
    type=int,
    show_default=True,
    help='Max # of orf used for computing taxonomic features; if # of orf in a seq exceeds the max limit, it is sub-sampled to this # to reduce computation; to turn off this, set it to -1; this option must be used together with --provirus-off option'
)
@click.option(
    '--min-length',
    default=0,
    type=int,
    show_default=True,
    help='minimal seq length required; all seqs shorter than this will be removed',
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
    help='Stop using the conda envs (vs2.yaml) that come with this package and use what are installed in current system; Only useful when you want to install dependencies on your own with your own prefer versions',
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
def run_workflow(workflow, working_dir, db_dir, seqfile, include_groups, jobs,  min_score, provirus_off, hallmark_required_on_short, max_orf_per_seq, min_length, tmpdir, rm_tmpdir, verbose, profile, dryrun, use_conda_off, snakemake_args):
    ''' Runs the virsorter main function: to classify viral sequences

    By default all steps are executed. The "classify" rerun classify
    step without previous steps that are computationally heavy. Most
    snakemake arguments can be appended to the command for more info see
    'snakemake --help'. For more details, see: github
    '''

    # hard coded, need to change all "iter-0" to Tmpdir in smk
    tmpdir = 'iter-0'

    os.makedirs(working_dir, exist_ok=True)
    config_f = os.path.join(working_dir,'config.yaml')

    if min_score > 1 or min_score < 0:
        logging.critical('--min-score needs to be between 0 and 1')
    if min_length < 0:
        logging.critical('--min-length needs to be >= 0')
    if jobs < 0:
        logging.critical('--jobs needs to be >= 0')

    if provirus_off:
        provirus = False
    else:
        provirus = True
        max_orf_per_seq = -1

    if workflow == 'classify':
        target_f = '{working_dir}/{tmpdir}/all-fullseq-proba.tsv'.format(
                working_dir=working_dir,
                tmpdir=tmpdir,
        )
        try:
            subprocess.run(['touch', target_f], check=True)
        except subprocess.CalledProcessError as e:
            # removes the traceback
            logging.critical(e)

    make_config(
            db_dir=db_dir, seqfile=seqfile, include_groups=include_groups,
            threads=jobs, config_f=config_f, provirus=provirus,
            hallmark_required_on_short=hallmark_required_on_short,
            max_orf_per_seq=max_orf_per_seq, 
            tmpdir=tmpdir, min_length=min_length, min_score=min_score,
    )
    config = load_configfile(config_f)

    if db_dir == None:
        db_dir = config['DBDIR']

    cmd = (
        'snakemake --snakefile {snakefile} --directory {working_dir} '
        '--jobs {jobs} '
        '--configfile {config_file} --conda-prefix {conda_prefix} '
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
        conda_prefix=os.path.join(db_dir,'conda_envs')
    )
    logging.info('Executing: %s' % cmd)
    try:
        subprocess.run(cmd, check=True, shell=True)
    except subprocess.CalledProcessError as e:
        # removes the traceback
        logging.critical(e)
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
@click.argument('snakemake_args', nargs=-1, type=click.UNPROCESSED)
def run_setup(db_dir,jobs, snakemake_args):
    '''Setup databases and install dependencies.
    
    Executes a snakemake workflow to download reference database files
    and validate based on their MD5 checksum, and install dependencies
    to {conda_prefix} path.
    '''
    cmd = (
        'snakemake --snakefile {snakefile} '
        '--directory {db_dir} --quiet '
        '--jobs {jobs} --rerun-incomplete --latency-wait 600 '
        '--nolock  --use-conda --conda-prefix {conda_prefix} '
        '{add_args} {args}'
    )
    cmd_str = cmd.format(
        snakefile=get_snakefile('rules/setup.smk'),
        db_dir=db_dir,
        jobs=jobs,
        conda_prefix=os.path.join(db_dir,'conda_envs'),
        add_args='' if snakemake_args and snakemake_args[0].startswith('-') else '--',
        args=' '.join(snakemake_args),
    )

    logging.info('Setting up VirSorter2 database; this might takes ~10 mins and only needs to be done once.')
    #logging.info('Executing: %s' % cmd_str)
    try:
        subprocess.run(cmd_str, check=True, shell=True)
    except subprocess.CalledProcessError as e:
        # removes the traceback
        #logging.critical(e)
        #retry
        try: 
            logging.info('Setting up VirSorter2 database failed due to server not responding or bad internet connection; retrying..; please be patient.')
            cmd_str = cmd.format(
                snakefile=get_snakefile('rules/setup-retry.smk'),
                db_dir=db_dir,
                jobs=jobs,
                conda_prefix=os.path.join(db_dir,'conda_envs'),
                add_args='' if snakemake_args and snakemake_args[0].startswith('-') else '--',
                args=' '.join(snakemake_args),
            )
            subprocess.run(cmd_str, check=True, shell=True)
        except subprocess.CalledProcessError as e:
            # remove the traceback
            #logging.critical(e)
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
    help='genome sequence file for training; for file pattern globbing, put quotes around the pattern, eg. "fasta-dir/*.fa"',
    type=str,
    required=True,
    multiple=True,
)
@click.option(
    '--hmm',
    help='customized viral HMMs for training; default to the one used in VirSorter2',
    type=click.Path(resolve_path=True),
)
@click.option(
    '--hallmark',
    help='hallmark gene hmm list from -hmm for training (a tab separated file with three columns: 1. hmm name 2. gene name of hmm 3. hmm bit score cutoff); default to the one used for dsDNAphage in VirSorter2',
    type=click.Path(resolve_path=True),
)
@click.option(
    '--prodigal-train',
    help='customized training db from prodigal; default to the one used in prodigal --meta mode',
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
    help='Max # of orf used for computing taxonomic features; if # of orf in a seq exceeds the max limit, it is sub-sampled to this # to reduce computation; to turn off this, set it to -1'
)
@click.option(
    '--genome-as-bin',
    default=False,
    is_flag=True,
    show_default=True,
    help='if applied, each file (genome bin) is a genome in --seqfile, else each sequence is a genome',
)
@click.option(
    '--use-conda-off',
    is_flag=True,
    default=False,
    show_default=True,
    help='Stop using the conda envs (vs2.yaml) that come with this package and use what are installed in current system; Only useful when you want to install dependencies on your own with your own prefer versions',
)
@click.argument('snakemake_args', nargs=-1, type=click.UNPROCESSED)
def train_feature(working_dir, seqfile, hmm, hallmark, prodigal_train, frags_per_genome, min_length, max_orf_per_seq, genome_as_bin, jobs, use_conda_off, snakemake_args):
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
        '--nolock  {use_conda_off} --quiet --conda-prefix {conda_prefix} '
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
        conda_prefix=os.path.join(DEFAULT_CONFIG['DBDIR'],'conda_envs'),
        add_args='' if snakemake_args and snakemake_args[0].startswith('-') else '--',
        args=' '.join(snakemake_args),
    )
    logging.info('Executing: %s' % cmd)
    try:
        subprocess.run(cmd, check=True, shell=True)
    except subprocess.CalledProcessError as e:
        # removes the traceback
        logging.critical(e)
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
    help='Stop using the conda envs (vs2.yaml) that come with this package and use what are installed in current system; Only useful when you want to install dependencies on your own with your own prefer versions',
)
@click.argument('snakemake_args', nargs=-1, type=click.UNPROCESSED)
def train_model(working_dir, viral_ftrfile, nonviral_ftrfile, balanced, jobs, use_conda_off, snakemake_args):
    '''Training customized classifier model.
    '''

    DEFAULT_CONFIG = get_default_config

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
        '--nolock {use_conda_off} --quiet --conda-prefix {conda_prefix} '
        '{add_args} {args}'
    ).format(
        snakefile=get_snakefile('rules/train-model.smk'),
        working_dir=working_dir,
        viral_ftrfile=viral_ftrfile,
        nonviral_ftrfile=nonviral_ftrfile,
        balanced=balanced,
        jobs=jobs,
        use_conda_off='' if use_conda_off else '--use-conda',
        conda_prefix=os.path.join(DEFAULT_CONFIG['DBDIR'],'conda_envs'),
        add_args='' if snakemake_args and snakemake_args[0].startswith('-') else '--',
        args=' '.join(snakemake_args),
    )
    logging.info('Executing: %s' % cmd)
    try:
        subprocess.run(cmd, check=True, shell=True)
    except subprocess.CalledProcessError as e:
        # removes the traceback
        logging.critical(e)
        exit(1)
