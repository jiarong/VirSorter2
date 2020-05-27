import sys
import os
import logging
import multiprocessing
import subprocess
import click

from virsorter import __version__
from virsorter.config import (set_logger, make_config, 
        load_configfile, validate_config)

set_logger()

def log_exception(msg):
    logging.critical(msg)
    logging.info("Documentation is available at: https://project.readthedocs.io")
    logging.info("Issues can be raised at: https://github.com/project/project/issues")
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
    type=click.Path(dir_okay=True,writable=True,resolve_path=True),
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
    help='classifiers of viral groups to included'
)
@click.option(
    '-j',
    '--jobs',
    default=multiprocessing.cpu_count(),
    type=int,
    show_default=True,
    help='use at most this many jobs in parallel (see cluster submission for mor details).',
)
@click.option(
    '--min-score',
    default=0.5,
    type=float,
    show_default=True,
    help='minimal score to be identified as viral',
)
@click.option(
    '--provirus',
    default=False,
    is_flag=True,
    show_default=True,
    help='To include extracting provirus after classifying full contig seqs; Turning this on can significantly reduce the speed',
)
@click.option(
    '--hallmark-required-on-short',
    default=False,
    is_flag=True,
    show_default=True,
    help='require hallmark gene short seqs (length cutoff were set in template-config.yaml file)',
)
@click.option(
    '--max-orf-per-seq',
    default=20,
    type=int,
    show_default=True,
    help='Max # of orf used for computing taxonomic features; if # of orf in a seq exceeds the max limit, it is sub-sampled to this # to reduce computation; to turn off this, set it to -1; this option is ignored when --provirus is enabled'
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
@click.argument(
    'snakemake_args', 
    nargs=-1, 
    type=click.UNPROCESSED, 
)
def run_workflow(workflow, working_dir, db_dir, seqfile, include_groups, jobs,  min_score, provirus, hallmark_required_on_short, max_orf_per_seq, min_length, tmpdir, verbose, profile, dryrun, snakemake_args):
    ''' Runs the virsorter main function: to classify viral sequences

    By default all steps are executed. The "classify" rerun classify
    step without previous steps that are computationally heavy. Most
    snakemake arguments can be appended to the command for more info see
    'snakemake --help'. For more details, see: github
    '''
    os.makedirs(working_dir, exist_ok=True)
    config_f = os.path.join(working_dir,'config.yaml')

    if min_score > 1 or min_score < 0:
        logging.critical('--min-score needs to be between 0 and 1')
    if min_length < 0:
        logging.critical('--min-length needs to be >= 0')
    if jobs < 0:
        logging.critical('--jobs needs to be >= 0')

    if provirus:
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
    validate_config(config_f, workflow)
    config = load_configfile(config_f)

    if db_dir == None:
        db_dir = config['DBDIR']

    cmd = (
        'snakemake --snakefile {snakefile} --directory {working_dir} '
        '--jobs {jobs} '
        '--configfile {config_file} --conda-prefix {conda_prefix} '
        '--rerun-incomplete --use-conda --nolock --latency-wait 600'
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
        exit(1)

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
        '--directory {db_dir} '
        '--jobs {jobs} --rerun-incomplete '
        '--nolock  --use-conda  --conda-prefix {conda_prefix} '
        '{add_args} {args}'
    ).format(
        snakefile=get_snakefile('rules/setup.smk'),
        db_dir=db_dir,
        jobs=jobs,
        conda_prefix=os.path.join(db_dir,'conda_envs'),
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
