import sys
import os
from ruamel.yaml import YAML
from snakemake.utils import min_version

### set minimum snakemake version ###
min_version('5.8.1')

Viral_ftrfile = config['Viral_ftrfile']
Nonviral_ftrfile = config['Nonviral_ftrfile'] 
Balanced = config['Balanced'] # true or false 
Jobs = config['Jobs'] # true or false 

# load other deault setting
# load template-config.yaml (not in the same dir as setup.smk)
#   need to go up 2 levels
src_config_dir = os.path.dirname(srcdir('.'))
src_config_dir = os.path.dirname(src_config_dir)

#print(srcdir('.'))
#print(workflow.basedir)
#print(workflow.snakefile)

src_template = os.path.join(src_config_dir, 'template-config.yaml')
yaml = YAML()
with open(src_template) as fp:
    config_default = yaml.load(fp)

Classify_threads = config_default['CLASSIFY_THREADS']
Proba_cutoff = config_default['PROBA_CUTOFF']

Contig_bp_per_split = config_default['CONTIG_BP_PER_SPLIT']
Faa_bp_per_split = config_default['FAA_BP_PER_SPLIT']
Gff_seqnum_per_split = config_default['GFF_SEQNUM_PER_SPLIT']

#Tmpdir=config['TMPDIR']
#Local_scratch=config['LOCAL_SCRATCH']
Tmpdir='iter-0'
Local_scratch='/tmp'

Srcdir = src_config_dir
Scriptdir='{}/scripts'.format(Srcdir)
Wkdir=os.path.abspath(os.getcwd())
Conda_yaml_dir = '../envs' # relative to rules/*.smk files

Groups_ftr = ['viral', 'nonviral']



wildcard_constraints:
    domain='[A-Za-z]+',

rule all:
    input: 'model'

rule make_classifer_model:
    input: 
        viral={Viral_ftrfile},
        nonviral={Nonviral_ftrfile},
    output: 'model'
    conda: f'{Conda_yaml_dir}/vs2.yaml'
    threads: Jobs
    shell:
        """
        if [ {Balanced} = "True" ]; then
            python {Scriptdir}/train-model.py {input.viral} {input.nonviral} --jobs {threads} --balanced
        else
            python {Scriptdir}/train-model.py {input.viral} {input.nonviral} --jobs {threads}
        fi
        """

container: "docker://continuumio/miniconda3"    

#onsuccess:
#    print('virsorter run finished')
#onerror:
#    print('Issues can be raised at: https://github.com/jiarong/VirSorter2/issues')



