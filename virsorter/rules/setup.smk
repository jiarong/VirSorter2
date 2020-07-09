import sys
import os
import shutil
import hashlib
import logging
from ruamel.yaml import YAML
from virsorter.config import DEFAULT_CONFIG, set_logger

set_logger()

ENV_YAML_DIR = '../envs'

# update DBDIR in template-config.yaml
user_config_dir = os.path.join(os.path.expanduser('~'), '.virsorter')
user_template = os.path.join(user_config_dir, 'template-config.yaml')

# following DONOT work in Snakefile!!
#src_config_dir = os.path.dirname(os.path.abspath(__file__)) 

# not in the same dir as setup.smk, need to go up 2 levels
src_config_dir = os.path.dirname(os.path.dirname(workflow.snakefile))
src_template_ori = os.path.join(src_config_dir, 'template-config-original.yaml')
src_template = os.path.join(src_config_dir, 'template-config.yaml')
Scriptdir=os.path.join(src_config_dir, 'scripts')

if os.access(src_template_ori, os.W_OK):
    # check .virsorter in user home direcory first
    template = src_template
    #os.makedirs(user_config_dir, exist_ok=True)
else:
    os.makedirs(user_config_dir, exist_ok=True)
    mes = ('Attention: can not write template-config.yaml in source directory:\n'
            f'{src_config_dir}\n'
            'makeing a copy to user home direcotry:\n'
            f'{user_template}\n')

    logging.info(mes)
    template = user_template

#print(srcdir('.'))
#print(workflow.basedir)
#print(workflow.snakefile)

dbdir = os.path.abspath(os.getcwd())
yaml = YAML()
with open(src_template_ori) as fp:
    config = yaml.load(fp)
    config['DBDIR'] = dbdir
    logging.info(f'saving to {dbdir} as DBDIR to config file {template}')

with open(template, 'w') as fw:
    yaml.dump(config, fw)

def md5(fname):
    # https://stackoverflow.com/questions/3431825/generating-an-md5-checksum-of-a-file
    hash_md5 = hashlib.md5()
    if not os.path.exists(fname):
        return None
    with open(fname, "rb") as f:
        for chunk in iter(lambda: f.read(4096), b""):
            hash_md5.update(chunk)
    return hash_md5.hexdigest()

D_FILE2MD5 = {
        'db.tgz': '727ed8addef563c14ca2c374f36c0c55',
        'combined.hmm': '9bb6d4e968aede068aa9fc414496473f',
        'Pfam-A-Archaea.hmm': '4ac89aba0a383b4e12957f260762820c',
        'Pfam-A-Bacteria.hmm': '13a47f8966737d456d27381ad0983cb0',
        'Pfam-A-Eukaryota.hmm': '3675151b16e675a4f979c76e7e6f5d18',
        'Pfam-A-Mixed.hmm': '450e93a2aa29463c8bee9ca5b1c55c79',
}

rule all:
    input: 'Done_all_setup'

rule download_hmm_viral:
    output: temp('combined.hmm.gz.split_{index}')
    shell:
        """
        for i in {{1..5}}; do
            Ret=$(wget --tries 2 --retry-connrefused --waitretry=60 --timeout=60 -q -O {output} https://zenodo.org/record/3823805/files/{output}?download=1 || echo "404 Not Found")
            if echo $Ret | grep -i -q "404 Not Found"; then
                if [ $i = 5 ]; then
                    exit 1
                fi
                echo "Server not responding for {output}; try $i/5.." | python {Scriptdir}/echo.py
                continue
            else
                break
            fi
        done
        """


rule download_hmm_pfam:
    output: temp('Pfam-A-{domain}.hmm')
    shell:
        """
        for i in {{1..5}}; do
            Ret=$(wget --tries 2 --retry-connrefused --waitretry=60 --timeout=60 -q -O {output}.gz https://zenodo.org/record/3823805/files/{output}.gz?download=1 || echo "404 Not Found")
            if echo $Ret | grep -i -q "404 Not Found"; then
                if [ $i = 5 ]; then
                    exit 1
                fi
                echo "Server not responding for {output}; try $i/5.." | python {Scriptdir}/echo.py
                continue
            else
                break
            fi
        done
        gunzip Pfam-A-{wildcards.domain}.hmm.gz
        """

rule download_other_db:
    output: temp('db.tgz')
    shell:
        """
        for i in {{1..5}}; do
            Ret=$(wget --tries 2 --retry-connrefused --waitretry=60 --timeout=60 -q -O db.tgz https://zenodo.org/record/3823805/files/db.tgz?download=1 || echo "404 Not Found")
            if echo $Ret | grep -i -q "404 Not Found"; then
                if [ $i = 5 ]; then
                    exit 1
                fi
                echo "*** Server not responding for {output}; try $i/5.." | python {Scriptdir}/echo.py
                continue
            else
                break
            fi
        done
        """

rule download_plamid_genomes:
    output: 'genomes-plasmid.fna.gz'
    shell:
        """
        for i in {{1..5}}; do
            Ret=$(wget --tries 2 --retry-connrefused --waitretry=60 --timeout=60 -q -O {output} https://zenodo.org/record/3823805/files/{output}?download=1 || echo "404 Not Found")
            if echo $Ret | grep -i -q "404 Not Found"; then
                if [ $i = 5 ]; then
                    exit 1
                fi
                echo "Server not responding for {output}; try $i/5.." | python {Scriptdir/echo.py
                continue
            else
                break
            fi
        done
        """

rule install_dependencies:
    output:
        temp(touch('Done-install-dependencies'))
    conda:
        '{}/vs2.yaml'.format(ENV_YAML_DIR)
    shell:
        """
        echo "Dependencies installed" | python {Scriptdir}/echo.py
        """

rule setup:
    input:
        'Done-install-dependencies',
        'db.tgz',
        expand('combined.hmm.gz.split_{index}', index=['00', '01', '02']),
        expand('Pfam-A-{domain}.hmm', 
                domain=['Archaea', 'Bacteria', 'Eukaryota', 'Mixed'])
    output:
        touch('Done_all_setup')
    run:
        shell(
        """
        rm -rf group hmm rbs
        tar -xzf db.tgz
        mv db/* .
        rmdir db
        mv Pfam-A-*.hmm hmm/pfam
        cat combined.hmm.gz.split_* | gunzip -c > hmm/viral/combined.hmm
        """
        )
        if md5('db.tgz') != D_FILE2MD5['db.tgz']:
            logging.error('Invalid checksum in for db.tgz')
            sys.exit(1)
        if md5('hmm/viral/combined.hmm') != D_FILE2MD5['combined.hmm']:
            logging.error('Invalid checksum in for combined.hmm')
            sys.exit(1)
        for domain in ['Archaea', 'Bacteria', 'Eukaryota', 'Mixed']:
            f = 'hmm/pfam/Pfam-A-{}.hmm'.format(domain)
            bname = 'Pfam-A-{}.hmm'.format(domain)
            if md5(f) != D_FILE2MD5[bname]:
                logging.error('Invalid checksum in for {}'.format(bname))
                sys.exit(1)

onerror:
    # clean up 
    if not os.path.exists('Done-install-dependencies'):
        shutil.rmtree('conda-envs')
        os.makedirs('conda-envs')

    fs = glob.glob('combined.hmm.gz.split*')
    fs.extend(glob.glob('Pfam-A-*.hmm'))
    fs.extend(['db.tgz'])
    for f in fs:
        if os.path.exists(f):
            os.remove(f)
    for di in ['group', 'hmm', 'rbs', 'db']:
        if os.path.exists(di):
            shutil.rmtree(di)

onsuccess:
    shell("""
        echo "All setup finished" | python {Scriptdir}/echo.py
        """)
