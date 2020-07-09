import hashlib
import os
import shutil
import glob
from ruamel.yaml import YAML


ENV_YAML_DIR = '../envs'
Srcdir = os.path.dirname(os.path.dirname(workflow.snakefile))
Scriptdir='{}/scripts'.format(Srcdir)


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
        'db.tgz': '9703c2d4f17a9714b3304fabbdfae3b2',
}

rule all:
    input: 'Done_all_setup'


rule download_db:
    output: temp('db.tgz')
    shell:
        """
        rm -f db.tgz
        wget -nv -O db.tgz https://osf.io/v46sc/download 
        echo "Download from osf finished.." | python {Scriptdir}/echo.py
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
    output:
        touch('Done_all_setup')
    run:
        shell(
        """
        rm -rf db
        tar -xzf db.tgz
        mv db/group db/hmm db/rbs .
        rm -rf db
        echo "All setup finished.." | python {Scriptdir}/echo.py
        """
        )
        assert md5('db.tgz') == D_FILE2MD5['db.tgz'], \
                '*** Invalid checksum in for db.tgz'

onstart:
    if not os.path.exists('Done-install-dependencies'):
        shutil.rmtree('conda-envs', ignore_errors=True)
        os.makedirs('conda-envs')

    fs = glob.glob('combined.hmm.gz.split*')
    fs.extend(glob.glob('Pfam-A-*.hmm'))
    fs.extend(['db.tgz'])
    for f in fs:
        if os.path.exists(f):
            os.remove(f)
    for di in ['group', 'hmm', 'rbs']:
        if os.path.exists(di):
            shutil.rmtree(di)

onerror:
    dbdir=os.path.abspath(os.getcwd())
    print(
        ('*** Download database from server failed '
        '(due to server temporary not responding or internet issue); '
        'You can download on you own through this link:\n'
        'https://osf.io/v46sc/download\n'
        'then untar and copy directories in "db" to {}').format(dbdir)
    )
