import hashlib
from ruamel.yaml import YAML

ENV_YAML_DIR = '../envs'


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
        wget -nv -O db.tgz https://osf.io/v46sc/download 
        """

rule install_dependencies:
    output:
        temp(touch('Done-install-dependencies'))
    conda:
        '{}/vs2.yaml'.format(ENV_YAML_DIR)
    shell:
        """
        echo "*** Dependencies installed"
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
        rm -rf db/group db/hmm db/rbs
        tar -xzf db.tgz
        mv db/group db/hmm db/rbs .
        rm -rf db
        """
        )
        assert md5('db.tgz') == D_FILE2MD5['db.tgz'], \
                '*** Invalid checksum in for db.tgz'
        print('*** All setup finished..')

onerror:
    dbdir=os.path.abspath(os.getcwd())
    print(
        ('*** Download database from server failed '
        '(due to server temporary not responding or internet issue); '
        'You can download on you own through this link:\n'
        'https://osf.io/v46sc/download\n'
        'then untar and copy directories in "db" to {}').format(dbdir)
    )
