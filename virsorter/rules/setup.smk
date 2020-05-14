import hashlib
from ruamel.yaml import YAML

ENV_YAML_DIR = '../envs'

# update DBDIR in template-config.yaml (not in the same dir as setup.smk)
#   need to go up 2 levels

# following DONOT work in Snakefile!!
#src_config_dir = os.path.dirname(os.path.abspath(__file__)) 

src_config_dir = os.path.dirname(srcdir('.'))
src_config_dir = os.path.dirname(src_config_dir)

#print(srcdir('.'))
#print(workflow.basedir)
#print(workflow.snakefile)

src_template = os.path.join(src_config_dir, 'template-config.yaml')
yaml = YAML()
with open(src_template) as fp:
    config = yaml.load(fp)
    config['DBDIR'] = os.path.abspath(os.getcwd())
with open(src_template, 'w') as fw:
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
    output: 'combined.hmm.gz.split_{index}'
    shell:
        """
          wget -nv -O {output} https://zenodo.org/record/3823805/files/{output}?download=1 || curl https://zenodo.org/record/3823805/files/{output}?download=1 > {output}
        """

rule download_hmm_pfam:
    output: 'Pfam-A-{domain}.hmm'
    shell:
        """
        wget -nv -O {output}.gz https://zenodo.org/record/3823805/files/{output}.gz?download=1 
        gunzip Pfam-A-{wildcards.domain}.hmm.gz
        """

rule download_other_db:
    output: 'db.tgz'
    shell:
        """
        wget -nv -O db.tgz https://zenodo.org/record/3823805/files/db.tgz?download=1 
        """

rule download_plamid_genomes:
    output: 'genomes-plasmid.fna.gz'
    shell:
        """
        wget -nv -O {output} https://zenodo.org/record/3823805/files/{output}?download=1 
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
        expand('combined.hmm.gz.split_{index}', index=['00', '01', '02']),
        expand('Pfam-A-{domain}.hmm', 
                domain=['Archaea', 'Bacteria', 'Eukaryota', 'Mixed'])
    output:
        touch('Done_all_setup')
    run:
        shell(
        """
        tar -xzf db.tgz
        mv db/* .
        mv Pfam-A-*.hmm hmm/pfam
        cat combined.hmm.gz.split_* | gunzip -c > hmm/viral/combined.hmm
        """
        )
        assert md5('db.tgz') == D_FILE2MD5['db.tgz'], \
                '*** Invalid checksum in for db.tgz'
        assert md5('hmm/viral/combined.hmm') == D_FILE2MD5['combined.hmm'], \
                '*** Invalid checksum in for combined.hmm'
        for domain in ['Archaea', 'Bacteria', 'Eukaryota', 'Mixed']:
            f = 'hmm/pfam/Pfam-A-{}.hmm'.format(domain)
            bname = 'Pfam-A-{}.hmm'.format(domain)
            assert md5(f) == D_FILE2MD5[bname], \
                '*** Invalid checksum in for {}'.format(bname)
