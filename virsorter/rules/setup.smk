import sys
import os
import shutil
import hashlib
import logging
from ruamel.yaml import YAML
from virsorter.config import set_logger, init_config_template

set_logger()

ENV_YAML_DIR = '../envs'
#ZENODO_ID = 3823805 
ZENODO_ID = 4269607

# update DBDIR in template-config.yaml
user_config_dir = os.path.join(os.path.expanduser('~'), '.virsorter')

# following DONOT work in Snakefile!!
#src_config_dir = os.path.dirname(os.path.abspath(__file__)) 

# not in the same dir as setup.smk, need to go up 2 levels
src_config_dir = os.path.dirname(os.path.dirname(workflow.snakefile))
Scriptdir=os.path.join(src_config_dir, 'scripts')

db_dir = os.path.abspath(os.getcwd())

init_config_template(src_config_dir, user_config_dir, db_dir)

#print(srcdir('.'))
#print(workflow.basedir)
#print(workflow.snakefile)

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
        'Pfam-A-Viruses.hmm': '7c0901ac28b3ba86561d0da06e97ee7b',
        'Pfam-A-acc2desc.tsv': 'a67245759c1529aee6485825ae5c3912',
}

rule all:
    input: 'Done_all_setup'

rule download_hmm_viral:
    output: temp('combined.hmm.gz.split_{index}')
    shell:
        """
        success=0
        for i in {{1..5}}; do
            Ret=$(wget --tries 2 --retry-connrefused --waitretry=60 --timeout=60 -q -O {output} https://zenodo.org/record/{ZENODO_ID}/files/{output}?download=1 || echo "404 Not Found")
            if echo $Ret | grep -i -q "404 Not Found"; then
                echo "Server not responding for {output}; try $i/5.." | python {Scriptdir}/echo.py
                sleep 2
                continue
            else
                success=1
                break
            fi
        done
        if [ $success -ne 1 ]; then
            echo "zenodo server not reponsive" | python {Scriptdir}/echo.py
            exit 1
        fi
        """

rule download_hmm_pfam:
    output: temp('Pfam-A-{domain}.hmm')
    shell:
        """
        success=0
        for i in {{1..5}}; do
            Ret=$(wget --tries 2 --retry-connrefused --waitretry=60 --timeout=60 -q -O {output}.gz https://zenodo.org/record/{ZENODO_ID}/files/{output}.gz?download=1 || echo "404 Not Found")
            if echo $Ret | grep -i -q "404 Not Found"; then
                echo "Server not responding for {output}; try $i/5.." | python {Scriptdir}/echo.py
                sleep 2
                continue
            else
                success=1
                break
            fi
        done
        if [ $success -ne 1 ]; then
            echo "zenodo server not reponsive" | python {Scriptdir}/echo.py
            exit 1
        fi
        gunzip Pfam-A-{wildcards.domain}.hmm.gz
        """

rule download_other_db:
    output: temp('db.tgz')
    shell:
        """
        success=0
        for i in {{1..5}}; do
            Ret=$(wget --tries 2 --retry-connrefused --waitretry=60 --timeout=60 -q -O db.tgz https://zenodo.org/record/{ZENODO_ID}/files/db.tgz?download=1 || echo "404 Not Found")
            if echo $Ret | grep -i -q "404 Not Found"; then
                echo "Server not responding for {output}; try $i/5.." | python {Scriptdir}/echo.py
                sleep 2
                continue
            else
                success=1
                break
            fi
        done
        if [ $success -ne 1 ]; then
            echo "zenodo server not reponsive" | python {Scriptdir}/echo.py
            exit 1
        fi
        """

rule download_misc:
    output: 
        #plasmid='genomes-plasmid.fna.gz',
        pfammap='Pfam-A-acc2desc.tsv',
    shell:
        """
        Target={output.pfammap}
        success=0
        for i in {{1..5}}; do
            Ret=$(wget --tries 2 --retry-connrefused --waitretry=60 --timeout=60 -q -O $Target.gz https://zenodo.org/record/{ZENODO_ID}/files/$Target.gz?download=1 || echo "404 Not Found")
            if echo $Ret | grep -i -q "404 Not Found"; then
                echo "server not responding for $Target.gz; try $i/5.." | python {Scriptdir}/echo.py
                sleep 2
                continue
            else
                success=1
                break
            fi
        done
        if [ $success -ne 1 ]; then
            echo "zenodo server not reponsive for $Target.gz" | python {Scriptdir}/echo.py
            exit 1
        fi
        gunzip $Target.gz
        """

if not config['Skip_deps_install']:
    rule install_dependencies:
        output:
            temp(touch('Done-install-dependencies'))
        conda:
            '{}/vs2.yaml'.format(ENV_YAML_DIR)
        shell:
            """
            echo "Dependencies installed" | python {Scriptdir}/echo.py
            """
else:
    rule install_dependencies:
        output:
            temp(touch('Done-install-dependencies'))
        shell:
            """
            echo "Dependencies installation skipped; make sure dependencies are installed on your own as shown in development version installation" | python {Scriptdir}/echo.py
            """

rule setup:
    input:
        'Done-install-dependencies',
        'db.tgz',
        'Pfam-A-acc2desc.tsv',
        expand('combined.hmm.gz.split_{index}', index=['00', '01', '02']),
        expand('Pfam-A-{domain}.hmm', 
                domain=['Archaea', 'Bacteria', 'Eukaryota', 'Mixed', 'Viruses'])
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
        mv Pfam-A-acc2desc.tsv hmm/pfam
        cat combined.hmm.gz.split_* | gunzip -c > hmm/viral/combined.hmm
        """
        )
        if md5('db.tgz') != D_FILE2MD5['db.tgz']:
            logging.error('Invalid checksum in for db.tgz')
            sys.exit(1)
        if md5('hmm/viral/combined.hmm') != D_FILE2MD5['combined.hmm']:
            logging.error('Invalid checksum in for combined.hmm')
            sys.exit(1)
        if md5('hmm/pfam/Pfam-A-acc2desc.tsv') != D_FILE2MD5['Pfam-A-acc2desc.tsv']:
            logging.error('Invalid checksum in for Pfam-A-acc2desc.tsv')
            sys.exit(1)

        for domain in ['Archaea', 'Bacteria', 'Eukaryota', 'Mixed', 'Viruses']:
            f = 'hmm/pfam/Pfam-A-{}.hmm'.format(domain)
            bname = 'Pfam-A-{}.hmm'.format(domain)
            if md5(f) != D_FILE2MD5[bname]:
                logging.error('Invalid checksum in for {}'.format(bname))
                sys.exit(1)

onerror:
    # clean up 
    if not os.path.exists('Done-install-dependencies'):
        shutil.rmtree('conda_envs')
        os.makedirs('conda_envs')

    fs = glob.glob('combined.hmm.gz.split*')
    fs.extend(glob.glob('Pfam-A-*.hmm'))
    fs.extend(['db.tgz', 'Pfam-A-acc2desc.tsv',])
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
