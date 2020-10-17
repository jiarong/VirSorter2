import sys
import os
from ruamel.yaml import YAML
from snakemake.utils import min_version

### set minimum snakemake version ###
min_version('5.8.1')


Viral_seqfiles=config['Viral_seqfile']
Min_length=config['Min_length'] 
Max_orf_per_seq=config['Max_orf_per_seq'] 
Viral_genome_as_bin=config['Viral_genome_as_bin'] 
Fragments_per_genome = config['Fragments_per_genome']



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

Hmmsearch_threads = config_default['HMMSEARCH_THREADS']
General_threads = config_default['GENERAL_THREADS']
Hmmsearch_score_min = config_default['HMMSEARCH_SCORE_MIN']
Proba_cutoff = config_default['PROBA_CUTOFF']

Contig_bp_per_split = config_default['CONTIG_BP_PER_SPLIT']
Faa_bp_per_split = config_default['FAA_BP_PER_SPLIT']
Gff_seqnum_per_split = config_default['GFF_SEQNUM_PER_SPLIT']
Dbdir = config_default['DBDIR']

Rbs=config['Rbs']
if Rbs == 'NA':
    Rbs = 'None'
else:
    assert os.path.exists(Rbs)

Hmmdb=config['Hmm']
if Hmmdb == 'NA':
    Hmmdb = '{}/hmm/viral/combined.hmm'.format(Dbdir)
    mes = f'*** --hmm not set, using {Hmmdb}'
    print(mes)
else:
    assert os.path.exists(Hmmdb)

Hallmark=config['Hallmark']
if Hallmark == 'NA':
    Hallmark = '{}/group/dsDNAphage/hallmark-gene.list'.format(Dbdir)
    mes = f'*** --hallmark not set, using {Hallmark}'
    print(mes)
else:
    assert os.path.exists(Hmmdb)
#Tmpdir=config['TMPDIR']
#Local_scratch=config['LOCAL_SCRATCH']
Tmpdir='iter-0'
Local_scratch='/tmp'

Srcdir = src_config_dir
Scriptdir='{}/scripts'.format(Srcdir)
Wkdir=os.path.abspath(os.getcwd())
Conda_yaml_dir = '../envs' # relative to rules/*.smk files

wildcard_constraints:
    domain='[A-Za-z]+',

rule all:
    input: 'all.pdg.ftr'

localrules: prep_fragments_from_genome
rule prep_fragments_from_genome:
    #input: Viral_seqfiles
    output: f'{Tmpdir}/fragments.fasta'
    conda: f'{Conda_yaml_dir}/vs2.yaml'
    shell:
        """
        mkdir -p log
        if [ ! {Viral_genome_as_bin} = "True" ]; then
            ### each contig is a genome
            python {Scriptdir}/get-list-from-seqfile.py {Viral_seqfiles} | awk 'BEGIN{{OFS="\t"}} {{print $1,1}}' > {Tmpdir}/accession-with-cnt.list
            python {Scriptdir}/prepdb-train-sample-random-fragments-per-contig-for-env.py {Fragments_per_genome} {Tmpdir}/accession-with-cnt.list {Viral_seqfiles} > {Tmpdir}/fragments.fasta 2> {Tmpdir}/fragments.stats || {{ cat {Tmpdir}/fragments.stats; exit 1; }}
        else
            ### NCLDV has MAGs as bins with > 1 contigs
            ###   each bin is a .fna.gz file
            echo {Viral_seqfiles} | tr ' ' '\n' | awk 'BEGIN{{OFS="\t"}} {{print $1,1}}' > {Tmpdir}/accession-with-cnt.list
            python {Scriptdir}/prepdb-train-sample-random-fragments-per-genome-for-env.py {Fragments_per_genome} {Tmpdir}/accession-with-cnt.list {Viral_seqfiles} > {Tmpdir}/fragments.fasta 2> {Tmpdir}/fragments.stats || {{ cat {Tmpdir}/fragments.stats; exit 1; }}
        fi
        """

localrules: split_contig_file
checkpoint split_contig_file:
    input: f'{Tmpdir}/fragments.fasta'
    output: directory(f'{Tmpdir}/fragments.fasta.splitdir')
    conda: '{}/vs2.yaml'.format(Conda_yaml_dir)
    shell:
        """
        Log=log/split-contig-file.log
        Total=$(grep -v '^>' {input} | wc -c)
        Bname=$(basename {input})
        if [ $Total -gt {Contig_bp_per_split} ]; then
            python {Scriptdir}/split-seqfile-even-bp-per-file.py {input} {output} {Contig_bp_per_split} &> $Log || {{ echo "See error details in $Log" | python {Scriptdir}/echo.py --level error; exit 1; }}
        else
            mkdir -p {output}
            (cd {output} && ln -sf ../$Bname $Bname.0.split)
        fi
        """

rule gene_call: 
    input: f'{Tmpdir}/fragments.fasta.splitdir/fragments.fasta.{{i}}.split'
    output: 
        gff=temp(f'{Tmpdir}/fragments.fasta.splitdir/fragments.fasta.{{i}}.split.pdg.splitgff'),
        faa=temp(f'{Tmpdir}/fragments.fasta.splitdir/fragments.fasta.{{i}}.split.pdg.splitfaa'),
    log: f'{Tmpdir}/fragments.fasta.splitdir/fragments.fasta.{{i}}.split.pdg.log'
    conda: '{}/vs2.yaml'.format(Conda_yaml_dir)
    shell:
        """
        Rbs_pdg_db={Rbs}
        if [ -s $Rbs_pdg_db ]; then
            prodigal -t $Rbs_pdg_db -i {input} -a {output.faa} -o {output.gff} -f gff &> {log} || {{ echo "See error details in {log}" | python {Scriptdir}/echo.py --level error; exit 1; }}
        else
            prodigal -p meta -i {input} -a {output.faa} -o {output.gff} -f gff  &> {log} || {{ echo "See error details in {log}" | python {Scriptdir}/echo.py --level error; exit 1; }}
        fi
        """

def merge_split_faa_gff_input_agg(wildcards):
    # the key line to tell snakemake this depend on a checkpoint
    contig_split_dir = \
        checkpoints.split_contig_file.get(**wildcards).output[0]

    #fs = glob.glob('{}/circular.ext.fna.*.split'.format(cp_output))
    _s = 'fragments.fasta.{i}.split'
    splits = glob_wildcards(os.path.join(contig_split_dir, _s)).i

    _s = 'fragments.fasta.{i}.split.pdg.splitgff'
    _s = os.path.join(contig_split_dir, _s)
    gff = expand(_s, i=splits)

    _s = 'fragments.fasta.{i}.split.pdg.splitfaa'
    _s = os.path.join(contig_split_dir, _s)
    faa = expand(_s, i=splits)

    return {'gff': gff, 'faa': faa}

localrules: merge_split_faa_gff
rule merge_split_faa_gff:
    input: unpack(merge_split_faa_gff_input_agg)
    output:
        gff=f'{Tmpdir}/all.pdg.gff',
        faa=f'{Tmpdir}/all.pdg.faa',
    conda: '{}/vs2.yaml'.format(Conda_yaml_dir)
    shell:
        """
        printf "%s\0" {input.gff} | xargs -0 cat > {output.gff}
        printf "%s\0" {input.faa} | xargs -0 cat > {output.faa}
        """

rule gff_feature:
    input: f'{Tmpdir}/all.pdg.gff'
    output: f'{Tmpdir}/all.pdg.gff.ftr'
    log: 'log/extract-feature-from-gff.log'
    conda: '{}/vs2.yaml'.format(Conda_yaml_dir)
    shell:
        """
        python {Scriptdir}/extract-feature-from-gff.py {Dbdir}/rbs/rbs-catetory.tsv {input} {output} &> {log} || {{ echo "See error details in {log}" | python {Scriptdir}/echo.py --level error; exit 1; }}
        """

localrules: split_faa
checkpoint split_faa:
    input: f'{Tmpdir}/all.pdg.faa'
    output: directory(f'{Tmpdir}/all.pdg.faa.splitdir')
    conda: '{}/vs2.yaml'.format(Conda_yaml_dir)
    shell:
        """
        Log={Tmpdir}/log/split-faa.log
        Total=$(grep -v '^>' {input} | wc -c)
        Bname=$(basename {input})
        if [ {Max_orf_per_seq} -ne -1 ]; then
            echo "MAX_ORF_PER_SEQ set to {Max_orf_per_seq}; subsampling orf when orf number in a contig exceeds {Max_orf_per_seq} to speed up the run" | python {Scriptdir}/echo.py
            python {Scriptdir}/subsample-faa.py {Max_orf_per_seq} {input} > {Tmpdir}/$Bname.ss
        else
            (cd {Tmpdir} && ln -sf $Bname $Bname.ss)
        fi
        if [ $Total -gt {Faa_bp_per_split} ]; then
            python {Scriptdir}/split-seqfile-even-bp-per-file.py {Tmpdir}/all.pdg.faa.ss {output} {Faa_bp_per_split}  &> $Log || {{ echo "See error details in $Log" | python {Scriptdir}/echo.py --level error; exit 1; }}
        else
            mkdir -p {output}
            (cd {output} && ln -sf ../$Bname $Bname.0.split)
        fi
        """

rule hmmsearch:
    input: f'{Tmpdir}/all.pdg.faa.splitdir/all.pdg.faa.{{i}}.split'
    output: temp(f'{Tmpdir}/all.pdg.faa.splitdir/all.pdg.faa.{{i}}.split.{{domain}}.splithmmtbl')
    threads: Hmmsearch_threads
    log: f'{Tmpdir}/all.pdg.faa.splitdir/all.pdg.faa.{{i}}.split.{{domain}}.hmm.log'
    conda: '{}/vs2.yaml'.format(Conda_yaml_dir)
    shell:
        """
        Domain={wildcards.domain}
        if [ $Domain = "Viruses" ]; then
            Hmmdb={Hmmdb}
        else
            Hmmdb={Dbdir}/hmm/pfam/Pfam-A-{wildcards.domain}.hmm
        fi
        To_scratch=false
        # move the heavy IO of hmmsearch in local scratch
        if [ -d "{Local_scratch}" ]; then
            Tmp=$(mktemp -d {Local_scratch}/vs2-XXXXXXXXXXXX) && To_scratch=true
            Avail=$(df -P {Local_scratch} | awk 'END{{print $4}}')
            Fsize=$(du -k {input} | awk '{{print $1*5}}')
            if [ "$Avail" -gt "$Fsize" ] && [ "$To_scratch" = "true" ]; then
                Bname=$(basename {input})
                cp {input} $Tmp/$Bname || To_scratch=false
            fi
        fi
        if [ "$To_scratch" = false ]; then
            # local scratch not set or not enough space in local scratch
            hmmsearch -T {Hmmsearch_score_min} --tblout {output} --cpu {threads} --noali -o /dev/null $Hmmdb {input} &> {log} || {{ echo "See error details in {log}" | python {Scriptdir}/echo.py --level error; exit 1; }}
        else
            hmmsearch -T {Hmmsearch_score_min} --tblout {output} --cpu {threads} --noali -o /dev/null $Hmmdb $Tmp/$Bname &> {log} || {{ echo "See error details in {log}" | python {Scriptdir}/echo.py --level error; exit 1; }}
            rm -f $Tmp/$Bname && rmdir $Tmp
        fi
        """

def merge_split_hmmtbl_input_agg(wildcards):
    # the key line to tell snakemake this depend on a checkpoint
    split_dir = checkpoints.split_faa.get(**wildcards).output[0]

    splits = glob_wildcards(
        os.path.join(split_dir, 'all.pdg.faa.{i}.split')).i
    _s = 'all.pdg.faa.{{i}}.split.{domain}.splithmmtbl'.format(
        domain=wildcards.domain)
    _s = os.path.join(split_dir, _s)
    fs = expand(_s, i=splits)
    return fs

localrules: merge_split_hmmtbl
rule merge_split_hmmtbl:
    input: merge_split_hmmtbl_input_agg
    output: f'{Tmpdir}/all.pdg.{{domain}}.hmmtbl',
    shell:
        """
        printf "%s\0" {input} | xargs -0 cat > {output}
        """

rule hmm_sort_to_best_hit_taxon:
    input: 
        arc = f'{Tmpdir}/all.pdg.Archaea.hmmtbl',
        bac = f'{Tmpdir}/all.pdg.Bacteria.hmmtbl',
        euk = f'{Tmpdir}/all.pdg.Eukaryota.hmmtbl',
        mix = f'{Tmpdir}/all.pdg.Mixed.hmmtbl',
        vir = f'{Tmpdir}/all.pdg.Viruses.hmmtbl',
        faa = f'{Tmpdir}/all.pdg.faa',
    output: 
        tax = f'{Tmpdir}/all.pdg.hmm.tax',
        ftr = f'{Tmpdir}/all.pdg.hmm.ftr'
    log: 'log/extract-feature-from-hmmout.log'
    conda: '{}/vs2.yaml'.format(Conda_yaml_dir)
    shell:
        """
        python {Scriptdir}/extract-feature-from-hmmout.py {Hmmsearch_score_min} "{input.arc},{input.bac},{input.euk},{input.mix},{input.vir}" "arc,bac,euk,mixed,vir" > {output.tax} 2> {log} || {{ echo "See error details in {log}" | python {Scriptdir}/echo.py --level error; exit 1; }}
        Hallmark_list_f={Hallmark}
        if [ -s $Hallmark_list_f ]; then
            python {Scriptdir}/add-unaligned-to-hmm-featrues.py {input.faa} {output.tax} --hallmark $Hallmark_list_f > {output.ftr} 2>> {log} || {{ echo "See error details in {log}" | python {Scriptdir}/echo.py --level error; exit 1; }}
        else
            python {Scriptdir}/add-unaligned-to-hmm-featrues.py {input.faa} {output.tax} > {output.ftr} 2>> {log} || {{ echo "See error details in {log}" | python {Scriptdir}/echo.py --level error; exit 1; }}

        fi
        """

localrules: merge_hmm_gff_features
rule merge_hmm_gff_features:
    input:
        gff_ftr = f'{Tmpdir}/all.pdg.gff.ftr',
        hmm_ftr = f'{Tmpdir}/all.pdg.hmm.ftr'
    output: 
        merged_ftr = f'{Tmpdir}/all.pdg.ftr'
    conda: '{}/vs2.yaml'.format(Conda_yaml_dir)
    shell:
        """
        Log=log/merge-feature.log
        python {Scriptdir}/merge-hmm-gff-features.py {input.gff_ftr} {input.hmm_ftr} > {output.merged_ftr} 2>> $Log || {{ echo "See error details in $Log" | python {Scriptdir}/echo.py --level error; exit 1; }}
        """

rule finalize:
    input: 
        f'{Tmpdir}/all.pdg.ftr'
    output:
        f'all.pdg.ftr'
    shell:
        """
        cp {input} {output}
        echo -e "
        ===> Traiing feature finished..
        The output feature file: {Wkdir}/{output}
        <===
        " | python {Scriptdir}/echo.py
        """

for ru in workflow.rules:
    if not 'mem' in ru.resources:
        ru.resources['mem'] = config_default['DEFAULT_MEM']
    if not 'time' in ru.resources:
        ru.resources['time'] = config_default['DEFAULT_TIME']

container: "docker://continuumio/miniconda3"

#onsuccess:
#    print('virsorter run finished')
onerror:
    print('Refer the path to the log file of rules for troubleshooting')
    print('Issues can be raised at: https://github.com/jiarong/VirSorter2/issues')





