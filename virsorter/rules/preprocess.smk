
# detect circular, predict genes on contigs, and filter on size (# of genes) and/or if circular
localrules: circular_linear_split
checkpoint circular_linear_split:
    input: f'{Seqfile}'
    output: f'{Tmpdir}/pp-seqname-length.tsv',
    conda: '{}/vs2.yaml'.format(Conda_yaml_dir)
    shell: 
        """
        # prep_logdir
        mkdir -p log/{Tmpdir}/step1-pp log/{Tmpdir}/step2-extract-feature log/{Tmpdir}/step3-classify

        Cnt=$(grep -c '^>' {input})
        if [ ${{Cnt}} = 0 ]; then
            echo "No sequnences found in contig file; exiting" \
              | python {Scriptdir}/echo.py --level error
            exit 1
        fi 

        python {Scriptdir}/circular-linear-split.py \
          {input} \
          {Tmpdir}/pp-circular.fna.preext\
          {Tmpdir}/pp-linear.fna \
          {output[0]} \
          "||rbs:common" \
          {Min_length}

        if [ ! -s {Tmpdir}/pp-circular.fna.preext ]; then
            echo "No circular seqs found in contig file" \
              | python {Scriptdir}/echo.py
            rm {Tmpdir}/pp-circular.fna.preext
        else
            python {Scriptdir}/circular-extend.py \
              {Tmpdir}/pp-circular.fna.preext {Tmpdir}/pp-circular.fna
        fi

        if [ ! -s {Tmpdir}/pp-linear.fna ]; then
            echo "No linear seqs found in contig file" \
              | python {Scriptdir}/echo.py
            rm {Tmpdir}/pp-linear.fna
        fi
        """

localrules: circular_linear_split_by_group
rule circular_linear_split_by_group:
    input: f'{Tmpdir}/pp-{{shape}}.fna',
    output: f'{Tmpdir}/{{group}}/pp-{{shape}}.fna'
    shell: 
        """
        Rbs_pdg_db={Dbdir}/group/{wildcards.group}/rbs-prodigal-train.db
        if [ -s $Rbs_pdg_db ]; then
            sed 's/rbs:common/rbs:{wildcards.group}/g' {input} > {output}
        else
            (cd {Tmpdir}/{wildcards.group} && ln -sf ../pp-{wildcards.shape}.fna .)
        fi
        """

localrules: split_contig_file
checkpoint split_contig_file:
    input: f'{Tmpdir}/pp-{{shape}}.fna'
    output: directory(f'{Tmpdir}/pp-{{shape}}.fna.splitdir')
    conda: '{}/vs2.yaml'.format(Conda_yaml_dir)
    shell:
        """
        Log={Wkdir}/log/{Tmpdir}/step1-pp/split-contig-file-{wildcards.shape}-common.log
        Total=$(grep -v '^>' {input} | wc -c)
        Bname=$(basename {input})
        # clean up output from following steps to avoid bug with checkpoint
        rm -f {Tmpdir}/pp-{wildcards.shape}.fna.splitdir/pp-{wildcards.shape}.fna.*.split.pdg.splitgff
        rm -f {Tmpdir}/pp-{wildcards.shape}.fna.splitdir/pp-{wildcards.shape}.fna.*.split.pdg.splitfaa
        if [ $Total -gt {Contig_bp_per_split} ]; then
            python {Scriptdir}/split-seqfile-even-bp-per-file.py {input} {output} {Contig_bp_per_split} &> $Log || {{ echo "See error details in $Log" | python {Scriptdir}/echo.py --level error; exit 1; }}
        else
            mkdir -p {output}
            (cd {output} && ln -sf ../$Bname $Bname.0.split)
        fi
        echo "Finish spliting {wildcards.shape} contig file with common rbs" | python {Scriptdir}/echo.py
        """

rule gene_call: 
    input: f'{Tmpdir}/pp-{{shape}}.fna.splitdir/pp-{{shape}}.fna.{{i}}.split'
    output: 
        gff=temp(f'{Tmpdir}/pp-{{shape}}.fna.splitdir/pp-{{shape}}.fna.{{i}}.split.pdg.splitgff'),
        faa=temp(f'{Tmpdir}/pp-{{shape}}.fna.splitdir/pp-{{shape}}.fna.{{i}}.split.pdg.splitfaa'),
    conda: '{}/vs2.yaml'.format(Conda_yaml_dir)
    shell:
        """
        Log='{Tmpdir}/pp-{wildcards.shape}.fna.splitdir/pp-{wildcards.shape}.fna.{wildcards.i}.split.pdg.log'
        prodigal -p meta -i {input} -a {output.faa} -o {output.gff} -f gff  &> $Log || {{ echo "See error details in $Log" | python {Scriptdir}/echo.py --level error; exit 1; }}
        rm -f $Log
        """

def merge_split_faa_gff_input_agg(wildcards):
    # the key line to tell snakemake this depend on a checkpoint
    contig_split_dir = \
        checkpoints.split_contig_file.get(**wildcards).output[0]

    #fs = glob.glob('{}/circular.ext.fna.*.split'.format(cp_output))
    _s = 'pp-{shape}.fna.{{i}}.split'.format(shape=wildcards.shape)
    splits = glob_wildcards(os.path.join(contig_split_dir, _s)).i

    _s = os.path.join(contig_split_dir, _s)
    contig = expand(_s, i=splits)

    _s = 'pp-{shape}.fna.{{i}}.split.pdg.splitgff'.format(shape=wildcards.shape)
    _s = os.path.join(contig_split_dir, _s)
    gff = expand(_s, i=splits)

    _s = 'pp-{shape}.fna.{{i}}.split.pdg.splitfaa'.format(shape=wildcards.shape)
    _s = os.path.join(contig_split_dir, _s)
    faa = expand(_s, i=splits)

    return {'gff': gff, 'faa': faa, 'contig': contig}

localrules: merge_split_faa_gff
rule merge_split_faa_gff:
    input: unpack(merge_split_faa_gff_input_agg)
    output:
        gff=f'{Tmpdir}/pp-{{shape}}.gff',
        faa=f'{Tmpdir}/pp-{{shape}}.faa',
    conda: '{}/vs2.yaml'.format(Conda_yaml_dir)
    shell:
        """
        printf "%s\n" {input.gff} | xargs cat > {output.gff}
        printf "%s\n" {input.faa} | xargs cat > {output.faa}
        printf "%s\n" {input.contig} | xargs rm -f
        """

localrules: split_contig_file_by_group
checkpoint split_contig_file_by_group:
    input: f'{Tmpdir}/{{group}}/pp-{{shape}}.fna'
    output: directory(f'{Tmpdir}/{{group}}/pp-{{shape}}.fna.splitdir')
    conda: '{}/vs2.yaml'.format(Conda_yaml_dir)
    shell:
        """
        Log={Wkdir}/log/{Tmpdir}/step1-pp/split-contig-file-{wildcards.shape}-{wildcards.group}.log
        Bname=$(basename {input})
        Rbs_pdg_db={Dbdir}/group/{wildcards.group}/rbs-prodigal-train.db

        rm -f {Tmpdir}/{wildcards.group}/pp-{wildcards.shape}.fna.splitdir/pp-{wildcards.shape}.fna.*.split.pdg.splitgff
        rm -f {Tmpdir}/{wildcards.group}/pp-{wildcards.shape}.fna.splitdir/pp-{wildcards.shape}.fna.*.split.pdg.splitfaa

        if [ -s $Rbs_pdg_db ]; then
            Total=$(grep -v '^>' {input} | wc -c)
            if [ $Total -gt {Contig_bp_per_split} ]; then
                python {Scriptdir}/split-seqfile-even-bp-per-file.py {input} {output} {Contig_bp_per_split} &> $Log || {{ echo "See error details in $Log" | python {Scriptdir}/echo.py --level error; exit 1; }}
            else
                mkdir -p {output}
                (cd {output} && ln -sf ../$Bname $Bname.0.split)
            fi
            echo "Finish spliting {wildcards.shape} contig file with {wildcards.group} rbs" | python {Scriptdir}/echo.py
        else
            mkdir -p {output}
        fi
        """

rule gene_call_by_group_tmp:
    input: f'{Tmpdir}/{{group}}/pp-{{shape}}.fna.splitdir/pp-{{shape}}.fna.{{i}}.split'
    output: 
        gff=temp(f'{Tmpdir}/{{group}}/pp-{{shape}}.fna.splitdir/pp-{{shape}}.fna.{{i}}.split.pdg.splitgff'),
        faa=temp(f'{Tmpdir}/{{group}}/pp-{{shape}}.fna.splitdir/pp-{{shape}}.fna.{{i}}.split.pdg.splitfaa'),
    conda: '{}/vs2.yaml'.format(Conda_yaml_dir)
    shell:
        """
        Log={Wkdir}/{Tmpdir}/{wildcards.group}/pp-{wildcards.shape}.fna.splitdir/pp-{wildcards.shape}.fna.{wildcards.i}.split.pdg.log
        Rbs_pdg_db={Dbdir}/group/{wildcards.group}/rbs-prodigal-train.db
        if [ -s $Rbs_pdg_db ]; then
            prodigal -t $Rbs_pdg_db -i {input} -a {output.faa} -o {output.gff} -f gff &> $Log || {{ echo "See error details in $Log" | python {Scriptdir}/echo.py --level error; exit 1; }}
        else
            touch {output.gff} {output.faa}
        fi
        rm -f $Log
        """

def merge_split_faa_gff_by_group_input_agg(wildcards):
    # the key line to tell snakemake this depend on a checkpoint
    contig_split_dir = \
        checkpoints.split_contig_file_by_group.get(**wildcards).output[0]

    #fs = glob.glob('{}/circular.ext.fna.*.split'.format(cp_output))
    _s = 'pp-{shape}.fna.{{i}}.split'.format(shape=wildcards.shape)
    splits = glob_wildcards(os.path.join(contig_split_dir, _s)).i

    _s = 'pp-{shape}.fna.{{i}}.split.pdg.splitgff'.format(
        group=wildcards.group, 
        shape=wildcards.shape,
    )
    _s = os.path.join(contig_split_dir, _s)
    contig = expand(_s, i=splits)

    _s = 'pp-{shape}.fna.{{i}}.split.pdg.splitgff'.format(
        group=wildcards.group, 
        shape=wildcards.shape,
    )
    _s = os.path.join(contig_split_dir, _s)
    gff = expand(_s, i=splits)

    _s = 'pp-{shape}.fna.{{i}}.split.pdg.splitfaa'.format(
        group=wildcards.group, 
        shape=wildcards.shape,
    )
    _s = os.path.join(contig_split_dir, _s)
    faa = expand(_s, i=splits)

    return {'gff': gff, 'faa': faa, 'contig': contig}

localrules: merge_split_faa_gff_by_group
rule merge_split_faa_gff_by_group:
    input: 
        unpack(merge_split_faa_gff_by_group_input_agg),
        common_gff=f'{Tmpdir}/pp-{{shape}}.gff',
        common_faa=f'{Tmpdir}/pp-{{shape}}.faa',
    output:
        gff=f'{Tmpdir}/{{group}}/pp-{{shape}}.gff',
        faa=f'{Tmpdir}/{{group}}/pp-{{shape}}.faa',
    shell:
        """
        Rbs_pdg_db={Dbdir}/group/{wildcards.group}/rbs-prodigal-train.db
        if [ -s $Rbs_pdg_db ]; then
            printf "%s\n" {input.gff} | xargs cat > {output.gff}
            printf "%s\n" {input.faa} | xargs cat > {output.faa}
        else
            (cd {Tmpdir}/{wildcards.group} && ln -s ../pp-{wildcards.shape}.gff && ln -s ../pp-{wildcards.shape}.faa)
        fi
        printf "%s\n" {input.contig} | xargs rm -f
        """

localrules: remove_partial_gene
rule remove_partial_gene:
    input:
        gff=f'{Tmpdir}/pp-{{shape}}.gff',
        faa=f'{Tmpdir}/pp-{{shape}}.faa'
    output: 
        gff=f'{Tmpdir}/pp-{{shape}}-flt.gff',
        faa=f'{Tmpdir}/pp-{{shape}}-flt.faa'
    conda: '{}/vs2.yaml'.format(Conda_yaml_dir)
    shell:
        """
        Log={Wkdir}/log/{Tmpdir}/step1-pp/{wildcards.shape}-remove-partial-gene-common.log
        if [ {wildcards.shape} = "circular" ]; then
            python {Scriptdir}/circular-remove-partial-gene.py {input.gff} {output.gff} &> $Log || {{ echo "See error detail in $Log" | python {Scriptdir}/echo.py --level error; exit 1; }}
            python {Scriptdir}/filter-seqs-by-gff.py {output.gff} {input.faa} {output.faa}
        else
            (cd {Tmpdir} && ln -s pp-{wildcards.shape}.gff pp-{wildcards.shape}-flt.gff)
            (cd {Tmpdir} && ln -s pp-{wildcards.shape}.faa pp-{wildcards.shape}-flt.faa)
        fi
        """

localrules: remove_partial_gene_by_group
rule remove_partial_gene_by_group:
    input:
        gff=f'{Tmpdir}/pp-{{shape}}-flt.gff',
        faa=f'{Tmpdir}/pp-{{shape}}-flt.faa',
        group_gff=f'{Tmpdir}/{{group}}/pp-{{shape}}.gff',
        group_faa=f'{Tmpdir}/{{group}}/pp-{{shape}}.faa'
    output: 
        gff=f'{Tmpdir}/{{group}}/pp-{{shape}}-flt.gff',
        faa=f'{Tmpdir}/{{group}}/pp-{{shape}}-flt.faa'
    conda: '{}/vs2.yaml'.format(Conda_yaml_dir)
    shell:
        """
        Log={Wkdir}/log/{Tmpdir}/step1-pp/{wildcards.shape}-remove-partial-gene-{wildcards.group}.log
        Rbs_pdg_db={Dbdir}/group/{wildcards.group}/rbs-prodigal-train.db
        if [ -s $Rbs_pdg_db ]; then
            if [ {wildcards.shape} = "circular" ]; then
                python {Scriptdir}/circular-remove-partial-gene.py {input.group_gff} {output.gff} &> $Log || {{ echo "See error detail in $Log" | python {Scriptdir}/echo.py --level error; exit 1; }}
                python {Scriptdir}/filter-seqs-by-gff.py {output.gff} {input.group_faa} {output.faa}
            else
                (cd {Tmpdir}/{wildcards.group} && ln -s pp-{wildcards.shape}.gff pp-{wildcards.shape}-flt.gff)
                (cd {Tmpdir}/{wildcards.group} && ln -s pp-{wildcards.shape}.faa pp-{wildcards.shape}-flt.faa)
            fi
        else
            (cd {Tmpdir}/{wildcards.group} && ln -sf ../pp-{wildcards.shape}-flt.gff && ln -sf ../pp-{wildcards.shape}-flt.faa)
        fi
        """

def combine_linear_circular_input_agg(wildcards):
    # {Tmpdir}/seqname-length.tsv
    out = checkpoints.circular_linear_split.get(**wildcards).output[0]
    outdir = os.path.dirname(out)
    pat = os.path.join(outdir, 'pp-{shape}.fna')
    shapes = glob_wildcards(pat).shape
    faa_fs = expand('{outdir}/pp-{shape}-flt.faa', 
                        outdir=outdir, shape=shapes)
    gff_fs = expand('{outdir}/pp-{shape}-flt.gff', 
                        outdir=outdir, shape=shapes)
    # pp-circular.fna is extended by its seq (duplicate and concatenate)
    # so not the original sequence; use pp-circular.fna.preext
    # Maybe CHANGE: pp-circular.fna (extended) is needed for extract 
    # 	 viral segment adjusted for partial gene;
    # 	 for now just keep the original sequence since adjust partial gene
    # 	 on circular is not virsorter's job
    contig_fs = [f'{outdir}/pp-{shape}.fna.preext' if shape == 'circular' 
                    else f'{outdir}/pp-{shape}.fna'
                    for shape in shapes]
    return {'faa': faa_fs, 'gff': gff_fs, 'contig': contig_fs}

localrules: combine_linear_circular
rule combine_linear_circular:
    input:
        unpack(combine_linear_circular_input_agg)
    output:
        faa=f'{Tmpdir}/all.pdg.faa',
        gff=f'{Tmpdir}/all.pdg.gff',
        contig=f'{Tmpdir}/all.fna',
    shell:
        """
        # it's OK to cat gff directly
        cat {input.gff} > {output.gff}
        cat {input.faa} > {output.faa}
        cat {input.contig} > {output.contig}
        if [ ! -s {output.faa} ]; then
            printf "\n*** NOTE: No genes from the contigs are left in {output.faa} after preprocess; virsorter replies on features from these genes for prediction; check quality of your contigs (too short or strange sequence composition), or just discard this sample\n\n" | python {Scriptdir}/echo.py
            touch EARLY-EXIT
            exit 234
        fi
        """

def combine_linear_circular_by_group_input_agg(wildcards):
    # {Tmpdir}/seqname-length.tsv
    out = checkpoints.circular_linear_split.get(**wildcards).output[0]
    outdir = os.path.dirname(out)
    pat = os.path.join(outdir, 'pp-{shape}.fna')
    shapes = glob_wildcards(pat).shape
    faa_fs = expand('{outdir}/{group}/pp-{shape}-flt.faa', 
                        outdir=outdir, group=wildcards.group, shape=shapes)
    gff_fs = expand('{outdir}/{group}/pp-{shape}-flt.gff', 
                        outdir=outdir, group=wildcards.group, shape=shapes)

    return {'group_gff': gff_fs, 'group_faa': faa_fs}

localrules: combine_linear_circular_by_group
rule combine_linear_circular_by_group:
    input:
        unpack(combine_linear_circular_by_group_input_agg),
        faa=f'{Tmpdir}/all.pdg.faa',
        gff=f'{Tmpdir}/all.pdg.gff',
        #group_faa=combine_linear_circular_by_group_input_agg_faa,
        #group_faa=combine_linear_circular_by_group_input_agg_faa,
        #group_gff=combine_linear_circular_by_group_input_agg_gff,
    output: 
        faa=f'{Tmpdir}/{{group}}/all.pdg.faa',
        gff=f'{Tmpdir}/{{group}}/all.pdg.gff',
    shell:
        """
        Rbs_pdg_db={Dbdir}/group/{wildcards.group}/rbs-prodigal-train.db
        if [ -s $Rbs_pdg_db ]; then
            cat {input.group_faa} > {output.faa}
            cat {input.group_gff} > {output.gff}
            if [ ! -s {output.faa} ]; then
                printf "\n*** NOTE: No genes from the contigs are left in {output.faa} after preprocess; virsorter replies on features from these genes for prediction; check quality of your contigs (too short or strange sequence composition) or remove {wildcards.group} from Groups; or just discard this sample\n\n" | python {Scriptdir}/echo.py
                touch EARLY-EXIT
                exit 234
            fi
        else
            (cd {Tmpdir}/{wildcards.group} && ln -sf ../all.pdg.gff && ln -sf ../all.pdg.faa)
        fi
        """
