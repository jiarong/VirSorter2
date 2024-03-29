
rule gff_feature:
    input: f'{Tmpdir}/all.pdg.gff'
    output: f'{Tmpdir}/all.pdg.gff.ftr'
    log: f'log/{Tmpdir}/step2-extract-feature/extract-feature-from-gff-common.log'
    conda: '{}/vs2.yaml'.format(Conda_yaml_dir)
    shell:
        """
        python {Scriptdir}/extract-feature-from-gff.py {Dbdir}/rbs/rbs-catetory.tsv {input} {output} &> {log} || {{ echo "See error details in {Wkdir}/{log}" | python {Scriptdir}/echo.py --level error; exit 1; }}
        """

rule gff_feature_by_group:
    input:
        ftr=f'{Tmpdir}/all.pdg.gff.ftr',
        gff=f'{Tmpdir}/{{group}}/all.pdg.gff'
    output: f'{Tmpdir}/{{group}}/all.pdg.gff.ftr'
    conda: '{}/vs2.yaml'.format(Conda_yaml_dir)
    shell:
        """
        Log={Wkdir}/log/{Tmpdir}/step2-extract-feature/extract-feature-from-gff-{wildcards.group}.log
        Rbs_pdg_db={Dbdir}/group/{wildcards.group}/rbs-prodigal-train.db
        if [ -s $Rbs_pdg_db ]; then
            python {Scriptdir}/extract-feature-from-gff.py {Dbdir}/rbs/rbs-catetory.tsv {input.gff} {output} &> $Log || {{ echo "See error details in $Log" | python {Scriptdir}/echo.py --level error; exit 1; }}
        else
            (cd {Tmpdir}/{wildcards.group} && ln -sf ../all.pdg.gff.ftr)
        fi
        """

localrules: split_faa
checkpoint split_faa:
    input: f'{Tmpdir}/all.pdg.faa'
    output: directory(f'{Tmpdir}/all.pdg.faa.splitdir')
    conda: '{}/vs2.yaml'.format(Conda_yaml_dir)
    shell:
        """
        Log={Wkdir}/log/{Tmpdir}/step1-pp/split-faa-common.log
        Total=$(grep -v '^>' {input} | wc -c)
        Bname=$(basename {input})

        rm -f {Tmpdir}/all.pdg.faa.splitdir/all.pdg.faa.*.split.*.splithmmtbl

        if [ {Provirus} != True ] && [ {Max_orf_per_seq} -ne -1 ]; then
            echo "provirus mode is off; MAX_ORF_PER_SEQ set to {Max_orf_per_seq}; subsampling orf when orf number in a contig exceeds {Max_orf_per_seq} to speed up the run" | python {Scriptdir}/echo.py
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
    conda: '{}/vs2.yaml'.format(Conda_yaml_dir)
    shell:
        """
        Domain={wildcards.domain}
        if [ $Domain = "Viruses" ]; then
            Hmmdb={Dbdir}/hmm/viral/combined.hmm
        else
            Domain2=$Domain
            if [ $Domain2 = "Pfamviruses" ]; then
                Domain2=Viruses
            fi
            Hmmdb={Dbdir}/hmm/pfam/Pfam-A-"$Domain2".hmm
        fi

        Bname=$(basename {input})
        To_scratch=false
        # move the heavy IO of hmmsearch in local scratch if possible
        if [ -d "{Local_scratch}" ]; then
            # not sure df or du are compatible in all linux; use "||To_scratch=false" 
            #   to prevent imcompatibility in some linux distro
            Tmp=$(mktemp -d {Local_scratch}/vs2-XXXXXXXXXXXX) && To_scratch=true || To_scratch=false
            Avail=$(df -P {Local_scratch} | awk 'END{{print $4}}') || To_scratch=false
            Fsize=$(du -k {input} | awk '{{print $1*5}}') || To_scratch=false
            if [ "$Avail" -gt "$Fsize" ] && [ "$To_scratch" = "true" ]; then
                cp {input} $Tmp/$Bname || To_scratch=false
            else
                To_scratch=false
            fi
        fi

        hmmsearch -h | grep '^# HMMER' > {output}.log
        if [ "$To_scratch" = false ]; then
            # local scratch not set or not enough space in local scratch
            hmmsearch -T {Hmmsearch_score_min} --tblout {output} --cpu {threads} --noali $Hmmdb {input} > /dev/null 2>> {output}.log || {{ echo "See error details in {Wkdir}/{output}.log" | python {Scriptdir}/echo.py --level error; exit 1; }}
        else
            # when To_scratch is true, Tmp and Bname should have been defined successfully
            hmmsearch -T {Hmmsearch_score_min} --tblout {output} --cpu {threads} --noali $Hmmdb $Tmp/$Bname > /dev/null 2>> {output}.log || {{ echo "See error details in {Wkdir}/{output}.log" | python {Scriptdir}/echo.py --level error; exit 1; }}

            rm -f $Tmp/$Bname && rmdir $Tmp
        fi
        rm -f {output}.log
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
        printf "%s\n" {input} | xargs cat > {output}
        """

localrules: split_faa_by_group
checkpoint split_faa_by_group:
    input: f'{Tmpdir}/{{group}}/all.pdg.faa'
    output: directory(f'{Tmpdir}/{{group}}/all.pdg.faa.splitdir')
    conda: '{}/vs2.yaml'.format(Conda_yaml_dir)
    shell:
        """
        # make sure grep command below does not fail if input is empty
        set +o pipefail 

        Log={Wkdir}/log/{Tmpdir}/step1-pp/split-faa-{wildcards.group}.log
        Group_specific_hmmdb={Dbdir}/group/{wildcards.group}/customized.hmm
        Rbs_pdg_db={Dbdir}/group/{wildcards.group}/rbs-prodigal-train.db
        Bname=$(basename {input})

        rm -f {Tmpdir}/{wildcards.group}/all.pdg.faa.splitdir/all.pdg.faa.*.split.*.splithmmtbl

        if [ {Provirus} != "True" ] && [ {Max_orf_per_seq} -ne -1 ]; then
            python {Scriptdir}/subsample-faa.py {Max_orf_per_seq} {input} > {input}.ss
        else
            (cd {Tmpdir}/{wildcards.group} && ln -sf $Bname $Bname.ss)
        fi

        if [ -s $Rbs_pdg_db ] || [ -s $Group_specific_hmmdb ]; then
            Total=$(grep -v '^>' {input}.ss | wc -c)
            if [ $Total -gt {Faa_bp_per_split} ]; then
                python {Scriptdir}/split-seqfile-even-bp-per-file.py {input}.ss {output} {Faa_bp_per_split}  &> $Log || {{ echo "See error details in $Log" | python {Scriptdir}/echo.py --level error; exit 1; }}
            else
                # it's just small dataset, no need to split
                #echo "Dataset is smaller than {Faa_bp_per_split}bp, no need to split" | python {Scriptdir}/echo.py
                mkdir -p {output}
                (cd {output} && ln -sf ../$Bname $Bname.0.split)
            fi
        else
            # there is no group specific rbs/hmmdb 
            #echo "{wildcards.group} do not use group specific rbs or hmm DB, so just use the common hmmsearch annotation; skipping the faa split and hmmsearch" | python {Scriptdir}/echo.py
            mkdir -p {output}
            (cd {output} && ln -sf ../$Bname $Bname.0.split)
        fi
        """

rule hmmsearch_by_group:
    input: f'{Tmpdir}/{{group}}/all.pdg.faa.splitdir/all.pdg.faa.{{i}}.split'
    output: temp(f'{Tmpdir}/{{group}}/all.pdg.faa.splitdir/all.pdg.faa.{{i}}.split.{{domain}}.splithmmtbl')
    threads: Hmmsearch_threads
    conda: '{}/vs2.yaml'.format(Conda_yaml_dir)
    shell:
        """
        Domain={wildcards.domain}
        Group_specific_hmmdb={Dbdir}/group/{wildcards.group}/customized.hmm
        Rbs_pdg_db={Dbdir}/group/{wildcards.group}/rbs-prodigal-train.db

        if [ $Domain = "Viruses" ]; then
            if [ -s $Group_specific_hmmdb ]; then
                Hmmdb=$Group_specific_hmmdb
            else
                Hmmdb={Dbdir}/hmm/viral/combined.hmm
            fi
        else
            Domain2=$Domain
            if [ $Domain2 = "Pfamviruses" ]; then
                Domain2=Viruses
            fi
            Hmmdb={Dbdir}/hmm/pfam/Pfam-A-"$Domain2".hmm
        fi

        if [ -s $Rbs_pdg_db ] || [ -s $Group_specific_hmmdb ]; then
            To_scratch=false
            Bname=$(basename {input})
            # move the heavy IO of hmmsearch in local scratch
            if [ -d "{Local_scratch}" ]; then
                # not sure df or du are compatible in all linux; use "||To_scratch=false" 
                #   to prevent imcompatibility in some linux distro
                Tmp=$(mktemp -d {Local_scratch}/vs2-XXXXXXXXXXXX) && To_scratch=true || To_scratch=false
                Avail=$(df -P {Local_scratch} | awk 'END{{print $4}}') || To_scratch=false
                Fsize=$(du -k {input} | awk '{{print $1*5}}') || To_scratch=false
                if [ "$Avail" -gt "$Fsize" ] && [ "$To_scratch" = "true" ]; then
                    cp {input} $Tmp/$Bname || To_scratch=false
                else
                    To_scratch=false
                fi
            fi

            hmmsearch -h | grep '^# HMMER' > {output}.log
            if [ "$To_scratch" = false ]; then
                # local scratch not set or not enough space in local scratch
                hmmsearch -T {Hmmsearch_score_min} --tblout {output} --cpu {threads} --noali $Hmmdb {input} > /dev/null 2>> {output}.log || {{ echo "See error details in {Wkdir}/{output}.log" | python {Scriptdir}/echo.py --level error; exit 1; }}
            else
                # when To_scratch is true, Tmp and Bname should have been defined successfully
                hmmsearch -T {Hmmsearch_score_min} --tblout {output} --cpu {threads} --noali $Hmmdb $Tmp/$Bname > /dev/null 2>> {output}.log || {{ echo "See error details in {Wkdir}/{output}.log" | python {Scriptdir}/echo.py --level error; exit 1; }}
                rm -f $Tmp/$Bname && rmdir $Tmp
            fi
        else
            touch {output}
        fi
        """

def merge_split_hmmtbl_by_group_input_agg(wildcards):
    # the key line to tell snakemake this depend on a checkpoint
    split_dir = checkpoints.split_faa_by_group.get(**wildcards).output[0]

    splits = glob_wildcards(
        os.path.join(split_dir, 'all.pdg.faa.{i}.split')).i
    _s = 'all.pdg.faa.{{i}}.split.{domain}.splithmmtbl'.format(
        domain=wildcards.domain)
    _s = os.path.join(split_dir, _s)
    fs = expand(_s, i=splits)
    return fs

localrules: merge_split_hmmtbl_by_group_tmp
rule merge_split_hmmtbl_by_group_tmp:
    input: merge_split_hmmtbl_by_group_input_agg
    output: temp(f'{Tmpdir}/{{group}}/all.pdg.{{domain}}.hmmtbl.tmp'),
    shell:
        """
        Group_specific_hmmdb={Dbdir}/group/{wildcards.group}/customized.hmm
        Rbs_pdg_db={Dbdir}/group/{wildcards.group}/rbs-prodigal-train.db
        if [ -s $Rbs_pdg_db ] || [ -s $Group_specific_hmmdb ]; then
            printf "%s\n" {input} | xargs cat > {output}
        else
            touch {output}
        fi
        """

localrules: merge_split_hmmtbl_by_group
rule merge_split_hmmtbl_by_group:
    input:
        f'{Tmpdir}/{{group}}/all.pdg.{{domain}}.hmmtbl.tmp',
        f'{Tmpdir}/all.pdg.{{domain}}.hmmtbl',
    output: f'{Tmpdir}/{{group}}/all.pdg.{{domain}}.hmmtbl'
    shell:
        """
        Group_specific_hmmdb={Dbdir}/group/{wildcards.group}/customized.hmm
        Rbs_pdg_db={Dbdir}/group/{wildcards.group}/rbs-prodigal-train.db
        if [ -s $Rbs_pdg_db ] || [ -s $Group_specific_hmmdb ]; then
            cp {input[0]} {output}
        else
            (cd {Tmpdir}/{wildcards.group} && ln -s ../all.pdg.{wildcards.domain}.hmmtbl)
        fi
        """
        
if Prep_for_dramv:
    rule hmm_sort_to_best_hit_taxon:
        input: 
            arc = f'{Tmpdir}/all.pdg.Archaea.hmmtbl',
            bac = f'{Tmpdir}/all.pdg.Bacteria.hmmtbl',
            euk = f'{Tmpdir}/all.pdg.Eukaryota.hmmtbl',
            mix = f'{Tmpdir}/all.pdg.Mixed.hmmtbl',
            vir = f'{Tmpdir}/all.pdg.Viruses.hmmtbl',
            pfamvir = f'{Tmpdir}/all.pdg.Pfamviruses.hmmtbl',
            faa = f'{Tmpdir}/all.pdg.faa',
        output: 
            tax = f'{Tmpdir}/all.pdg.hmm.tax',
            taxpfam = f'{Tmpdir}/all.pdg.hmm.taxpfam',
            taxwhm = f'{Tmpdir}/all.pdg.hmm.taxwhm',
            ftr = f'{Tmpdir}/all.pdg.hmm.ftr'
        log: f'log/{Tmpdir}/step2-extract-feature/extract-feature-from-hmmout-common.log'
        conda: '{}/vs2.yaml'.format(Conda_yaml_dir)
        shell:
            """
            python {Scriptdir}/extract-feature-from-hmmout.py {Hmmsearch_score_min} "{input.arc},{input.bac},{input.euk},{input.mix},{input.vir}" "arc,bac,euk,mixed,vir" > {output.tax} 2> {log} || {{ echo "See error details in {Wkdir}/{log}" | python {Scriptdir}/echo.py --level error; exit 1; }}
            python {Scriptdir}/add-unaligned-to-hmm-featrues.py {input.faa} {output.tax} > {output.ftr}

            # pfam only annotation
            python {Scriptdir}/extract-feature-from-hmmout.py {Hmmsearch_score_min} "{input.arc},{input.bac},{input.euk},{input.mix},{input.pfamvir}" "arc,bac,euk,mixed,vir" > {output.tax}pfam 2>> {log} || {{ echo "See error details in {Wkdir}/{log}" | python {Scriptdir}/echo.py --level error; exit 1; }}
            python {Scriptdir}/add-unaligned-to-hmm-featrues.py {input.faa} {output.tax}pfam > {output.ftr}pfam

            # add hallmark info to .tax file for making affi-contigs.tab file
            python {Scriptdir}/add-hallmark-to-taxfile.py {output.tax} {output.tax}whm

            # remove split faa files
            echo {Tmpdir}/all.pdg.faa.*.split | xargs rm -f
            """
else:
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
            taxwhm = f'{Tmpdir}/all.pdg.hmm.taxwhm',
            ftr = f'{Tmpdir}/all.pdg.hmm.ftr'
        log: f'log/{Tmpdir}/step2-extract-feature/extract-feature-from-hmmout-common.log'
        conda: '{}/vs2.yaml'.format(Conda_yaml_dir)
        shell:
            """
            python {Scriptdir}/extract-feature-from-hmmout.py {Hmmsearch_score_min} "{input.arc},{input.bac},{input.euk},{input.mix},{input.vir}" "arc,bac,euk,mixed,vir" > {output.tax} 2> {log} || {{ echo "See error details in {Wkdir}/{log}" | python {Scriptdir}/echo.py --level error; exit 1; }}
            python {Scriptdir}/add-unaligned-to-hmm-featrues.py {input.faa} {output.tax} > {output.ftr}

            # add hallmark info to .tax file for making affi-contigs.tab file
            python {Scriptdir}/add-hallmark-to-taxfile.py {output.tax} {output.tax}whm
            # remove split faa files
            echo {Tmpdir}/all.pdg.faa.*.split | xargs rm -f
            """

if Prep_for_dramv:
    rule hmm_sort_to_best_hit_taxon_by_group:
        input: 
            tax = '{Tmpdir}/all.pdg.hmm.tax',
            faa = '{Tmpdir}/{group}/all.pdg.faa',
            arc = '{Tmpdir}/{group}/all.pdg.Archaea.hmmtbl',
            bac = '{Tmpdir}/{group}/all.pdg.Bacteria.hmmtbl',
            euk = '{Tmpdir}/{group}/all.pdg.Eukaryota.hmmtbl',
            mix = '{Tmpdir}/{group}/all.pdg.Mixed.hmmtbl',
            vir = '{Tmpdir}/{group}/all.pdg.Viruses.hmmtbl',
            pfamvir = '{Tmpdir}/{group}/all.pdg.Pfamviruses.hmmtbl',
        output: 
            tax = '{Tmpdir}/{group}/all.pdg.hmm.tax',
            taxpfam = '{Tmpdir}/{group}/all.pdg.hmm.taxpfam',
            taxwhm = '{Tmpdir}/{group}/all.pdg.hmm.taxwhm',
        conda: '{}/vs2.yaml'.format(Conda_yaml_dir)
        shell:
            """
            Log={Wkdir}/log/{Tmpdir}/step2-extract-feature/extract-feature-from-hmmout-{wildcards.group}.log
            Hallmark_list_f={Dbdir}/group/{wildcards.group}/hallmark-gene.list
            Group_specific_hmmdb={Dbdir}/group/{wildcards.group}/customized.hmm
            Rbs_pdg_db={Dbdir}/group/{wildcards.group}/rbs-prodigal-train.db
            if [ -s $Rbs_pdg_db ] || [ -s $Group_specific_hmmdb ]; then
                python {Scriptdir}/extract-feature-from-hmmout.py {Hmmsearch_score_min} "{input.arc},{input.bac},{input.euk},{input.mix},{input.vir}" "arc,bac,euk,mixed,vir" > {output.tax} 2> $Log || {{ echo "See error details in $Log" | python {Scriptdir}/echo.py --level error; exit 1; }}
                # pfam only annotation
                python {Scriptdir}/extract-feature-from-hmmout.py {Hmmsearch_score_min} "{input.arc},{input.bac},{input.euk},{input.mix},{input.pfamvir}" "arc,bac,euk,mixed,vir" > {output.tax}pfam 2>> $Log || {{ echo "See error details in $Log" | python {Scriptdir}/echo.py --level error; exit 1; }}
            else
                (cd {Tmpdir}/{wildcards.group} && ln -sf ../all.pdg.hmm.tax)
                (cd {Tmpdir}/{wildcards.group} && ln -sf ../all.pdg.hmm.taxpfam)
            fi

            if [ -s $Hallmark_list_f ]; then
                # add hallmark info to .tax file for making affi-contigs.tab file
                python {Scriptdir}/add-hallmark-to-taxfile.py {output.tax} {output.tax}whm --hallmark $Hallmark_list_f 2>> $Log || {{ echo "See error details in $Log" | python {Scriptdir}/echo.py --level error; exit 1; }}
            else
                python {Scriptdir}/add-hallmark-to-taxfile.py {output.tax} {output.tax}whm 2>> $Log || {{ echo "See error details in $Log" | python {Scriptdir}/echo.py --level error; exit 1; }}
            fi

            # remove split faa files
            echo {Tmpdir}/{wildcards.group}/all.pdg.faa.*.split | xargs rm -f
            """
else:
    rule hmm_sort_to_best_hit_taxon_by_group:
        input: 
            tax = '{Tmpdir}/all.pdg.hmm.tax',
            faa = '{Tmpdir}/{group}/all.pdg.faa',
            arc = '{Tmpdir}/{group}/all.pdg.Archaea.hmmtbl',
            bac = '{Tmpdir}/{group}/all.pdg.Bacteria.hmmtbl',
            euk = '{Tmpdir}/{group}/all.pdg.Eukaryota.hmmtbl',
            mix = '{Tmpdir}/{group}/all.pdg.Mixed.hmmtbl',
            vir = '{Tmpdir}/{group}/all.pdg.Viruses.hmmtbl',
        output: 
            tax = '{Tmpdir}/{group}/all.pdg.hmm.tax',
            taxwhm = '{Tmpdir}/{group}/all.pdg.hmm.taxwhm',
        conda: '{}/vs2.yaml'.format(Conda_yaml_dir)
        shell:
            """
            Log={Wkdir}/log/{Tmpdir}/step2-extract-feature/extract-feature-from-hmmout-{wildcards.group}.log
            Hallmark_list_f={Dbdir}/group/{wildcards.group}/hallmark-gene.list
            Group_specific_hmmdb={Dbdir}/group/{wildcards.group}/customized.hmm
            Rbs_pdg_db={Dbdir}/group/{wildcards.group}/rbs-prodigal-train.db
            if [ -s $Rbs_pdg_db ] || [ -s $Group_specific_hmmdb ]; then
                python {Scriptdir}/extract-feature-from-hmmout.py {Hmmsearch_score_min} "{input.arc},{input.bac},{input.euk},{input.mix},{input.vir}" "arc,bac,euk,mixed,vir" > {output.tax} 2> $Log || {{ echo "See error details in $Log" | python {Scriptdir}/echo.py --level error; exit 1; }}
            else
                (cd {Tmpdir}/{wildcards.group} && ln -sf ../all.pdg.hmm.tax)
            fi

            if [ -s $Hallmark_list_f ]; then
                # add hallmark info to .tax file for making affi-contigs.tab file
                python {Scriptdir}/add-hallmark-to-taxfile.py {output.tax} {output.tax}whm --hallmark $Hallmark_list_f 2>> $Log || {{ echo "See error details in $Log" | python {Scriptdir}/echo.py --level error; exit 1; }}
            else
                python {Scriptdir}/add-hallmark-to-taxfile.py {output.tax} {output.tax}whm 2>> $Log || {{ echo "See error details in $Log" | python {Scriptdir}/echo.py --level error; exit 1; }}
            fi
            # remove split faa files
            echo {Tmpdir}/{wildcards.group}/all.pdg.faa.*.split | xargs rm -f
            """

localrules: hmm_features_by_group
rule hmm_features_by_group:
    input:
        ftr = f'{Tmpdir}/all.pdg.hmm.ftr',
        tax = f'{Tmpdir}/{{group}}/all.pdg.hmm.tax',
        faa = f'{Tmpdir}/{{group}}/all.pdg.faa'
    output: f'{Tmpdir}/{{group}}/all.pdg.hmm.ftr'
    conda: '{}/vs2.yaml'.format(Conda_yaml_dir)
    shell:
        """
        Log={Wkdir}/log/{Tmpdir}/step2-extract-feature/merge-feature-{wildcards.group}.log
        Hallmark_list_f={Dbdir}/group/{wildcards.group}/hallmark-gene.list
        Group_specific_hmmdb={Dbdir}/group/{wildcards.group}/customized.hmm
        Rbs_pdg_db={Dbdir}/group/{wildcards.group}/rbs-prodigal-train.db

        if [ -s $Hallmark_list_f ]; then
            python {Scriptdir}/add-unaligned-to-hmm-featrues.py {input.faa} {input.tax} --hallmark $Hallmark_list_f > {output} 2> $Log || {{ echo "See error details in $Log" | python {Scriptdir}/echo.py --level error; exit 1; }}
            if [ {Prep_for_dramv} = True ]; then
                python {Scriptdir}/add-unaligned-to-hmm-featrues.py {input.faa} {input.tax}pfam --hallmark $Hallmark_list_f > {output}pfam 2> $Log || {{ echo "See error details in $Log" | python {Scriptdir}/echo.py --level error; exit 1; }}
            fi

        elif [ -s $Rbs_pdg_db ] || [ -s $Group_specific_hmmdb ]; then
            python {Scriptdir}/add-unaligned-to-hmm-featrues.py {input.faa} {input.tax} > {output} 2> $Log || {{ echo "See error details in $Log" | python {Scriptdir}/echo.py --level error; exit 1; }} 
            python {Scriptdir}/add-unaligned-to-hmm-featrues.py {input.faa} {input.tax}pfam > {output}pfam 2>> $Log || {{ echo "See error details in $Log" | python {Scriptdir}/echo.py --level error; exit 1; }} 
        else
            (cd {Tmpdir}/{wildcards.group} && ln -fs ../all.pdg.hmm.ftr)
            if [ {Prep_for_dramv} = True ]; then
                (cd {Tmpdir}/{wildcards.group} && ln -fs ../all.pdg.hmm.ftrpfam)
            fi
        fi
        """

localrules: merge_hmm_gff_features_by_group
rule merge_hmm_gff_features_by_group:
    input:
        gff_ftr = f'{Tmpdir}/{{group}}/all.pdg.gff.ftr',
        hmm_ftr = f'{Tmpdir}/{{group}}/all.pdg.hmm.ftr'
    output: 
        merged_ftr = f'{Tmpdir}/{{group}}/all.pdg.ftr'
    conda: '{}/vs2.yaml'.format(Conda_yaml_dir)
    shell:
        """
        Log={Wkdir}/log/{Tmpdir}/step2-extract-feature/merge-feature-{wildcards.group}.log
        python {Scriptdir}/merge-hmm-gff-features.py {input.gff_ftr} {input.hmm_ftr} > {output.merged_ftr} 2>> $Log || {{ echo "See error details in $Log" | python {Scriptdir}/echo.py --level error; exit 1; }}
        """
