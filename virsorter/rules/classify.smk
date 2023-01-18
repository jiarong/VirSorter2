
rule classify_by_group:
    input: f'{Tmpdir}/{{group}}/all.pdg.ftr'
    output: f'{Tmpdir}/{{group}}/all.pdg.clf'
    conda: '{}/vs2.yaml'.format(Conda_yaml_dir)
    shell:
        """
        Log={Wkdir}/log/{Tmpdir}/step3-classify/all-score-{wildcards.group}.log
        python {Scriptdir}/classify.py {input} {Dbdir}/group/{wildcards.group}/model {wildcards.group} {output} 2> $Log || {{ echo "See error details in $Log" | python {Scriptdir}/echo.py --level error; exit 1; }}
        """

localrules: check_point_for_reclassify
rule check_point_for_reclassify:
    input: expand(f'{Tmpdir}/{{group}}/all.pdg.clf', group=Groups)
    output: touch(f'{Tmpdir}/reclassify.trigger')

localrules: merge_classification
rule merge_classification:
    input: 
        expand(f'{Tmpdir}/{{group}}/all.pdg.clf', group=Groups),
        f'{Tmpdir}/reclassify.trigger',
    output: 
        clf=f'{Tmpdir}/all-fullseq-proba.tsv',
    conda: '{}/vs2.yaml'.format(Conda_yaml_dir)
    params:
        clf_fs_str = ' '.join(
            [f'{Tmpdir}/{group}/all.pdg.clf' for group in Groups]
        ),
    shell:
        """
        Log={Wkdir}/log/{Tmpdir}/step3-classify/all-score-merge.log
        ### merge clf from all groups
        python {Scriptdir}/merge-clf-from-groups.py {output.clf} {params.clf_fs_str} 2> $Log || {{ echo "See error details in $Log" | python {Scriptdir}/echo.py --level error; exit 1; }}
        """

localrules: pick_viral_fullseq
rule pick_viral_fullseq:
    input:
        clf=f'{Tmpdir}/all-fullseq-proba.tsv',
    output: 
        contig=f'{Tmpdir}/viral-fullseq.fa',
        hmk_cnt=f'{Tmpdir}/all-hallmark-cnt.tsv',
        lt2gene=f'{Tmpdir}/viral-lt2gene-w-hallmark.fa',
    conda: '{}/vs2.yaml'.format(Conda_yaml_dir)
    params:
        group_str = ','.join(
            [group for group in Groups if os.path.exists(f'{Dbdir}/group/{group}/hallmark-gene.list')]
        ),
        hmk_fs_str = ','.join(
            [f'{Dbdir}/group/{group}/hallmark-gene.list' for group in Groups if os.path.exists(f'{Dbdir}/group/{group}/hallmark-gene.list')]
        ),
        tax_fs_str = ','.join(
            [f'{Tmpdir}/{group}/all.pdg.hmm.tax' for group in Groups if os.path.exists(f'{Dbdir}/group/{group}/hallmark-gene.list')]
        ),
        all_group_str = ','.join(Groups),
        all_tax_fs_str = ','.join(
            [f'{Tmpdir}/{group}/all.pdg.hmm.tax' for group in Groups]
        ),
        all_gff_fs_str = ','.join(
            [f'{Tmpdir}/{group}/all.pdg.gff' for group in Groups]
        ),
        all_ftr_fs_str = ','.join(
            [f'{Tmpdir}/{group}/all.pdg.ftr' for group in Groups]
        ),
    shell:
        """
        Log={Wkdir}/log/{Tmpdir}/step3-classify/pick-viral-fullseq.log
        python {Scriptdir}/pick-viral-contig-from-clf.py {Proba_cutoff} {input.clf} {Tmpdir}/all.fna > {output.contig}.tmp 2> $Log || {{ echo "See error details in $Log" | python {Scriptdir}/echo.py --level error; exit 1; }}

        python {Scriptdir}/add-extra-to-fullseq-fasta-header.py {output.contig}.tmp {params.all_ftr_fs_str} {params.all_group_str} > {output.contig} 2>> $Log || {{ echo "See error details in $Log" | python {Scriptdir}/echo.py --level error; exit 1; }}
        rm {output.contig}.tmp

        python {Scriptdir}/get-hallmark-cnt-for-each-seq.py {output.hmk_cnt} "{params.group_str}" "{params.hmk_fs_str}" "{params.tax_fs_str}" 2>> $Log || {{ echo "See error details in $Log" | python {Scriptdir}/echo.py --level error; exit 1; }}

        python {Scriptdir}/get-seq-w-lt2gene-w-hallmark.py {output.hmk_cnt} {input.clf} {Tmpdir}/all.fna > {output.lt2gene}.tmp 2>> $Log || {{ echo "See error details in $Log" | python {Scriptdir}/echo.py --level error; exit 1; }}

        python {Scriptdir}/add-extra-to-lt2gene-fasta-header.py {output.lt2gene}.tmp {params.all_gff_fs_str} {params.all_tax_fs_str} {params.all_group_str} > {output.lt2gene} 2>> $Log || {{ echo "See error details in $Log" | python {Scriptdir}/echo.py --level error; exit 1; }}
        rm {output.lt2gene}.tmp

        """

if Provirus:
    checkpoint split_gff_by_group:
        input: 
            gff=f'{Tmpdir}/{{group}}/all.pdg.gff',
            clf=f'{Tmpdir}/all-fullseq-proba.tsv',
        output: directory(f'{Tmpdir}/{{group}}/all.pdg.gff.splitdir')
        conda: '{}/vs2.yaml'.format(Conda_yaml_dir)
        shell:
            """
            rm -f {Tmpdir}/{wildcards.group}/all.pdg.gff.splitdir/all.pdg.gff.*.split.prv.bdy
            rm -f {Tmpdir}/{wildcards.group}/all.pdg.gff.splitdir/all.pdg.gff.*.split.prv.ftr
            rm -f {Tmpdir}/{wildcards.group}/all.pdg.gff.splitdir/all.pdg.gff.*.anno
            rm -f {Tmpdir}/{wildcards.group}/all.pdg.gff.splitdir/all.pdg.gff.*.affi.tab
            python {Scriptdir}/split-gff-even-seqnum-per-file.py {input.gff} {output} {Gff_seqnum_per_split}
            """

    rule provirus_call_by_group_by_split:
        input: 
            gff=f'{Tmpdir}/{{group}}/all.pdg.gff.splitdir/all.pdg.gff.{{idx}}.split',
            tax=f'{Tmpdir}/{{group}}/all.pdg.hmm.tax',
            clf=f'{Tmpdir}/all-fullseq-proba.tsv',
        output: 
            boundry=temp(f'{Tmpdir}/{{group}}/all.pdg.gff.splitdir/all.pdg.gff.{{idx}}.split.prv.bdy'),
            ftr=temp(f'{Tmpdir}/{{group}}/all.pdg.gff.splitdir/all.pdg.gff.{{idx}}.split.prv.ftr'),
        conda: '{}/vs2.yaml'.format(Conda_yaml_dir)
        shell:
            """
            Log={Wkdir}/{input.gff}.prv.log
            Hallmark_list_f={Dbdir}/group/{wildcards.group}/hallmark-gene.list
            if [ -s $Hallmark_list_f ]; then
                python {Scriptdir}/provirus.py {input.gff} {input.tax} {Dbdir}/rbs/rbs-catetory.tsv {Dbdir}/group/{wildcards.group}/model {output.boundry} {output.ftr} --fullseq-clf {input.clf} --group {wildcards.group} --proba {Proba_cutoff} --hallmark {Dbdir}/group/{wildcards.group}/hallmark-gene.list 2> $Log || {{ echo "See error details in $Log" | python {Scriptdir}/echo.py --level error; exit 1; }}
            else
                python {Scriptdir}/provirus.py {input.gff} {input.tax} {Dbdir}/rbs/rbs-catetory.tsv {Dbdir}/group/{wildcards.group}/model {output.boundry} {output.ftr} --fullseq-clf {input.clf} --group {wildcards.group} --proba {Proba_cutoff} 2> $Log || {{ echo "See error details in $Log" | python {Scriptdir}/echo.py --level error; exit 1; }}
            fi
            rm -f $Log
            """

    def merge_provirus_call_by_group_by_split_input_agg(wildcards):
        split_dir = checkpoints.split_gff_by_group.get(
                        **wildcards).output[0]
        idx_lis = glob_wildcards(
                '{}/all.pdg.gff.{{idx}}.split'.format(split_dir)).idx
        bdy_str = '{}/all.pdg.gff.{{idx}}.split.prv.bdy'.format(split_dir)
        ftr_str = '{}/all.pdg.gff.{{idx}}.split.prv.ftr'.format(split_dir)

        bdy_lis = expand(bdy_str, idx=idx_lis)
        ftr_lis = expand(ftr_str, idx=idx_lis)
        return {'bdy': bdy_lis, 'ftr': ftr_lis}

    localrules: merge_provirus_call_by_group_by_split
    rule merge_provirus_call_by_group_by_split:
        input: unpack(merge_provirus_call_by_group_by_split_input_agg)
        output: 
            bdy=f'{Tmpdir}/{{group}}/all.pdg.prv.bdy',
            ftr=f'{Tmpdir}/{{group}}/all.pdg.prv.ftr',
        conda: '{}/vs2.yaml'.format(Conda_yaml_dir)
        shell:
            """
            printf "%s\n" {input.bdy} | xargs cat | awk '!/^seqname\t/ || !f++' > {output.bdy}
            printf "%s\n" {input.ftr} | xargs cat | awk '!/^seqname\t/ || !f++' > {output.ftr}
            """

    localrules: merge_provirus_call_from_groups
    rule merge_provirus_call_from_groups:
        input: expand(f'{Tmpdir}/{{group}}/all.pdg.prv.bdy', group=Groups)
        output:
            partial=temp(f'{Tmpdir}/viral.partseq.tsv.tmp'),
            full=temp(f'{Tmpdir}/viral.fullseq.tsv.tmp'),
        conda: '{}/vs2.yaml'.format(Conda_yaml_dir)
        params:
            prv_fs_str = ','.join(
                [f'{Tmpdir}/{group}/all.pdg.prv.bdy' for group in Groups]
            ),
            groups_str = ','.join(Groups),
        shell:
            """
            Log={Wkdir}/log/{Tmpdir}/step3-classify/provirus-score-merge.log
            python {Scriptdir}/merge-provirus-from-groups.py {params.prv_fs_str} {params.groups_str} {output.partial} {output.full} 2> $Log || {{ echo "See error details in $Log" | python {Scriptdir}/echo.py --level error; exit 1; }}
            """

    localrules: extract_provirus_seqs
    rule extract_provirus_seqs:
        input:
            #contig={Seqfile},
            contig=f'{Tmpdir}/all.fna',
            full=f'{Tmpdir}/viral.fullseq.tsv.tmp',
            partial=f'{Tmpdir}/viral.partseq.tsv.tmp',
            hmk_cnt=f'{Tmpdir}/all-hallmark-cnt.tsv',
            lt2gene=f'{Tmpdir}/viral-lt2gene-w-hallmark.fa',
        output:
            fullseq=f'{Tmpdir}/viral-fullseq-trim.fa',
            partial=f'{Tmpdir}/viral-partseq.fa',
            combined=f'{Tmpdir}/viral-combined.fa',
            fulltab=f'{Tmpdir}/viral-fullseq.tsv',
            parttab=f'{Tmpdir}/viral-partseq.tsv',
        conda: '{}/vs2.yaml'.format(Conda_yaml_dir)
        shell:
            """
            Log={Wkdir}/log/{Tmpdir}/step3-classify/provirus-seq-extract.log
            python {Scriptdir}/extract-provirus-seqs.py {input.contig} {input.full} {input.partial} {output.fullseq} {output.partial} {output.fulltab} {output.parttab} 2> $Log || {{ echo "See error details in $Log" | python {Scriptdir}/echo.py --level error; exit 1; }}
            cat {output.fullseq} {output.partial} > {output.combined}
            python {Scriptdir}/add-suffix-seqname-keep-desc.py {input.lt2gene} "||lt2gene" >> {output.combined} 2>> $Log || {{ echo "See error details in $Log" | python {Scriptdir}/echo.py --level error; exit 1; }}
            """

    ##### add score for non-best-scoring groups
    rule classify_full_and_part_by_group:
        input: 
            gff=f'{Tmpdir}/{{group}}/all.pdg.gff',
            tax=f'{Tmpdir}/{{group}}/all.pdg.hmm.tax',
            seqfile=f'{Tmpdir}/viral-combined.fa',
        output: f'{Tmpdir}/{{group}}/viral.trim.clf'
        conda: '{}/vs2.yaml'.format(Conda_yaml_dir)
        shell:
            """
            Log={Wkdir}/log/{Tmpdir}/step3-classify/provirus-score-{wildcards.group}.log
            Hallmark_list_f={Dbdir}/group/{wildcards.group}/hallmark-gene.list
            if [ -s $Hallmark_list_f ]; then
                python {Scriptdir}/classify-trimed.py {input.gff} {input.tax} {Dbdir}/rbs/rbs-catetory.tsv {Dbdir}/group/{wildcards.group}/model {input.seqfile} {output} --group {wildcards.group} --hallmark {Dbdir}/group/{wildcards.group}/hallmark-gene.list 2> $Log || {{ echo "See error details in $Log" | python {Scriptdir}/echo.py --level error; exit 1; }}
            else
                python {Scriptdir}/classify-trimed.py {input.gff} {input.tax} {Dbdir}/rbs/rbs-catetory.tsv {Dbdir}/group/{wildcards.group}/model {input.seqfile} {output} --group {wildcards.group} 2> $Log || {{ echo "See error details in $Log" | python {Scriptdir}/echo.py --level error; exit 1; }}
            fi

            """

    localrules: merge_full_and_part_classification
    rule merge_full_and_part_classification:
        input: expand(f'{Tmpdir}/{{group}}/viral.trim.clf', group=Groups)
        output: 
            clf=f'{Tmpdir}/viral-combined-proba.tsv',
        conda: '{}/vs2.yaml'.format(Conda_yaml_dir)
        params:
            clf_fs_str = ' '.join(
                [f'{Tmpdir}/{group}/viral.trim.clf' for group in Groups]
            ),
        shell:
            """
            Log={Wkdir}/log/{Tmpdir}/step3-classify/provirus-score-merge.log
            ### merge clf from all groups
            python {Scriptdir}/merge-clf-trim-from-groups.py {output.clf} {params.clf_fs_str} 2> $Log || {{ echo "See error details in $Log" | python {Scriptdir}/echo.py --level error; exit 1; }}
            """

    rule make_annotation_table_by_group_by_split:
        input: 
            gff=f'{Tmpdir}/{{group}}/all.pdg.gff.splitdir/all.pdg.gff.{{idx}}.split',
            taxwhm=f'{Tmpdir}/{{group}}/all.pdg.hmm.taxwhm',
            taxpfam=f'{Tmpdir}/{{group}}/all.pdg.hmm.taxpfam',
            seqfile=f'{Tmpdir}/viral-combined.fa',
        output: 
            anno=temp(f'{Tmpdir}/{{group}}/all.pdg.gff.splitdir/all.pdg.gff.{{idx}}.anno'),
            affi=temp(f'{Tmpdir}/{{group}}/all.pdg.gff.splitdir/all.pdg.gff.{{idx}}.affi.tab'),
        conda: '{}/vs2.yaml'.format(Conda_yaml_dir)
        shell:
            """
            python {Scriptdir}/make-affi-contigs-tabfile-on-split-gff.py --pfamtax {input.taxpfam} {input.seqfile} {output.anno} {output.affi} {input.gff} {input.taxwhm} {wildcards.group}
            """

    def merge_annotation_table_by_group_from_split_input_agg(wildcards):
        split_dir = checkpoints.split_gff_by_group.get(
                        **wildcards).output[0]
        idx_lis = glob_wildcards(
                '{}/all.pdg.gff.{{idx}}.split'.format(split_dir)).idx
        anno_str = '{}/all.pdg.gff.{{idx}}.anno'.format(split_dir)
        affi_str = '{}/all.pdg.gff.{{idx}}.affi.tab'.format(split_dir)
        anno_lis = expand(anno_str, idx=idx_lis)
        affi_lis = expand(affi_str, idx=idx_lis)
        return {'anno': anno_lis, 'affi': affi_lis}

    localrules: merge_annotation_table_by_group_from_split
    rule merge_annotation_table_by_group_from_split:
        input: unpack(merge_annotation_table_by_group_from_split_input_agg)
        output: 
            anno=f'{Tmpdir}/{{group}}/all.pdg.anno',
            affi=f'{Tmpdir}/{{group}}/all.pdg.affi',
        conda: '{}/vs2.yaml'.format(Conda_yaml_dir)
        shell:
            """
            printf "%s\n" {input.anno} | xargs cat | awk '!/^seqname\t/ || !f++' > {output.anno}
            printf "%s\n" {input.affi} | xargs cat > {output.affi}
            """

    localrules: merge_annotation_table_from_groups
    rule merge_annotation_table_from_groups:
        input: 
            anno=expand(f'{Tmpdir}/{{group}}/all.pdg.anno', group=Groups),
            affi=expand(f'{Tmpdir}/{{group}}/all.pdg.affi', group=Groups),
        output:
            anno=f'{Tmpdir}/viral.anno',
            affi=f'{Tmpdir}/viral-affi-contigs-for-dramv.tab',
        conda: '{}/vs2.yaml'.format(Conda_yaml_dir)
        shell:
            """
            printf "%s\n" {input.anno} | xargs cat | awk '!/^seqname\t/ || !f++' > {output.anno}
            printf "%s\n" {input.affi} | xargs cat > {output.affi}
            """

    def finalize_input_agg(wildcards):
        if Prep_for_dramv:
            l = [
                    f'{Tmpdir}/viral-fullseq-trim.fa', 
                    f'{Tmpdir}/viral-partseq.fa', 
                    f'{Tmpdir}/viral-combined-proba.tsv', 
                    f'{Tmpdir}/viral-affi-contigs-for-dramv.tab',
            ]
        else:
            l = [
                    f'{Tmpdir}/viral-fullseq-trim.fa',
                    f'{Tmpdir}/viral-partseq.fa',
                    f'{Tmpdir}/viral-combined-proba.tsv',
            ]
        return l

    localrules: finalize
    rule finalize:
        input:
            finalize_input_agg
        output: 
            score=f'{Label}final-viral-score.tsv',
            fa=f'{Label}final-viral-combined.fa',
            boundary=f'{Label}final-viral-boundary.tsv',
        conda: '{}/vs2.yaml'.format(Conda_yaml_dir)
        shell:
            """
            echo {Tmpdir}/*/all.pdg.gff.splitdir/all.pdg.gff.*.split | xargs rm -f
            python {Scriptdir}/add-extra-to-table.py {Tmpdir}/viral-combined-proba.tsv {Tmpdir}/viral-combined.fa {Tmpdir}/viral-combined-proba-more-cols.tsv
            python {Scriptdir}/filter-score-table.py config.yaml {Tmpdir}/viral-combined-proba-more-cols.tsv {Tmpdir}/viral-combined.fa {output.score} {output.fa}.trim
            python {Scriptdir}/keep-original-seq.py {output.fa}.trim {Seqfile} > {output.fa}.original
            cp {Tmpdir}/viral-fullseq.tsv {output.boundary}.tmp
            tail -n +2 {Tmpdir}/viral-partseq.tsv >> {output.boundary}.tmp 
            python {Scriptdir}/finalize-boundary-file.py {output.boundary}.tmp {output.score} > {output.boundary}
	    rm -f {output.boundary}.tmp

            if [ {Keep_original_seq} = "True" ]; then
                cp {output.fa}.original {output.fa}
            else
                cp {output.fa}.trim {output.fa}
            fi

            if [ {Prep_for_dramv} = "True" ]; then
                mkdir -p {Label}for-dramv
                python {Scriptdir}/modify-seqname-for-dramv.py {output.fa}.original {output.score} -o {Label}for-dramv/final-viral-combined-for-dramv.fa
                cp {Tmpdir}/viral-affi-contigs-for-dramv.tab {Label}for-dramv
            fi
            rm -f {output.fa}.trim {output.fa}.original

            N_lt2gene=$(grep -c '^>.*||lt2gene$' {output.fa} || :)
            N_lytic=$(grep -c '^>.*||full$' {output.fa} || :)
            N_lysogenic=$(grep -c '^>.*||.*_partial$' {output.fa} || :)
            if [ {Prep_for_dramv} = True ]; then
                Dramv_notes="{Label}for-dramv                   ==> dir with input files for dramv"
                Dramv_notes2="For seqnames in files for dramv, 
                    | is replaced with _ to be compatible with DRAMv
                "
            else
                Dramv_notes=""
                Dramv_notes2=""
            fi

            if [ {Seqname_suffix_off} = True ]; then
                sed -i -E 's/(\|\|full([[:space:]]+)|\|\|[0-9]+_partial([[:space:]]+)|\|\|lt2gene([[:space:]]+))/\\2\\3\\4/;' {output.score}
                sed -i -E 's/(\|\|full$|\|\|[0-9]+_partial$|\|\|lt2gene$)//;' {output.fa} {output.boundary} 
                if [ {Prep_for_dramv} = True ]; then
                    sed -i -E 's/(__full(\|[0-9]+\|(c|l)$)|__[0-9]+_partial(\|[0-9]+\|(c|l)$)|__lt2gene(\|[0-9]+\|(c|l)$))/\\2\\4\\6/;'  {Label}for-dramv/viral-affi-contigs-for-dramv.tab
                    sed -i -E 's/(__full(__[0-9]+\|)|__[0-9]+_partial(__[0-9]+\|)|__lt2gene(__[0-9]+\|))/\\2\\3\\4/;' {Label}for-dramv/viral-affi-contigs-for-dramv.tab
                    sed -i -E 's/(__full(-cat_[1-6]$)|__[0-9]+_partial(-cat_[1-6]$)|__lt2gene(-cat_[1-6]$))/\\2\\3\\4/;' {Label}for-dramv/final-viral-combined-for-dramv.fa 
                fi
                Suffix_notes=""
            else
                Suffix_notes="
                Suffix is added to seq names in final-viral-combined.fa:
                full    seqs (>=2 genes) as viral:\t||full
                partial seqs (>=2 genes) as viral:\t||partial
                short   seqs (< 2 genes) as viral:\t||lt2gene
                $Dramv_notes2
                "
            fi

            printf "
            ====> VirSorter run (provirus mode) finished.
            # of full    seqs (>=2 genes) as viral:\t$N_lytic
            # of partial seqs (>=2 genes) as viral:\t$N_lysogenic
            # of short   seqs (< 2 genes) as viral:\t$N_lt2gene

            Useful output files:
                {Label}final-viral-score.tsv       ==> score table
                {Label}final-viral-combined.fa     ==> all viral seqs
                {Label}final-viral-boundary.tsv    ==> table with boundary info
                $Dramv_notes
            $Suffix_notes
            NOTES:
            Users can further screen the results based on the following 
                columns in {Label}final-viral-score.tsv:
                - contig length (length) 
                - hallmark gene count (hallmark)
                - viral gene %% (viral) 
                - cellular gene %% (cellular)
            The "group" field in {Label}final-viral-score.tsv should NOT be used
                as reliable taxonomy info
            We recommend this SOP/tutorial for quality control 
                (make sure to use the lastest version):
                https://dx.doi.org/10.17504/protocols.io.bwm5pc86

            <====
            " | python {Scriptdir}/echo.py
            """
# provirus off
else:
    localrules: get_viral_combined
    rule get_viral_combined:
        input: 
            fullseq=f'{Tmpdir}/viral-fullseq.fa',
            proba=f'{Tmpdir}/all-fullseq-proba.tsv',
            lt2gene=f'{Tmpdir}/viral-lt2gene-w-hallmark.fa',
        output: 
            combined=f'{Tmpdir}/viral-combined.fa',
            proba=f'{Tmpdir}/viral-combined-proba-more-cols.tsv',
        conda: '{}/vs2.yaml'.format(Conda_yaml_dir)
        shell:
            """
            python {Scriptdir}/filter-table-and-add-suffix-to-seqname.py --suffix "||full" {Tmpdir}/all-fullseq-proba.tsv {Tmpdir}/viral-fullseq.fa {Tmpdir}/viral-combined-proba.tsv
            python {Scriptdir}/add-suffix-seqname-keep-desc.py {Tmpdir}/viral-fullseq.fa "||full" > {Tmpdir}/viral-combined.fa
            python {Scriptdir}/add-suffix-seqname-keep-desc.py {input.lt2gene} "||lt2gene" >> {Tmpdir}/viral-combined.fa
            python {Scriptdir}/add-extra-to-table.py {Tmpdir}/viral-combined-proba.tsv {Tmpdir}/viral-combined.fa {Tmpdir}/viral-combined-proba-more-cols.tsv
            """


    checkpoint split_gff_by_group:
        input: f'{Tmpdir}/{{group}}/all.pdg.gff'
        output: directory(f'{Tmpdir}/{{group}}/all.pdg.gff.splitdir')
        conda: '{}/vs2.yaml'.format(Conda_yaml_dir)
        shell:
            """
            rm -f {Tmpdir}/{wildcards.group}/all.pdg.gff.splitdir/all.pdg.gff.*.anno
            rm -f {Tmpdir}/{wildcards.group}/all.pdg.gff.splitdir/all.pdg.gff.*.affi.tab
            python {Scriptdir}/split-gff-even-seqnum-per-file.py {input} {output} {Gff_seqnum_per_split}
            """

    rule make_annotation_table_by_group_by_split:
        input: 
            gff=f'{Tmpdir}/{{group}}/all.pdg.gff.splitdir/all.pdg.gff.{{idx}}.split',
            taxwhm=f'{Tmpdir}/{{group}}/all.pdg.hmm.taxwhm',
            taxpfam=f'{Tmpdir}/{{group}}/all.pdg.hmm.taxpfam',
            seqfile=f'{Tmpdir}/viral-combined.fa'
        output: 
            anno=temp(f'{Tmpdir}/{{group}}/all.pdg.gff.splitdir/all.pdg.gff.{{idx}}.anno'),
            affi=temp(f'{Tmpdir}/{{group}}/all.pdg.gff.splitdir/all.pdg.gff.{{idx}}.affi.tab'),
        conda: '{}/vs2.yaml'.format(Conda_yaml_dir)
        shell:
            """
            python {Scriptdir}/make-affi-contigs-tabfile-on-split-gff.py --pfamtax {input.taxpfam} {input.seqfile} {output.anno} {output.affi} {input.gff} {input.taxwhm} {wildcards.group}
            """

    def merge_annotation_table_by_group_from_split_input_agg(wildcards):
        split_dir = checkpoints.split_gff_by_group.get(
                        **wildcards).output[0]
        idx_lis = glob_wildcards(
                '{}/all.pdg.gff.{{idx}}.split'.format(split_dir)).idx
        anno_str = '{}/all.pdg.gff.{{idx}}.anno'.format(split_dir)
        affi_str = '{}/all.pdg.gff.{{idx}}.affi.tab'.format(split_dir)
        anno_lis = expand(anno_str, idx=idx_lis)
        affi_lis = expand(affi_str, idx=idx_lis)
        return {'anno': anno_lis, 'affi': affi_lis}

    localrules: merge_annotation_table_by_group_from_split
    rule merge_annotation_table_by_group_from_split:
        input: unpack(merge_annotation_table_by_group_from_split_input_agg)
        output: 
            anno=f'{Tmpdir}/{{group}}/all.pdg.anno',
            affi=f'{Tmpdir}/{{group}}/all.pdg.affi',
        conda: '{}/vs2.yaml'.format(Conda_yaml_dir)
        shell:
            """
            printf "%s\n" {input.anno} | xargs cat | awk '!/^seqname\t/ || !f++' > {output.anno}
            printf "%s\n" {input.affi} | xargs cat | awk '!/^seqname\t/ || !f++' > {output.affi}
            """

    localrules: merge_annotation_table_from_groups
    rule merge_annotation_table_from_groups:
        input: 
            anno=expand(f'{Tmpdir}/{{group}}/all.pdg.anno', group=Groups),
            affi=expand(f'{Tmpdir}/{{group}}/all.pdg.affi', group=Groups),
        output:
            anno=f'{Tmpdir}/viral.anno',
            affi=f'{Tmpdir}/viral-affi-contigs-for-dramv.tab',
        conda: '{}/vs2.yaml'.format(Conda_yaml_dir)
        shell:
            """
            printf "%s\n" {input.anno} | xargs cat | awk '!/^seqname\t/ || !f++' > {output.anno}
            printf "%s\n" {input.affi} | xargs cat | awk '!/^seqname\t/ || !f++' > {output.affi}
            """

    #localrules: make_annotation_table
    #rule make_annotation_table:
    #    input: 
    #        f'{Tmpdir}/viral-combined.fa',
    #        expand(f'{Tmpdir}/{{group}}/all.pdg.gff', group=Groups),
    #        expand(f'{Tmpdir}/{{group}}/all.pdg.hmm.taxwhm', group=Groups),
    #        expand(f'{Tmpdir}/{{group}}/all.pdg.hmm.taxpfam', group=Groups),
    #    output: 
    #        anno=f'{Tmpdir}/viral.anno',
    #        affi=f'{Tmpdir}/viral-affi-contigs-for-dramv.tab',
    #    conda: '{}/vs2.yaml'.format(Conda_yaml_dir)
    #    params:
    #        gff_fs_str = ','.join(
    #            [f'{Tmpdir}/{group}/all.pdg.gff' for group in Groups]
    #        ),
    #        taxwhm_fs_str = ','.join(
    #            [f'{Tmpdir}/{group}/all.pdg.hmm.taxwhm' for group in Groups]
    #        ),
    #        taxpfam_fs_str = ','.join(
    #            [f'{Tmpdir}/{group}/all.pdg.hmm.taxpfam' for group in Groups]
    #        ),
    #        groups_str = ','.join(Groups),
    #    shell:
    #        """
    #        python {Scriptdir}/make-affi-contigs-tabfile.py --pfamtax-list-str {params.taxpfam_fs_str} {Tmpdir}/viral-combined.fa {output.anno} {output.affi} {params.gff_fs_str} {params.taxwhm_fs_str} {params.groups_str}
    #        """

    def finalize_input_agg(wildcards):
        if Prep_for_dramv:
            d = dict(
                    seqfile=f'{Tmpdir}/viral-combined.fa',
                    proba=f'{Tmpdir}/viral-combined-proba-more-cols.tsv',
                    affi=f'{Tmpdir}/viral-affi-contigs-for-dramv.tab',
            )
        else:
            d = dict(
                    seqfile=f'{Tmpdir}/viral-combined.fa',
                    proba=f'{Tmpdir}/viral-combined-proba-more-cols.tsv',
            )
        return d

    localrules: finalize
    rule finalize:
        input: 
            unpack(finalize_input_agg)
        output: 
            fa=f'{Label}final-viral-combined.fa',
            score=f'{Label}final-viral-score.tsv',
        conda: '{}/vs2.yaml'.format(Conda_yaml_dir)
        shell:
            """
            echo {Tmpdir}/*/all.pdg.gff.splitdir/all.pdg.gff.*.split | xargs rm -f
            python {Scriptdir}/filter-score-table.py config.yaml {Tmpdir}/viral-combined-proba-more-cols.tsv {Tmpdir}/viral-combined.fa {output.score} {output.fa}

            if [ {Prep_for_dramv} = "True" ]; then
                mkdir -p {Label}for-dramv
                python {Scriptdir}/modify-seqname-for-dramv.py {output.fa} {output.score} -o {Label}for-dramv/final-viral-combined-for-dramv.fa
                cp {Tmpdir}/viral-affi-contigs-for-dramv.tab {Label}for-dramv
            fi

            N_viral_fullseq=$(grep -c '^>.*||full$' {output.fa} || :)
            N_viral_lt2gene=$(grep -c '^>.*||lt2gene$' {output.fa} || :)
            if [ {Prep_for_dramv} = True ]; then
                Dramv_notes="{Label}for-dramv                  ==> dir with input files for dramv"
                Dramv_notes2="For seqnames in files for dramv, 
                    | is replaced with _ to be compatible with DRAMv"
            else
                Dramv_notes=""
                Dramv_notes2=""
            fi
            if [ {Seqname_suffix_off} = True ]; then
                sed -i -E 's/(\|\|full([[:space:]]+)|\|\|[0-9]+_partial([[:space:]]+)|\|\|lt2gene([[:space:]]+))/\\2\\3\\4/;' {output.score}
                sed -i -E 's/(\|\|full$|\|\|[0-9]+_partial$|\|\|lt2gene$)//;' {output.fa}
                if [ {Prep_for_dramv} = True ]; then
                    sed -i -E 's/(__full(\|[0-9]+\|(c|l)$)|__[0-9]+_partial(\|[0-9]+\|(c|l)$)|__lt2gene(\|[0-9]+\|(c|l)$))/\\2\\4\\6/;'  {Label}for-dramv/viral-affi-contigs-for-dramv.tab
                    sed -i -E 's/(__full(__[0-9]+\|)|__[0-9]+_partial(__[0-9]+\|)|__lt2gene(__[0-9]+\|))/\\2\\3\\4/;' {Label}for-dramv/viral-affi-contigs-for-dramv.tab
                    sed -i -E 's/(__full(-cat_[1-6]$)|__[0-9]+_partial(-cat_[1-6]$)|__lt2gene(-cat_[1-6]$))/\\2\\3\\4/;' {Label}for-dramv/final-viral-combined-for-dramv.fa 
                fi
                Suffix_notes=""
            else
                Suffix_notes="
                Suffix is added to seq names in {output.fa}:
                contigs (>=2 genes) as viral:\t||full
                contigs (< 2 genes) as viral:\t||lt2gene
                $Dramv_notes2
                "
            fi
                
            printf "
            ====> VirSorter run (non-provirus mode) finished.
            # of contigs w/ >=2 genes as viral:\t$N_viral_fullseq
            # of contigs w/ < 2 genes as viral:\t$N_viral_lt2gene

            Useful output files:
            {output.score}      ==> score table
            {output.fa}    ==> all viral seqs
            $Dramv_notes
            $Suffix_notes

            NOTES: 
            Users can further screen the results based on the 
                following columns in {output.score}
                - contig length (length) 
                - hallmark gene count (hallmark)
                - viral gene %% (viral) 
                - cellular gene %% (cellular)
            The "group" field in {output.score} should NOT be used
                as reliable taxonomy info
            We recommend this SOP/tutorial for quality control 
                (make sure to use the lastest version):
                https://dx.doi.org/10.17504/protocols.io.bwm5pc86

            <====
            " | python {Scriptdir}/echo.py
            

            """
