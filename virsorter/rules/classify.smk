
rule classify_by_group:
    input: 'iter-0/{group}/all.pdg.ftr'
    output: 'iter-0/{group}/all.pdg.clf'
    conda: '{}/vs2.yaml'.format(Conda_yaml_dir)
    shell:
        """
        Log={Wkdir}/log/iter-0/step3-classify/all-score-{wildcards.group}.log
        python {Scriptdir}/classify.py {input} {Dbdir}/group/{wildcards.group}/model {wildcards.group} {output} 2> $Log || {{ echo "See error details in $Log" | python {Scriptdir}/echo.py --level error; exit 1; }}
        """

localrules: merge_classification
rule merge_classification:
    input: expand('iter-0/{group}/all.pdg.clf', group=Groups)
    output: 
        clf='iter-0/all-fullseq-proba.tsv',
    conda: '{}/vs2.yaml'.format(Conda_yaml_dir)
    params:
        clf_fs_str = ' '.join(
            ['iter-0/{}/all.pdg.clf'.format(group) for group in Groups]
        ),
    shell:
        """
        Log={Wkdir}/log/iter-0/step3-classify/all-score-merge.log
        ### merge clf from all groups
        python {Scriptdir}/merge-clf-from-groups.py {output.clf} {params.clf_fs_str} 2> $Log || {{ echo "See error details in $Log" | python {Scriptdir}/echo.py --level error; exit 1; }}
        """

localrules: pick_viral_fullseq
rule pick_viral_fullseq:
    input:
        clf='iter-0/all-fullseq-proba.tsv',
    output: 
        contig='iter-0/viral-fullseq.fa',
        hmk_cnt='iter-0/all-hallmark-cnt.tsv',
        lt2gene='iter-0/viral-lt2gene-w-hallmark.fa',
    conda: '{}/vs2.yaml'.format(Conda_yaml_dir)
    params:
        group_str = ','.join(
            [group for group in Groups if os.path.exists('{}/group/{}/hallmark-gene.list'.format(Dbdir, group))]
        ),
        hmk_fs_str = ','.join(
            ['{}/group/{}/hallmark-gene.list'.format(Dbdir, group) for group in Groups if os.path.exists('{}/group/{}/hallmark-gene.list'.format(Dbdir, group))]
        ),
        tax_fs_str = ','.join(
            ['iter-0/{}/all.pdg.hmm.tax'.format(group) for group in Groups if os.path.exists('{}/group/{}/hallmark-gene.list'.format(Dbdir, group))]
        ),
        all_group_str = ','.join(Groups),
        all_tax_fs_str = ','.join(
            ['iter-0/{}/all.pdg.hmm.tax'.format(group) for group in Groups]
        ),
        all_gff_fs_str = ','.join(
            ['iter-0/{}/all.pdg.gff'.format(group) for group in Groups]
        ),
        all_ftr_fs_str = ','.join(
            ['iter-0/{}/all.pdg.ftr'.format(group) for group in Groups]
        ),
    shell:
        """
        Log={Wkdir}/log/iter-0/step3-classify/pick-viral-fullseq.log
        python {Scriptdir}/pick-viral-contig-from-clf.py {Proba_cutoff} {input.clf} iter-0/all.fna > {output.contig}.tmp 2> $Log || {{ echo "See error details in $Log" | python {Scriptdir}/echo.py --level error; exit 1; }}

        python {Scriptdir}/add-extra-to-fullseq-fasta-header.py {output.contig}.tmp {params.all_ftr_fs_str} {params.all_group_str} > {output.contig} 2>> $Log || {{ echo "See error details in $Log" | python {Scriptdir}/echo.py --level error; exit 1; }}
        rm {output.contig}.tmp

        python {Scriptdir}/get-hallmark-cnt-for-each-seq.py {output.hmk_cnt} "{params.group_str}" "{params.hmk_fs_str}" "{params.tax_fs_str}" 2>> $Log || {{ echo "See error details in $Log" | python {Scriptdir}/echo.py --level error; exit 1; }}

        python {Scriptdir}/get-seq-w-lt2gene-w-hallmark.py {output.hmk_cnt} {input.clf} iter-0/all.fna > {output.lt2gene}.tmp 2>> $Log || {{ echo "See error details in $Log" | python {Scriptdir}/echo.py --level error; exit 1; }}

        python {Scriptdir}/add-extra-to-lt2gene-fasta-header.py {output.lt2gene}.tmp {params.all_gff_fs_str} {params.all_tax_fs_str} {params.all_group_str} > {output.lt2gene} 2>> $Log || {{ echo "See error details in $Log" | python {Scriptdir}/echo.py --level error; exit 1; }}
        rm {output.lt2gene}.tmp

        """

if Provirus:
    checkpoint split_gff_by_group:
        input: 'iter-0/{group}/all.pdg.gff'
        output: directory('iter-0/{group}/all.pdg.gff.splitdir')
        shell:
            """
            python {Scriptdir}/split-gff-even-seqnum-per-file.py {input} {output} {Gff_seqnum_per_split}
            """

    rule provirus_call_by_group_by_split:
        input: 
            gff='iter-0/{group}/all.pdg.gff.splitdir/all.pdg.gff.{idx}.split',
            tax='iter-0/{group}/all.pdg.hmm.tax',
            clf='iter-0/all-fullseq-proba.tsv',
        output: 
            boundry=temp('iter-0/{group}/all.pdg.gff.splitdir/all.pdg.gff.{idx}.split.prv.bdy'),
            ftr=temp('iter-0/{group}/all.pdg.gff.splitdir/all.pdg.gff.{idx}.split.prv.ftr'),
        conda: '{}/vs2.yaml'.format(Conda_yaml_dir)
        shell:
            """
            Log={Wkdir}/{input.gff}.prv.log
            Hallmark_list_f={Dbdir}/group/{wildcards.group}/hallmark-gene.list
            if [ -s $Hallmark_list_f ]; then
                python {Scriptdir}/provirus.py {input.gff} {input.tax} {Dbdir}/rbs/rbs-catetory.tsv {Dbdir}/group/{wildcards.group}/model {output.boundry} {output.ftr} --fullseq-clf {input.clf} --group {wildcards.group} --proba {Proba_cutoff} --hallmark {Dbdir}/group/{wildcards.group}/hallmark-gene.list 2> $Log || {{ echo "See error details in $Log" | python {Scriptdir}/echo.py --level error; exit 1; }}
            else
                python {Scriptdir}/prophage.py {input.gff} {input.tax} {Dbdir}/rbs/rbs-catetory.tsv {Dbdir}/group/{wildcards.group}/model {output.boundry} {output.ftr} --fullseq-clf {input.clf} --group {wildcards.group} --proba {Proba_cutoff} 2> $Log || {{ echo "See error details in $Log" | python {Scriptdir}/echo.py --level error; exit 1; }}
            fi
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
            bdy='iter-0/{group}/all.pdg.prv.bdy',
            ftr='iter-0/{group}/all.pdg.prv.ftr',
        conda: '{}/vs2.yaml'.format(Conda_yaml_dir)
        shell:
            """
            cat {input.bdy} | awk '!/^seqname\t/ || !f++' > {output.bdy}
            cat {input.ftr} | awk '!/^seqname\t/ || !f++' > {output.ftr}
            """

    localrules: merge_provirus_call_from_groups
    rule merge_provirus_call_from_groups:
        input: expand('iter-0/{group}/all.pdg.prv.bdy', group=Groups)
        output:
            partial='iter-0/viral-partseq.tsv',
            full='iter-0/viral-fullseq.tsv',
        conda: '{}/vs2.yaml'.format(Conda_yaml_dir)
        params:
            prv_fs_str = ','.join(
                ['iter-0/{}/all.pdg.prv.bdy'.format(group) for group in Groups]
            ),
            groups_str = ','.join(Groups),
        shell:
            """
            Log={Wkdir}/log/iter-0/step3-classify/provirus-score-merge.log
            python {Scriptdir}/merge-provirus-from-groups.py {params.prv_fs_str} {params.groups_str} {output.partial} {output.full} 2> $Log || {{ echo "See error details in $Log" | python {Scriptdir}/echo.py --level error; exit 1; }}
            """

    localrules: extract_provirus_seqs
    rule extract_provirus_seqs:
        input:
            #contig={Seqfile},
            contig='iter-0/all.fna',
            full='iter-0/viral-fullseq.tsv',
            partial='iter-0/viral-partseq.tsv',
            hmk_cnt='iter-0/all-hallmark-cnt.tsv',
            lt2gene='iter-0/viral-lt2gene-w-hallmark.fa',
        output:
            fullseq='iter-0/viral-fullseq-trim.fa',
            partial='iter-0/viral-partseq.fa',
            combined='iter-0/viral-combined.fa',
        conda: '{}/vs2.yaml'.format(Conda_yaml_dir)
        shell:
            """
            Log={Wkdir}/log/iter-0/step3-classify/provirus-seq-extract.log
            python {Scriptdir}/extract-provirus-seqs.py {input.contig} {input.full} {input.partial} {output.fullseq} {output.partial} 2> $Log || {{ echo "See error details in $Log" | python {Scriptdir}/echo.py --level error; exit 1; }}

            cat {output.fullseq} {output.partial} > {output.combined}
            python {Scriptdir}/add-suffix-seqname-keep-desc.py {input.lt2gene} "||lt2gene" >> {output.combined} 2>> $Log || {{ echo "See error details in $Log" | python {Scriptdir}/echo.py --level error; exit 1; }}

            """

    ##### add score for non-best-scoring groups
    rule classify_full_and_part_by_group:
        input: 
            gff='iter-0/{group}/all.pdg.gff',
            tax='iter-0/{group}/all.pdg.hmm.tax',
            seqfile='iter-0/viral-combined.fa',
        output: 'iter-0/{group}/viral.trim.clf'
        conda: '{}/vs2.yaml'.format(Conda_yaml_dir)
        shell:
            """
            Log={Wkdir}/log/iter-0/step3-classify/provirus-score-{wildcards.group}.log
            Hallmark_list_f={Dbdir}/group/{wildcards.group}/hallmark-gene.list
            if [ -s $Hallmark_list_f ]; then
                python {Scriptdir}/classify-trimed.py {input.gff} {input.tax} {Dbdir}/rbs/rbs-catetory.tsv {Dbdir}/group/{wildcards.group}/model {input.seqfile} {output} --group {wildcards.group} --hallmark {Dbdir}/group/{wildcards.group}/hallmark-gene.list 2> $Log || {{ echo "See error details in $Log" | python {Scriptdir}/echo.py --level error; exit 1; }}
            else
                python {Scriptdir}/classify-trimed.py {input.gff} {input.tax} {Dbdir}/rbs/rbs-catetory.tsv {Dbdir}/group/{wildcards.group}/model {input.seqfile} {output} --group {wildcards.group} 2> $Log || {{ echo "See error details in $Log" | python {Scriptdir}/echo.py --level error; exit 1; }}
            fi

            """

    localrules: merge_full_and_part_classification
    rule merge_full_and_part_classification:
        input: expand('iter-0/{group}/viral.trim.clf', group=Groups)
        output: 
            clf='iter-0/viral-combined-proba.tsv',
        conda: '{}/vs2.yaml'.format(Conda_yaml_dir)
        params:
            clf_fs_str = ' '.join(
                ['iter-0/{}/viral.trim.clf'.format(group) for group in Groups]
            ),
        shell:
            """
            Log={Wkdir}/log/iter-0/step3-classify/provirus-score-merge.log
            ### merge clf from all groups
            python {Scriptdir}/merge-clf-trim-from-groups.py {output.clf} {params.clf_fs_str} 2> $Log || {{ echo "See error details in $Log" | python {Scriptdir}/echo.py --level error; exit 1; }}
            """

    localrules: finalize
    rule finalize:
        input:
            'iter-0/viral-fullseq-trim.fa',
            'iter-0/viral-partseq.fa',
            'iter-0/viral-combined-proba.tsv',
        output: 
            'final-viral-score.tsv',
            'final-viral-combined.fa',
            'final-viral-boundary.tsv',
        conda: '{}/vs2.yaml'.format(Conda_yaml_dir)
        shell:
            """
            python {Scriptdir}/add-extra-to-table.py iter-0/viral-combined-proba.tsv iter-0/viral-combined.fa iter-0/viral-combined-proba-more-cols.tsv
            python {Scriptdir}/filter-score-table.py config.yaml iter-0/viral-combined-proba-more-cols.tsv iter-0/viral-combined.fa final-viral-score.tsv final-viral-combined.fa
            cp iter-0/viral-fullseq.tsv final-viral-boundary.tsv
            grep -v '^seqname' iter-0/viral-partseq.tsv >> final-viral-boundary.tsv || : 
            N_lt2gene=$(grep -c '^>.*||lt2gene' final-viral-combined.fa || :)
            N_lytic=$(grep -c '^>.*||full' final-viral-combined.fa || :)
            N_lysogenic=$(grep -c '^>.+||.*_partial' final-viral-combined.fa || :)
            printf "
            ====> VirSorter run (provirus mode) finished.
            # of full    seqs (>=2 genes) as viral:\t$N_lytic
            # of partial seqs (>=2 genes) as viral:\t$N_lysogenic
            # of short   seqs (< 2 genes) as viral:\t$N_lt2gene

            Useful output files:
            final-viral-score.tsv       ==> score table
            final-viral-combined.fa     ==> all viral seqs
            final-viral-boundary.tsv    ==> table with boundary info
            
            Suffix is added to seq names in final-viral-combined.fa:
            full    seqs (>=2 genes) as viral:\t||full
            partial seqs (>=2 genes) as viral:\t||partial
            short   seqs (< 2 genes) as viral:\t||lt2gene

            NOTES:
            Users can further screen the results based on the following 
                columns in final-viral-score.tsv:
                - contig length (length) 
                - hallmark gene count (hallmark)
                - viral gene %% (viral) 
                - cellular gene %% (cellular)
            The "group" field in final-viral-score.tsv should NOT be used
                as reliale taxonomy info

            <====
            " | python {Scriptdir}/echo.py
            """
# provirus off
else:
    localrules: finalize
    rule finalize:
        input: 
            fullseq='iter-0/viral-fullseq.fa',
            proba='iter-0/all-fullseq-proba.tsv',
            lt2gene='iter-0/viral-lt2gene-w-hallmark.fa',
        output: 
            combined='final-viral-combined.fa',
            proba='final-viral-score.tsv',
        conda: '{}/vs2.yaml'.format(Conda_yaml_dir)
        shell:
            """
            python {Scriptdir}/filter-table-and-add-suffix-to-seqname.py "||full" iter-0/all-fullseq-proba.tsv iter-0/viral-fullseq.fa iter-0/viral-combined-proba.tsv
            python {Scriptdir}/add-suffix-seqname-keep-desc.py iter-0/viral-fullseq.fa "||full" > iter-0/viral-combined.fa
            python {Scriptdir}/add-suffix-seqname-keep-desc.py {input.lt2gene} "||lt2gene" >> iter-0/viral-combined.fa

            python {Scriptdir}/add-extra-to-table.py iter-0/viral-combined-proba.tsv iter-0/viral-combined.fa iter-0/viral-combined-proba-more-cols.tsv
            python {Scriptdir}/filter-score-table.py config.yaml iter-0/viral-combined-proba-more-cols.tsv iter-0/viral-combined.fa final-viral-score.tsv final-viral-combined.fa
            N_viral_fullseq=$(grep -c '^>.*||full' final-viral-combined.fa || :)
            N_viral_lt2gene=$(grep -c '^>.*||lt2gene' final-viral-combined.fa || :)
            printf "
            ====> VirSorter run (non-provirus mode) finished.
            # of contigs w/ >=2 genes as viral:\t$N_viral_fullseq
            # of contigs w/ < 2 genes as viral:\t$N_viral_lt2gene

            Useful output files:
            final-viral-score.tsv      ==> score table
            final-viral-combined.fa    ==> all viral seqs

            Suffix is added to seq names in final-viral-combined.fa:
            contigs (>=2 genes) as viral:\t||full
            contigs (< 2 genes) as viral:\t||lt2gene

            NOTES: 
            Users can further screen the results based on the 
                following columns in final-viral-score.tsv
                - contig length (length) 
                - hallmark gene count (hallmark)
                - viral gene %% (viral) 
                - cellular gene %% (cellular)
            The "group" field in final-viral-score.tsv should NOT be used
                as reliale taxonomy info

            <====
            " | python {Scriptdir}/echo.py
            """
