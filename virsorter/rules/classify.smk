
rule classify_by_group:
    input: 'iter-0/{group}/all.pdg.ftr'
    output: 'iter-0/{group}/all.pdg.clf'
    ##conda: 'envs/sklearn0212.yml'
    #conda: '{}/vs2.yaml'.format(Conda_yaml_dir)
    shell:
        """
        Log={Wkdir}/log/iter-0/step3-classify/classify-{wildcards.group}.log
        python {Scriptdir}/classify.py {input} {Dbdir}/group/{wildcards.group}/model {wildcards.group} {output} 2> $Log || {{ echo "See error details in $Log" | python {Scriptdir}/echo.py --level error; exit 1; }}
        """

localrules: merge_classification
rule merge_classification:
    input: expand('iter-0/{group}/all.pdg.clf', group=Groups)
    output: 
        clf='iter-0/all-fullseq-proba.tsv',
    #conda: '{}/vs2.yaml'.format(Conda_yaml_dir)
    params:
        clf_fs_str = ' '.join(
            ['iter-0/{}/all.pdg.clf'.format(group) for group in Groups]
        ),
    shell:
        """
        Log={Wkdir}/log/iter-0/step3-classify/classify-merged.log
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
    #conda: '{}/vs2.yaml'.format(Conda_yaml_dir)
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
    shell:
        """
        Log={Wkdir}/log/iter-0/step3-classify/classify-merged.log
        python {Scriptdir}/pick-viral-contig-from-clf.py {Proba_cutoff} {input.clf} iter-0/all.fna > iter-0/viral-fullseq-contig.fa.tmp 2>> $Log || {{ echo "See error details in $Log" | python {Scriptdir}/echo.py --level error; exit 1; }}

        python {Scriptdir}/get-hallmark-cnt-for-each-seq.py {output.hmk_cnt} "{params.group_str}" "{params.hmk_fs_str}" "{params.tax_fs_str}" 2>> $Log || {{ echo "See error details in $Log" | python {Scriptdir}/echo.py --level error; exit 1; }}

        python {Scriptdir}/get-seq-w-lt2genes-w-hallmark.py {output.hmk_cnt} {input.clf} iter-0/all.fna > {output.lt2gene} 2>> $Log || {{ echo "See error details in $Log" | python {Scriptdir}/echo.py --level error; exit 1; }}
        if [ {Hallmark_required_on_short} = "True" ]; then
            python {Scriptdir}/remove-short-seq-wo-hallmark.py {output.hmk_cnt} iter-0/viral-fullseq-contig.fa.tmp > {output.contig} 2>> $Log || {{ echo "See error details in $Log" | python {Scriptdir}/echo.py --level error; exit 1; }}
            rm iter-0/viral-fullseq-contig.fa.tmp
        else
            mv iter-0/viral-fullseq-contig.fa.tmp {output.contig}
        fi
        """

if Provirus:
    rule provirus_call_by_group:
        input: 
            gff='iter-0/{group}/all.pdg.gff',
            tax='iter-0/{group}/all.pdg.hmm.tax',
            clf='iter-0/all-fullseq-proba.tsv',
        output: 
            boundry='iter-0/{group}/all.pdg.prv.bdy',
            ftr='iter-0/{group}/all.pdg.prv.ftr',
        #conda: '{}/vs2.yaml'.format(Conda_yaml_dir)
        shell:
            """
            Log={Wkdir}/log/iter-0/step3-classify/classify-{wildcards.group}.log
            Hallmark_list_f={Dbdir}/group/{wildcards.group}/hallmark-gene.list
            if [ -s $Hallmark_list_f ]; then
                python {Scriptdir}/provirus.py {input.gff} {input.tax} {Dbdir}/rbs/rbs-catetory.tsv {Dbdir}/group/{wildcards.group}/model {output.boundry} {output.ftr} --fullseq-clf {input.clf} --group {wildcards.group} --proba {Proba_cutoff} --hallmark {Dbdir}/group/{wildcards.group}/hallmark-gene.list 2>> $Log || {{ echo "See error details in $Log" | python {Scriptdir}/echo.py --level error; exit 1; }}
            else
                python {Scriptdir}/prophage.py {input.gff} {input.tax} {Dbdir}/rbs/rbs-catetory.tsv {Dbdir}/group/{wildcards.group}/model {output.boundry} {output.ftr} --fullseq-clf {input.clf} --group {wildcards.group} --proba {Proba_cutoff} 2>> $Log || {{ echo "See error details in $Log" | python {Scriptdir}/echo.py --level error; exit 1; }}
            fi
            """

    localrules: merge_provirus_call_from_groups
    rule merge_provirus_call_from_groups:
        input: expand('iter-0/{group}/all.pdg.prv.bdy', group=Groups)
        output:
            partial='iter-0/viral-partseq.tsv',
            full='iter-0/viral-fullseq.tsv',
        #conda: '{}/vs2.yaml'.format(Conda_yaml_dir)
        params:
            prv_fs_str = ','.join(
                ['iter-0/{}/all.pdg.prv.bdy'.format(group) for group in Groups]
            ),
            groups_str = ','.join(Groups),
        shell:
            """
            Log={Wkdir}/log/iter-0/step3-classify/classify-merged.log
            python {Scriptdir}/merge-provirus-from-groups.py {params.prv_fs_str} {params.groups_str} {output.partial} {output.full} 2>> $Log || {{ echo "See error details in $Log" | python {Scriptdir}/echo.py --level error; exit 1; }}
            """

    localrules: extract_provirus_seqs
    rule extract_provirus_seqs:
        input:
            #contig={Seqfile},
            contig='iter-0/all.fna',
            full='iter-0/viral-fullseq.tsv',
            partial='iter-0/viral-partseq.tsv',
            hmk_cnt='iter-0/all-hallmark-cnt.tsv',
        output:
            fullseq='iter-0/viral-fullseq-trim.fa',
            partial='iter-0/viral-partseq.fa',
        #conda: '{}/vs2.yaml'.format(Conda_yaml_dir)
        shell:
            """
            Log={Wkdir}/log/iter-0/step3-classify/classify-merged.log
            python {Scriptdir}/extract-provirus-seqs.py {input.contig} {input.full} {input.partial} iter-0/viral-fullseq-trim.fa.tmp {output.partial} 2>> $Log || {{ echo "See error details in $Log" | python {Scriptdir}/echo.py --level error; exit 1; }}

            if [ {Hallmark_required_on_short} = 'true' ]; then
                python {Scriptdir}/remove-short-seq-wo-hallmark.py {input.hmk_cnt} iter-0/viral-fullseq-trim.fa.tmp > {output.fullseq} 2>> $Log || {{ echo "See error details in $Log" | python {Scriptdir}/echo.py --level error; exit 1; }}
                rm iter-0/viral-fullseq-trim.fa.tmp
            else
                mv iter-0/viral-fullseq-trim.fa.tmp {output.fullseq}
            fi
            """

    localrules: finalize
    rule finalize:
        input:
            'iter-0/viral-fullseq-trim.fa',
            'iter-0/viral-partseq.fa',
        output: 
            'final-viral-fullseq-trim.fa',
            'final-viral-partseq.fa',
            'final-viral-boundary.tsv',
        shell:
            """
            cp iter-0/viral-fullseq.fa final-viral-fullseq.fa
            cp iter-0/all-fullseq-proba.tsv final-all-fullseq-proba.tsv
            cp iter-0/viral-fullseq-trim.fa final-viral-fullseq-trim.fa
            cp iter-0/viral-partseq.fa final-viral-partseq.fa
            cp iter-0/viral-fullseq.tsv final-viral-boundary.tsv
            cat iter-0/viral-fullseq-trim.fa iter-0/viral-partseq.fa > final-viral-combined.fa
            grep -v '^seqname' iter-0/viral-partseq.tsv >> final-viral-boundary.tsv || : 
            N_lytic=$(grep -c '^>' final-viral-fullseq-trim.fa || :)
            N_lysogenic=$(grep -c '^>' final-viral-partseq.fa || :)
            echo -e "
            ====> VirSorter run (provirus mode) finished.
            # of full    seqs as viral:\t$N_lytic
            # of partial seqs as viral:\t$N_lysogenic

            Useful output files:
            final-viral-fullseq-trim.fa
            final-viral-partial.fa
            final-viral-boundary.tsv
            <====
            " | python {Scriptdir}/echo.py
            """
else:
    localrules: finalize
    rule finalize:
        input: 
            'iter-0/viral-fullseq.fa',
            'iter-0/all-fullseq-proba.tsv',
        output: 
            'final-viral-fullseq.fa',
            'final-all-fullseq-proba.tsv',
        shell:
            """
            cp iter-0/viral-fullseq.fa final-viral-fullseq.fa
            cp iter-0/all-fullseq-proba.tsv final-all-fullseq-proba.tsv
            ln -s iter-0/viral-fullseq.fa final-viral-combined.fa
            N_viral_fullseq=$(grep -c '^>' final-viral-fullseq.fa || :)
            echo -e "
            ====> VirSorter run (non-provirus mode) finished.
            # of viral contigs (whole seq classified as viral):\t$N_viral_fullseq

            Useful output files:
            final-viral-fullseq.fa
            final-all-fullseq-proba.tsv
            <====
            " | python {Scriptdir}/echo.py
            """
