
rule classify_by_group:
    input: 'iter-0/{group}/all.pdg.ftr'
    output: 'iter-0/{group}/all.pdg.clf'
    conda: '{}/vs2.yaml'.format(Conda_yaml_dir)
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
    conda: '{}/vs2.yaml'.format(Conda_yaml_dir)
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
    shell:
        """
        Log={Wkdir}/log/iter-0/step3-classify/classify-merged.log
        python {Scriptdir}/pick-viral-contig-from-clf.py {Proba_cutoff} {input.clf} iter-0/all.fna > iter-0/viral-fullseq-contig.fa.tmp 2>> $Log || {{ echo "See error details in $Log" | python {Scriptdir}/echo.py --level error; exit 1; }}

        python {Scriptdir}/get-hallmark-cnt-for-each-seq.py {output.hmk_cnt} "{params.group_str}" "{params.hmk_fs_str}" "{params.tax_fs_str}" 2>> $Log || {{ echo "See error details in $Log" | python {Scriptdir}/echo.py --level error; exit 1; }}

        python {Scriptdir}/get-seq-w-lt2gene-w-hallmark.py {output.hmk_cnt} {input.clf} iter-0/all.fna > {output.lt2gene} 2>> $Log || {{ echo "See error details in $Log" | python {Scriptdir}/echo.py --level error; exit 1; }}
        if [ {Hallmark_required_on_short} = "true" ]; then
            python {Scriptdir}/remove-short-seq-wo-hallmark.py {output.hmk_cnt} iter-0/viral-fullseq-contig.fa.tmp > {output.contig} 2>> $Log || {{ echo "See error details in $Log" | python {Scriptdir}/echo.py --level error; exit 1; }}
            rm iter-0/viral-fullseq-contig.fa.tmp
        else
            mv iter-0/viral-fullseq-contig.fa.tmp {output.contig}
        fi
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
            Log={input.gff}.prv.log
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
            lt2gene='iter-0/viral-lt2gene-w-hallmark.fa',
        output:
            fullseq='iter-0/viral-fullseq-trim.fa',
            partial='iter-0/viral-partseq.fa',
            combined='iter-0/viral-combined.fa',
        conda: '{}/vs2.yaml'.format(Conda_yaml_dir)
        shell:
            """
            Log={Wkdir}/log/iter-0/step3-classify/classify-merge.log
            python {Scriptdir}/extract-provirus-seqs.py {input.contig} {input.full} {input.partial} iter-0/viral-fullseq-trim.fa.tmp {output.partial} 2>> $Log || {{ echo "See error details in $Log" | python {Scriptdir}/echo.py --level error; exit 1; }}

            if [ {Hallmark_required_on_short} = "true" ]; then
                python {Scriptdir}/remove-short-seq-wo-hallmark.py {input.hmk_cnt} iter-0/viral-fullseq-trim.fa.tmp > {output.fullseq} 2>> $Log || {{ echo "See error details in $Log" | python {Scriptdir}/echo.py --level error; exit 1; }}
                rm iter-0/viral-fullseq-trim.fa.tmp
            else
                mv iter-0/viral-fullseq-trim.fa.tmp {output.fullseq}
            fi
            cat {output.fullseq} {output.partial} > {output.combined}
            python {Scriptdir}/add-suffix-seqname-keep-desc.py {input.lt2gene} "||lt2gene" >> {output.combined}

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
            Log={Wkdir}/log/iter-0/step3-classify/classify-trimmed-{wildcards.group}.log
            Hallmark_list_f={Dbdir}/group/{wildcards.group}/hallmark-gene.list
            if [ -s $Hallmark_list_f ]; then
                python {Scriptdir}/classify-trimed.py {input.gff} {input.tax} {Dbdir}/rbs/rbs-catetory.tsv {Dbdir}/group/{wildcards.group}/model {input.seqfile} {output} --group {wildcards.group} --hallmark {Dbdir}/group/{wildcards.group}/hallmark-gene.list 2> $Log || {{ echo "See error details in $Log" | python {Scriptdir}/echo.py --level error; exit 1; }}
            else
                python {Scriptdir}/classify-trimed.py {input.gff} {input.tax} {Dbdir}/rbs/rbs-catetory.tsv {Dbdir}/group/{wildcards.group}/model {input.seqfile} {output} --group {wildcards.group} 2> $Log || {{ echo "See error details in $Log" | python {Scriptdir}/echo.py --level error; exit 1; }}
            fi

            """

    localrules: merge_classification
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
            Log={Wkdir}/log/iter-0/step3-classify/classify-trimmed-merge.log
            ### merge clf from all groups
            python {Scriptdir}/merge-clf-from-groups.py {output.clf} {params.clf_fs_str} 2> $Log || {{ echo "See error details in $Log" | python {Scriptdir}/echo.py --level error; exit 1; }}
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
        shell:
            """
            cp iter-0/viral-combined-proba.tsv final-viral-score.tsv
            cp iter-0/viral-combined.fa final-viral-combined.fa
            cp iter-0/viral-fullseq.tsv final-viral-boundary.tsv
            grep -v '^seqname' iter-0/viral-partseq.tsv >> final-viral-boundary.tsv || : 
            N_lt2gene=$(grep -c '^>' iter-0/viral-lt2gene-w-hallmark.fa || :)
            N_lytic=$(grep -c '^>' iter-0/viral-fullseq-trim.fa || :)
            N_lysogenic=$(grep -c '^>' iter-0/viral-partseq.fa || :)
            echo -e "
            ====> VirSorter run (provirus mode) finished.
            # of full     seqs as viral:\t$N_lytic
            # of partial  seqs as viral:\t$N_lysogenic
            # of short (<2 genes) seqs as viral:\t$N_lt2gene

            Useful output files:
            final-viral-score.tsv       ==> score table
            final-viral-combined.fa     ==> all viral seqs
            final-viral-boundary.tsv    ==> table with boudary info
            
            Suffix are added to seqname in final-viral-combined.fa:
            full seqs as viral:               ||full 
            partial seqs as viral:            ||partial
            short (<2 genes) seqs as viral:   ||lt2gene
            <====
            " | python {Scriptdir}/echo.py
            """
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
        shell:
            """
            python {Scriptdir}/add-suffix-seqname-keep-desc.py iter-0/viral-fullseq.fa ||full > {output.combined}
            python {Scriptdir}/add-suffix-seqname-keep-desc.py {input.lt2gene} ||lt2gene >> {output.combined}
            cp iter-0/all-fullseq-proba.tsv {output.proba}
            N_viral_fullseq=$(grep -c '^>' iter-0/viral-fullseq.fa || :)
            N_viral_lt2gene=$(grep -c '^>' iter-0/viral-lt2gene-w-hallmark.fa || :)
            echo -e "
            ====> VirSorter run (non-provirus mode) finished.
            # of viral contigs:\t$N_viral_fullseq
            # of short viral contigs (<2 genes):\t$N_viral_lt2gene

            Useful output files:
            final-viral-score.tsv
            final-viral-combined.fa
            <====
            " | python {Scriptdir}/echo.py
            """
