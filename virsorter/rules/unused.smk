
#rule make_diamond_db_on_faa:
#    input: 'iter-0/all.pdg.faa'
#    output: 'iter-0/all.pdg.faa.dmnd'
#    shell:
#        """
#        # output add .dmnd suffix to --db
#        diamond makedb --in {input} --db {input}
#        """
#
#rule make_diamond_db_on_faa_by_group:
#    input: 
#        dmnd='iter-0/all.pdg.faa.dmnd',
#        faa='iter-0/{group}/all.pdg.faa',
#    output: 'iter-0/{group}/all.pdg.faa.dmnd'
#    shell:
#        """
#        Rbs_pdg_db={Dbdir}/group/{wildcards.group}/rbs-prodigal-train.db
#        if [ -s $Rbs_pdg_db ]; then
#            diamond makedb --in {input.faa} --db {input.faa}
#        else
#            (cd iter-0/{wildcards.group} && ln -sf ../all.pdg.faa.dmnd)
#        fi
#        """
#
####
## loop version of hmm update and classify iteration 
####
#rule update_hmm_and_classify_loop:
#    input:
#        'iter-0/viral-contig-proba.tsv',
#        'iter-0/viral-contig.fa',
#        expand('iter-0/{group}/all.pdg.faa.dmnd', group=Groups)
#    output:
#        clf='viral-contig-proba-final.tsv',
#        contig='viral-contig-final.fa'
#    threads: General_threads
#    #conda: 'env/vfam.yml'
#    shell:
#        """
#        if [ {Iter_max} -eq 0 ]; then
#            cp iter-0/viral-contig.fa viral-contig-final.fa
#            cp iter-0/viral-contig-proba.tsv viral-contig-proba-final.tsv
#            echo "*** Max # of iteration is set to 1; exiting..."
#            exit 0
#        fi
#
#        if grep -q "^>" iter-0/viral-contig.fa; then
#            Cnt_prev=$(grep "^>" -c iter-0/viral-contig.fa)
#        else
#            Cnt_prev=0
#        fi
#
#        if [ $Cnt_prev -eq 0 ]; then
#            echo "*** No viral seqs are detected at first iteration..."
#            cp iter-0/viral-contig.fa viral-contig-final.fa 
#            cp iter-0/viral-contig-proba.tsv viral-contig-proba-final.tsv
#            exit 0
#        fi
#
#        for group in {Groups_str}; do
#            cat iter-0/$group/all.pdg.hmm.tax | cut -f1 | sort | uniq > iter-0/$group/prot-hits.list
#            # get protein that does not have hits in current db
#            python {Scriptdir}/pick-seq-not-in-list.py iter-0/$group/viral-contig-prot.faa iter-0/$group/prot-hits.list > iter-0/$group/unannotated.faa
#            python {Scriptdir}/get-list-from-seqfile.py iter-0/$group/unannotated.faa > iter-0/$group/unannotated.list
#        done
#
#        ####
#        # exit on no more new viral protein
#        ####
#        cat iter-0/*/unannotated.faa > iter-0/unannotated.faa
#
#        if grep -q "^>" iter-0/unannotated.faa; then
#            Unanno_cnt=$(grep -c "^>" iter-0/unannotated.faa)
#        else
#            Unanno_cnt=0
#        fi
#
#        #if [ $Unanno_cnt -eq 0 -o {Iter_max} -eq 0 ]; then
#        #    cp iter-0/viral-contig.fa viral-contig-final.fa
#        #    cp iter-0/viral-contig-proba.tsv viral-contig-proba-final.tsv
#        #    echo "*** No more viral proteins found after iteration 0; exiting..."
#        #    exit 0
#        #fi
#
#        for ((i=1; i<={Iter_max}; i++)); do
#
#            echo "=================================="
#            echo " Starting iter-${{i}}"
#            echo "=================================="
#            mkdir -p iter-${{i}}
#            ( #################################### start of subshell
#            cd iter-${{i}}
#
#            # vfam to build hmm, concatenated hmm as clustered.hmm,
#            #  a. clusterd seqs as clustered_inProfiles.fasta 
#            #  b. deduped seqs as fastaFiles/unanno_prev_noDupes_minL1_c100.fa
#            #  c. difference between a. and b. are singletons
#            # -m 1 --> minL1
#            # -p 1 --> c100 
#            mkdir -p vfam
#            (cd vfam && ln -sf ../../iter-$((i-1))/unannotated.faa unanno_prev.faa && {Scriptdir}/vfam.py -f unanno_prev.faa -n 2 -m 1 -p 1 -a {threads} -o clustered)
#            python {Scriptdir}/diff-seq.py vfam/fastaFiles/unanno_prev_noDupes_minL1_c100.faa vfam/clustered_inProfiles.fasta > vfam/unclustered.faa
#            # hmmsearch
#            hmmsearch -T {Hmmsearch_score_min} --cpu {threads} --tblout new-viral-hits-hmmsearch.tbl vfam/clustered.hmm ../iter-0/all.pdg.faa
#            # diamond
#            diamond blastp --min-score {Diamond_score_min} --threads {threads} --max-target-seqs 500 --block-size 1  --db ../iter-0/all.pdg.faa --query vfam/unclustered.faa --more-sensitive --out new-viral-hits-diamond.tbl --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen
#            # merge output
#            grep -v "^#" new-viral-hits-hmmsearch.tbl | awk '{{print $1}}' > new-viral-hits-hmmsearch.list
#            grep -v "^#" new-viral-hits-diamond.tbl | awk '{{print $2}}' > new-viral-hits-diamond.list
#
#            # choose to do loop instead of parallel 
#            #   since only NCLDV has different Rbs_pdg_db and thus all.pdg.faa
#            ##################################
#            Groups_str="{Groups_str}"
#            for group in $Groups_str; do
#                Rbs_pdg_db={Dbdir}/group/$group/rbs-prodigal-train.db
#                Hallmark_list_f={Dbdir}/group/$group/hallmark-gene.list
#                mkdir -p $group
#                if [ -s $Rbs_pdg_db ]; then
#                    #2. hmmsearch
#                    hmmsearch -T {Hmmsearch_score_min} -o /dev/null --cpu {threads} --tblout $group/new-viral-hits-hmmsearch.tbl vfam/clustered.hmm ../iter-0/$group/all.pdg.faa
#                    #3. diamond
#                    diamond blastp --min-score {Diamond_score_min} --threads {threads} --max-target-seqs 500 --block-size 1  --db ../iter-0/$group/all.pdg.faa --query vfam/unclustered.faa --more-sensitive --out $group/new-viral-hits-diamond.tbl --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen
#                    #4. merge output
#                    grep -v "^#" $group/new-viral-hits-hmmsearch.tbl | awk '{{print $1}}' > $group/new-viral-hits-hmmsearch.list
#                    grep -v "^#" $group/new-viral-hits-diamond.tbl | awk '{{print $2}}' > $group/new-viral-hits-diamond.list
#
#                    cat ../iter-$((i-1))/$group/unannotated.list $group/new-viral-hits-hmmsearch.list $group/new-viral-hits-diamond.list | sort | uniq > $group/new-viral-hits.list
#                    #5. update hmm.tax; add seq_name, tax(vir), -, -
#                    python {Scriptdir}/update-taxfile.py $group/new-viral-hits.list ../iter-$((i-1))/$group/all.pdg.hmm.tax > $group/all.pdg.hmm.tax
#
#                    if [ -s $Hallmark_list_f ]; then
#                        python {Scriptdir}/add-unaligned-to-hmm-featrues-w-hallmark.py ../iter-0/$group/all.pdg.faa $group/all.pdg.hmm.tax $Hallmark_list_f > $group/all.pdg.hmm.ftr
#                    else
#                        python {Scriptdir}/add-unaligned-to-hmm-featrues.py ../iter-0/$group/all.pdg.faa $group/all.pdg.hmm.tax > $group/all.pdg.hmm.ftr
#                    fi
#                else
#                    (cd $group && ln -sf ../new-viral-hits-hmmsearch.list && ln -sf ../new-viral-hits-diamond.list)
#                    cat ../iter-$((i-1))/$group/unannotated.list $group/new-viral-hits-hmmsearch.list $group/new-viral-hits-diamond.list | sort | uniq > $group/new-viral-hits.list
#                    #5. update hmm.tax; add seq_name, tax(vir), -, -
#                    python {Scriptdir}/update-taxfile.py $group/new-viral-hits.list ../iter-$((i-1))/$group/all.pdg.hmm.tax > $group/all.pdg.hmm.tax
#
#                    if [ -s $Hallmark_list_f ]; then
#                        python {Scriptdir}/add-unaligned-to-hmm-featrues-w-hallmark.py ../iter-0/all.pdg.faa $group/all.pdg.hmm.tax $Hallmark_list_f > $group/all.pdg.hmm.ftr
#                    else
#                        python {Scriptdir}/add-unaligned-to-hmm-featrues.py ../iter-0/all.pdg.faa $group/all.pdg.hmm.tax > $group/all.pdg.hmm.ftr
#                    fi
#                fi
#
#                python {Scriptdir}/merge-hmm-gff-features.py ../iter-0/$group/all.pdg.gff.ftr $group/all.pdg.hmm.ftr > $group/all.pdg.ftr
#
#                # classify
#                if [ -s $Hallmark_list_f ]; then
#                    python {Scriptdir}/classify.py $group/all.pdg.ftr {Dbdir}/group/$group/scaler-w-hm.scl {Dbdir}/group/$group/model-w-hm.clf $group/all.pdg.clf
#                else
#                    python {Scriptdir}/classify.py $group/all.pdg.ftr {Dbdir}/group/$group/scaler.scl {Dbdir}/group/$group/model.clf $group/all.pdg.clf
#                fi
#
#            # finish group loop
#            done #######################################
#
#            # merge clf from all groups
#            python {Scriptdir}/merge-clf-from-groups.py viral-contig-proba.tsv */all.pdg.clf
#            python {Scriptdir}/pick-viral-contig-from-clf.py {Proba_cutoff} viral-contig-proba.tsv {Seqfile} > viral-contig.fa
#
#            Arr=($Groups_str)
#            Arr2=()
#            Arr3=()
#
#            for i in ${{Arr[@]}}; do
#                Arr2+=(../iter-0/$i/all.pdg.faa)
#                Arr3+=($i/viral-contig-prot.faa)
#            done
#
#            Groups_str_comma_sep=$(echo $Groups_str | tr " " ",")
#            Faa_fs_str=$(echo ${{Arr2[@]}} | tr " " ",")
#            Vir_faa_fs_str=$(echo ${{Arr3[@]}} | tr " " ",")
#            python {Scriptdir}/pick-viral-faa-from-clf-by-group.py {Proba_cutoff} viral-contig-proba.tsv "$Groups_str_comma_sep" "$Faa_fs_str" "$Vir_faa_fs_str"
#            cat $(echo $Vir_faa_fs_str | tr "," " ") > viral-contig-prot-combined.faa
#
#            ##################
#            for group in $Groups_str; do
#                cat $group/all.pdg.hmm.tax | cut -f1 | sort | uniq > $group/prot-hits.list
#                # get protein that does not have hits in current db
#                python {Scriptdir}/pick-seq-not-in-list.py $group/viral-contig-prot.faa $group/prot-hits.list > $group/unannotated.faa
#                python {Scriptdir}/get-list-from-seqfile.py $group/unannotated.faa > $group/unannotated.list
#            done
#            #################
#            ) ############################# end of subshell
#
#            ####
#            # exit on no more new viral protein
#            ####
#            cat iter-${{i}}/*/unannotated.faa > iter-${{i}}/unannotated.faa
#            if grep -q "^>" iter-${{i}}/unannotated.faa; then
#                Unanno_cnt=$(grep -c "^>" iter-${{i}}/unannotated.faa)
#            else
#                Unanno_cnt=0
#            fi
#
#            if [ $Unanno_cnt -eq 0 ]; then
#                echo "*** No more viral proteins found after iteration ${{i}}; exiting..."
#                cp iter-${{i}}/viral-contig.fa viral-contig-final.fa
#                cp iter-${{i}}/viral-contig-proba.tsv viral-contig-proba-final.tsv
#                exit 0
#            fi
#
#            ###
#            # exit on max iteration limit
#            ###
#            if [ $i -eq {Iter_max} ]; then
#                cp iter-${{i}}/viral-contig.fa viral-contig-final.fa
#                cp iter-${{i}}/viral-contig-proba.tsv viral-contig-proba-final.tsv
#                echo "*** Reaching maximum iteration # of $Iter_max; exiting..."
#                exit 0
#            fi
#
#            # end of iter i loop
#            echo "============================================="
#            echo "Finishing iter-${{i}}"
#            echo "============================================="
#
#        done #####################
#        """
#
#####################################################################
# recursive rule version of hmm update and classify

# Start iterative viral gene call and viral contig indentification
# phage protein catalog 
# merge  annotation from hmm and blast (simon's step 2)
# detect (simon's step 3)
# summary (simon's step 4)
# generate output (simon's step 5)
###########################################################

#def prep_viral_db_recursive_input(wc):
#    n = int(wc.n)
#    if n > Iteration_max:
#        raise ValueError('Iteration number exceeded "Iteration_max": {}'.format(Iteration_max))
#    elif n == 1:
#        return '{}/{}/{}-PC.hmm'.format(Vdbdir, wc.group, Db_label)
#    else:
#        # multiple groups in Vdbdir folder
#        return 'iteration-{}/global-summary.tsv'.format(n-1)
#
#
##---------------
#### possible way of not using loop
#rule prep_viral_db:
#    input: prep_viral_db_recursive_input
#    output: 
#        hmm='iteration-{n}/{group}/db/PC.hmm',
#        unclustered='iteration-{n}/{group}/db/unclustered.faa',
#    wildcard_constraints:
#        n = '[0-9]+'
#    threads: Mmseqs_threads
#    params:
#        m = lambda wildcards: int(wildcards.n) - 1,
#        prefix = lambda wildcards: 'iteration-{}/{}'.format(wildcards.n, wildcards.group)
#    shell:
#        """
#        # no seq in new_viral_protein.faa from previous iteration
#        #   just touch output files to skip rest of iterations
#        if [ -f 'Done' ]; then
#            touch {output.hmm} {output.unclustered}
#            exit 0
#        fi
#
#        mkdir -p {params.prefix}/db
#
#        if [ {wildcards.n} -eq 1]; then
#            if [ ! -f {Vdbdir}/{group}/{Db_label}-PC.hmm ]; then
#                if [ -f {Vdbdir}/{group}/customized.faa ]; then
#                    Absdir=$(cd {Vdbdir}/{group} && pwd)
#                    ln -s $Absdir/customized.faa \
#                        {params.prefix}/db/pooled.faa
#                    mmseqs createdb \
#                        {params.prefix}/db/pooled.faa \
#                        {params.prefix}/db/pooled.faa.mmdb
#                    mmseqs cluster \
#                        --threads {threads} \
#                        {params.prefix}/db/pooled.faa.mmdb \
#                        {params.prefix}/db/pooled.faa.clu
#                    ln -s {params.prefix}/db/pooled.faa.mmdb \
#                        {params.prefix}/db/pooled.faa.mmdb.update
#                    ####
#                    # update the PC.hmm and unclustered.faa
#                    ###
#                else
#                    echo -e "*** Error: No {Db_lable}.faa or costomized.faa found in {Vdbdir}/{group}/db; at least one is needed to proceed\n"
#                    exit 1
#                fi
#            else
#                if [ -f {Vdbdir}/{group}/customized.faa ]; then
#                    cat {Vdbdir}/{group}/unclustered.faa \
#                        {Vdbdir}/{group}/customized.faa \
#                        > {params.prefix}/db/pooled.faa
#                    mmseqs createdb \
#                        {params.prefix}/db/pooled.faa \
#                        {params.prefix}/db/pooled.faa.mmdb
#                    mmseqs cluster \
#                        --threads {threads} \
#                        {params.prefix}/db/pooled.faa.mmdb \
#                        {params.prefix}/db/pooled.faa.clu
#                    ln -s {params.prefix}/db/pooled.faa.mmdb \
#                        {params.prefix}/db/pooled.faa.mmdb.update
#                    ####
#                    # update the PC.hmm and unclustered.faa
#                    ###
#                else
#                    Absdir=$(cd {Vdbdir}/{group} && pwd)
#                    ln -s $Absdir/unclustered.faa \
#                        iteration_{n}/{group}/db/pooled.faa
#                    mmseqs createdb {params.prefix}/db/pooled.faa \
#                        {params.prefix}/db/pooled.faa.mmdb
#                    ln -s {params.prefix}/db/pooled.faa.mmdb \
#                        {params.prefix}/db/pooled.faa.mmdb.update
#                    ln -s $Absdir/{Db_label}-PC.hmm \
#                        {output.hmm}
#                    ln -s $Absdir/{Db_label}-unclustered.faa \
#                        {output.unclustered}
#                    # no customized protein seqs provided
#                    #   no need to rebuild PCs at interation 1
#                    exit 0
#                fi
#            fi
#        else
#            # after first iteration n >= 2
#            cat iteration_{m}/{group}/db/new-viral-protein.faa \
#                iteration_{m}/{group}/db/pooled.faa \
#                > {params.prefix}/db/pooled.fa
#            # run mmseq for making new protein clusters
#            mmseqs createdb \
#                {params.prefix}/db/pooled.faa \
#                {params.prefix}/db/pooled.faa.mmdb
#            # mmseqs clusterupdate <old-db> <new-db> <old-cluster> <new-db.updated> <new-cluster> <tmpdir>
#            mmseqs clusterupdate \
#                --threads {threads} \
#                iteration_{m}/{group}/db/pooled.faa.mmdb.updated \
#                {params.prefix}/db/pooled.faa.mmdb \
#                iteration_{m}/{group}/db/pooled.faa.clu \
#                {params.prefix}/db/pooled.faa.mmdb.updated \
#                {params.prefix}/db/pooled.faa.clu \
#                {params.prefix}/db/mmseqs-tmp 
#            # work around a bug on multi-threading
#            #   each thread create cluster file
#            if [ ! -f {params.prefix}/db/pooled.faa.clu ]; then
#                cat {params.prefix}/db/pooled.faa.clu.+([0-9]) \
#                    > {params.prefix}/db/pooled.faa.clu
#                rm {params.prefix}/db/pooled.faa.clu.+([0-9])
#            fi
#        fi
# 
#        ###
#        # update the PC.hmm and unclustered.faa
#        # e.g.
#        # mmseqs createseqfiledb DB clu clu_seq 
#        # mmseqs result2flat DB DB clu_seq clu_seq.fasta
#        ### 
#        mmseqs createseqfiledb \
#            --threads {threads}
#            {params.prefix}/db/pooled.faa.mmdb.updated \
#            {params.prefix}/db/pooled.faa.clu \
#            {params.prefix}/db/pooled.faa.cluseq
#        mmseqs result2flat \
#            {params.prefix}/db/pooled.faa.mmdb.updated \
#            {params.prefix}/db/pooled.faa.mmdb.updated \
#            {params.prefix}/db/pooled.faa.cluseq \
#            {params.prefix}/db/pooled.faa.cluseq.mmfa
#        python scripts/parse-mmfa.py \
#            {params.prefix}/db/pooled.faa.cluseq.mmfa \
#            {params.prefix}/db/clustered.outdir \
#            {params.prefix}/db/unclustered.faa
#
#        # realign with clustalo and make hmm
#        #mmseqs apply --threads {threads} \
#        #    {params.prefix}/db/pooled.faa.cluseq \
#        #    {params.prefix}/db/pooled.faa.cluseq.msa \
#        #    -- \
#        #    clustalo --threads=1 -i - -o -
#        #scripts/split-mmmsa-by-cluster.py \
#        #    {params.prefix}/db/pooled.faa.cluseq.msa \
#        #    {params.prefix}/db/cluster.dir
#
#        ls {params.prefix}/db/clustered.outdir/*.fa \
#            | xargs -n 1 -P {threads} -I File bash -c \
#            'Bname=$(basename $File .fa); clustalo --threads=1 -i $File -o {params.prefix}/db/clustered.outdir/$Banme.afa; hmmbuild --cpu 1 {params.prefix}/db/clustered.outdir/$Bname.hmm {params.prefix}/db/clustered.outdir/$Bname.afa'
#        cat {params.prefix}/db/clustered.outdir/*.hmm \
#            > {params.prefix}/db/PC.hmm
#
#        """
#
#
#def search_viral_protein_recursive_input(wc):
#    n = int(wc.n)
#    if n > Iteration_max:
#        raise ValueError('Iteration number exceeded "Iteration_max": {}'.format(Iteration_max))
#    else:
#        return 'iteration-{}/{}/db/PC.hmm'.format(n, wc.group), 'iteration-{}/{}/db/unclustered.faa'.format(n, wc.group)
#
#
#rule search_viral_protein:
#    input: search_viral_protein_recursive_input
#    output:
#        hmmout='iteration-{n}/{group}/hmm.tblout',
#        mmseqout='iteration-{n}/{group}/mmseq.tblout',
#    shell:
#        """
#        # no seq in new_viral_protein.faa from previous iteration
#        #   just touch output files to skip rest of iterations
#        if [ -f 'Done' ]; then
#            touch {output.hmmout} {output.mmseqout}
#            exit 0
#        fi
#        """
#
#def merge_tbl_recursive_input(wc):
#    n = int(wc.n)
#    if n > Iteration_max:
#        raise ValueError('Iteration number exceeded "Iteration_max": {}'.format(Iteration_max))
#    else:
#        return 'iteration-{}/{}/hmm.tblout'.format(n, wc.group), 'iteration-{}/{}/mmseq.tblout'.format(n, wc.group)
#
#rule merge_tbl:
#    input: merge_tbl_recursive_input
#    output: 'iteration-{n}/{group}/merged.tblout'
#    shell:
#        """
#        # no seq in new_viral_protein.faa from previous iteration
#        #   just touch output files to skip rest of iterations
#        if [ -f 'Done' ]; then
#            touch {output}
#            exit 0
#        fi
#        touch {output}
#        """
#
#def detect_virus_recursive_input(wc):
#    n = int(wc.n)
#    if n > Iteration_max:
#        raise ValueError('Iteration number exceeded "Iteration_max": {}'.format(Iteration_max))
#    else:
#        return 'iteration-{}/{}/merged.tblout'.format(n, wc.group)
#
#rule detect_virus:
#    input: detect_virus_recursive_input
#    output:
#        'iteration-{n}/{group}/db/new-viral-protein.faa'
#    shell:
#        """ 
#        # update "all.pdf.hmm.features" here and then rerun the classifer
#        # no seq in new_viral_protein.faa from previous iteration
#        #   just touch output files to skip rest of iterations
#        if [ -f 'Done' ]; then
#            touch {output}
#            exit 0
#        fi
#        touch {output}
#        """
#
#def summarize_across_groups_recursive_input(wc):
#    n = int(wc.n)
#    if n > Iteration_max:
#        raise ValueError('Iteration number exceeded "Iteration_max": {}'.format(Iteration_max))
#    else:
#        return ['iteration-{}/{}/db/new-viral-protein.faa'.format(n, group) for group in Groups]
#
#rule summarize_across_groups:
#    input: summarize_across_groups_recursive_input
#    output: 'iteration-{n}/global-summary.tsv'
#    shell:
#        """
#        # no seq in new_viral_protein.faa from previous iteration
#        #   just touch output files to skip rest of iterations
#        if [ -f 'Done' ]; then
#            touch {output}
#            exit 0
#        fi
#        touch {output}
#        """
