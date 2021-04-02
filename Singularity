Bootstrap: docker
From: ubuntu:18.04

%environment
    export PATH=/miniconda3/bin:$PATH

%runscript
    exec virsorter "$@"

%post
    apt-get update && apt-get install -y automake build-essential bzip2 wget git unzip

    export PATH=/miniconda3/bin:$PATH

    # Install miniconda to save dependency nightmares
    wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
    bash Miniconda3-latest-Linux-x86_64.sh -b -f -p /miniconda3/

    . /miniconda3/etc/profile.d/conda.sh  # Only activates conda, but don't need to "activate base"
    conda install -y conda-build

    # Use the same as development version installation
    conda install -y -c conda-forge -c bioconda "python>=3.6" scikit-learn=0.22.1 imbalanced-learn pandas seaborn hmmer==3.3 prodigal screed ruamel.yaml "snakemake>=5.18,<=5.26" click mamba
    git clone https://github.com/jiarong/VirSorter2.git
    cd VirSorter2
    pip install .

    # generate template-config.yaml;  db_dir ONLY for cyverse  app
    virsorter config --init-source --db-dir /work/projects/Cyverse/iVirus/VS2/db
 
    # TACC's Stampede compliant,for iVirus/CyVerse
    mkdir /home1 && mkdir /scratch && mkdir /work && mkdir /trigger-rebuild
