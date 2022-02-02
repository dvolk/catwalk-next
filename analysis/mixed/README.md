# Mixed sample analysis

## Creating mixed samples

### Copy and unarchive `benchmark/sim/sim0.fasta.gz`

    cp ../../benchmark/sim/sim0.fasta.gz .
    tar xf sim0.fasta.gz

### Generate combined unmixed and mixed fasta file

eg.:

    python3 mixed_utils.py -r 1 -o 1050randomNs sim0.fasta --number-of-random-ns 50 --number-of-primer-gaps 2

this will create the files `1050randomNs.fasta` and `1050randomNs.csv`

## mixed sample analysis process description

### Start catwalk

    ./cw_server --instance-name=test_sim --reference-filepath=reference/nc_045512.fasta --mask-filepath=reference/covid-exclude.txt

### Load samples into catwalk

    python3
    import requests
    requests.post("http://localhost:5000/add_samples_from_mfsl", json={"filepath": "analysis/mixed/1050randomNs.fasta"})

### Generate graphs and mixed analysis csv

    python3 mixed_analysis.py make-csv-and-graphs 3 1050randomNs

### Remove samples with fewer than 3 neighbours

    cat 1050randomNs-mixanalysis.csv | grep -v ,, > 1050randomNs-mixanalysis-cleaned.csv

### Generate ROC curve graphs and roc csv

    python3 mixed_analysis.py make-roc-curve-graphs 1050randomNs-mixanalysis-cleaned.csv 1050randomNs

The output files are:

    1050randomNs.csv
    1050randomNs.fasta
    1050randomNs-mixanalysis-cleaned.csv
    1050randomNs-mixanalysis.csv
    1050randomNs-roc.csv
    1050randomNs-igraph_pngs/*
    1050randomNs-roc_curves/*
