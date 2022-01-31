# mixed sample analysis process description

## Unarchive samples

    tar xf mix2.fa.gz

## Start catwalk

    ./cw_server --instance-name=test_sim --reference-filepath=reference/nc_045512.fasta --mask-filepath=reference/covid-exclude.txt

## Load mix2.fa samples into catwalk

    python3
    import requests
    requests.post("http://localhost:5000/add_samples_from_mfsl", json={"filepath": "analysis/mixed/mix2.fa"})

## Generate graphs and mixed analysis csv

    python3 mixed_analysis.py make-csv-and-graphs 3

## Remove samples with fewer than 3 neighbours

    cat mixanalysis.csv | grep -v ,, > mixanalysisclean.csv

## Generate ROC curve graphs and roc csv

    python3 mixed_analysis.py make-roc-curve-graphs mixanalysisclean.csv
