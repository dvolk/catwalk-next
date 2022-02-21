# Mixed sample analysis

## Step 1

Creates mixed samples based on existing samples and parameters, then inserts into
catwalk, queries neighbourhood and creates graphs and saves statistics in csv

eg.:

    python3 doanalysis.py --number-of-random-ns 0 --number-of-primer-gaps 0

in this example, the output files would be written to the directory `sim0.fasta_mixed50_randomNs0_primergaps0/`

parameters:

    usage: doanalysis.py [-h] [-i INPUT_FASTA] [--k-start K_START] [--k-end K_END] [--number-of-mixed NUMBER_OF_MIXED] [--number-of-random-ns NUMBER_OF_RANDOM_NS] [--number-of-primer-gaps NUMBER_OF_PRIMER_GAPS]

    optional arguments:
      -h, --help            show this help message and exit
      -i INPUT_FASTA, --input-fasta INPUT_FASTA
                            'sim0.fasta'
      --k-start K_START     0
      --k-end K_END         10
      --number-of-mixed NUMBER_OF_MIXED
                            50
      --number-of-random-ns NUMBER_OF_RANDOM_NS
                            50
      --number-of-primer-gaps NUMBER_OF_PRIMER_GAPS
                            2

## Step 2

Creates PNG graphs of variation of AUC based on cut-off from statistics in csv

eg.:

    python3 doanalysis2.py sim0.fasta_mixed50_randomNs0_primergaps0/
