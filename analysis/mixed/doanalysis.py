import os
import time
import logging
import json
import shlex
import pathlib

import argh

logging.basicConfig(level=logging.INFO)


def run(cmd):
    logging.info(f"command: {cmd}")
    os.system(cmd)


def make_mixed_samples(
    input_fasta, outprefix, number_of_mixed, number_of_random_ns, number_of_primer_gaps
):
    run(
        f"python3 mixed_utils.py {input_fasta} --number-of-mixed {number_of_mixed} -r 1 -o {outprefix} --number-of-random-ns {number_of_random_ns} --number-of-primer-gaps {number_of_primer_gaps}"
    )


def go(
    input_fasta="sim0.fasta",
    k_start=0,
    k_end=10,
    number_of_mixed=50,
    number_of_random_ns=50,
    number_of_primer_gaps=2,
):
    run_prefix = f"{input_fasta}_mixed{number_of_mixed}_randomNs{number_of_random_ns}_primergaps{number_of_primer_gaps}"
    pathlib.Path(run_prefix).mkdir()
    print(f"{run_prefix=}")

    ## create samples
    make_mixed_samples(
        input_fasta,
        run_prefix,
        number_of_mixed,
        number_of_random_ns,
        number_of_primer_gaps,
    )

    for k in range(k_start, k_end + 1):
        loop_prefix = f"{run_prefix}/k{k}"
        print("{loop_prefix=}")
        pathlib.Path(loop_prefix).mkdir()

        ## load them into catwalk
        run(
            "nohup ../../cw_server --instance-name=test_sim --reference-filepath=../../reference/nc_045512.fasta --mask-filepath=../../reference/covid-exclude.txt &"
        )
        time.sleep(2)
        post_data = json.dumps({"filepath": f"{run_prefix}/mixed.fasta"})
        run(
            f"curl -X POST http://localhost:5000/add_samples_from_mfsl -H 'Content-Type: application/json' -d {shlex.quote(post_data)}"
        )
        time.sleep(2)
        ## do analysis
        run(f"python3 mixed_analysis.py make-csv-and-graphs {k} {loop_prefix}")

        ## remove samples with fewer than 2 neighbours
        run(
            f"cat {loop_prefix}/mixanalysis.csv | grep -v ,, > {loop_prefix}/mixanalysis-cleaned.csv"
        )
        ## make roc curve images and csv
        run(f"python3 mixed_analysis.py make-roc-curve-graphs {loop_prefix}")
        #
        ## stop catwalk
        run("pkill cw_server")


if __name__ == "__main__":
    argh.dispatch_command(go)
