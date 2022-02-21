"""Generate mixed samples from existing samples by picking
pairs that are close together in time, and making a new sample
where the different bases are N.

This script can only be used with the sequences with the header
format ">sampleid:yyyy-mm-dd:etc"

usage: mixed_utils.py [-h] [-n NUMBER_OF_MIXED] [-d DAYS] [-o OUTPREFIX] mfsl-filepath

parameters:

[-n NUMBER_OF_MIXED] - number of mixed samples to create
[-d DAYS] - only accept sample pairs for mixing if they are < DAYS
apart based on the sequence header
[-o OUTPREFIX] - the file prefix for the output fasta and csv files

mixed samples are written to outprefix.fasta and a csv with sample data
is written to outprefix.csv"""

from dataclasses import dataclass
import datetime
import random
import pathlib

import argh


@dataclass
class Sample:
    header: str
    sequence: str
    epochtime: int
    name: str
    mix_debug: str
    count_ns_mixture: int
    count_ns_random: int
    count_ns_primer_dropout: int
    mix_parent1: str
    mix_parent2: str


def format_sample(sample):
    return sample.header + "\n" + sample.sequence + "\n"


outcsv_header = (
    ",".join(
        [
            "header",
            "name",
            "epochtime",
            "count_ns_mixture",
            "count_ns_random",
            "count_ns_primer_dropout",
            "mix_parent1",
            "mix_parent2",
        ]
    )
    + "\n"
)


def format_sample_csvline(sample):
    row = [
        sample.header,
        sample.name,
        str(int(sample.epochtime)),
        str(sample.count_ns_mixture),
        str(sample.count_ns_random),
        str(sample.count_ns_primer_dropout),
        sample.mix_parent1,
        sample.mix_parent2,
    ]
    return ",".join(row) + "\n"


def parse_header(sample):
    try:
        elems = sample.header[1:].split(":")
        d = datetime.datetime.strptime(elems[1], "%d-%m-%Y")
        sample.epochtime = int(d.strftime("%s"))
        sample.name = elems[0]
    except Exception:
        pass


def insert_random_ns(sequence, number_of_ns):
    # create list of eligable positions
    # this is inefficient on large sequences
    candidates = [i for i, b in enumerate(sequence) if b in "ACGT"]
    # choose positions from candidates
    number_of_ns = min(number_of_ns, len(candidates))
    mask_idx = random.sample(candidates, k=number_of_ns)
    # mask sequence
    newseq = list(sequence)
    for i in mask_idx:
        newseq[i] = "N"
    return "".join(newseq)


def sample_insert_random_ns(sample, number_of_ns):
    sample.sequence = insert_random_ns(sample.sequence, number_of_ns)
    sample.count_ns_random = number_of_ns


primer_gaps = list()


def mask_primer_gaps(sequence, number_of_primer_gaps):
    count = 0
    ranges_to_mask = random.sample(primer_gaps, k=number_of_primer_gaps)
    newseq = list(sequence)
    for range_to_mask in ranges_to_mask:
        for x in range_to_mask:
            newseq[x] = "N"
            count += 1
    return "".join(newseq), count


def sample_mask_primer_gaps(sample, number_of_primer_gaps):
    sample.sequence, n = mask_primer_gaps(sample.sequence, number_of_primer_gaps)
    sample.count_ns_primer_dropout = n


def parse_mfsl(mfsl_filepath):
    mfsl = open(mfsl_filepath).readlines()
    samples = list()
    for linenum, _ in enumerate(mfsl):
        if linenum % 2 == 1:
            s = Sample(
                sequence=mfsl[linenum].strip().upper(),
                header=mfsl[linenum - 1].strip(),
                name="",
                mix_debug="",
                epochtime=0,
                count_ns_mixture=0,
                count_ns_random=0,
                count_ns_primer_dropout=0,
                mix_parent1="",
                mix_parent2="",
            )
            parse_header(s)
            samples.append(s)
    return samples


def generate_all_sample_pairs(samples, days):
    days_in_seconds = days * 24 * 60 * 60
    pairs = list()
    for sam1 in samples:
        for sam2 in samples:
            if abs(sam1.epochtime - sam2.epochtime) < days_in_seconds:
                if sam1.sequence != sam2.sequence:
                    pairs.append((sam1, sam2))
    return pairs


def mix_pair(pair, number_of_random_ns, number_of_primer_gaps):
    sam1, sam2 = pair
    new_seq = list()
    count_mixture_ns = 0
    for i, (base1, base2) in enumerate(zip(sam1.sequence, sam2.sequence)):
        if base1 in "ACGT" and base2 in "ACGT" and base1 != base2:
            new_seq.append("N")
            count_mixture_ns += 1
        else:
            if i % 2 == 0:
                new_seq.append(base1)
            else:
                new_seq.append(base2)
    s = Sample(
        header=f">{sam1.name}+{sam2.name}",
        sequence="".join(new_seq),
        epochtime=(sam1.epochtime + sam2.epochtime) / 2,
        name=f"{sam1.name}+{sam2.name}",
        mix_debug=f"{sam1.sequence}+{sam2.sequence}",
        count_ns_mixture=count_mixture_ns,
        # these two counts are approximate?
        count_ns_random=(sam1.count_ns_random + sam2.count_ns_random) // 2,
        count_ns_primer_dropout=(
            sam1.count_ns_primer_dropout + sam2.count_ns_primer_dropout
        )
        // 2,
        mix_parent1=sam1.name,
        mix_parent2=sam2.name,
    )

    return s


def parse_covid_bed_file():
    """Populate global variable list primer_gaps with ranges from
    the covid bed file that aren't covered by more than 1 primer.

    Assumptions are made about the base range overlaps."""

    ls = open("SARS-CoV-2.insert.bed").readlines()
    rs = [x.split("\t")[1:3] for x in ls]
    global primer_gaps
    for i, _ in enumerate(rs):
        if i == 0:
            begin, end = rs[0][0], rs[1][0]
        elif i == len(rs) - 1:
            begin, end = rs[i - 1][1], rs[i][1]
        else:
            begin, end = rs[i - 1][1], rs[i + 1][0]
        primer_gaps.append(range(int(begin), int(end)))


def print_mixed_samples(
    mfsl_filepath,
    number_of_mixed=10,
    days=7.0,
    number_of_random_ns=0,
    number_of_primer_gaps=0,
    random_seed=None,
    outprefix="mix",
):
    """Add mixed samples generated by taking two random samples
    and turning their differences into ns.

    Output samples to outprefix.fasta and corresponding csv to
    outprefix.csv"""

    if number_of_primer_gaps:
        parse_covid_bed_file()

    if random_seed:
        random.seed(int(random_seed))

    all_samples = parse_mfsl(mfsl_filepath)
    all_ok_sample_pairs = generate_all_sample_pairs(all_samples, days)
    picked_pairs = random.sample(all_ok_sample_pairs, k=number_of_mixed)

    pathlib.Path(outprefix).mkdir(exist_ok=True)

    outsamples = open(f"{outprefix}/mixed.fasta", "w")
    outcsv = open(f"{outprefix}/mixed.csv", "w")

    outcsv.write(outcsv_header)

    for sam in all_samples:
        sample_insert_random_ns(sam, number_of_random_ns)
        sample_mask_primer_gaps(sam, number_of_primer_gaps)
        outsamples.write(format_sample(sam))
        outcsv.write(format_sample_csvline(sam))

    for picked_pair in picked_pairs:
        mix = mix_pair(picked_pair, number_of_random_ns, number_of_primer_gaps)
        outsamples.write(format_sample(mix))
        outcsv.write(format_sample_csvline(mix))


if __name__ == "__main__":
    argh.dispatch_command(print_mixed_samples)
