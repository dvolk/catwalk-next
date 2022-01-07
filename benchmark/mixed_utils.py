"""Generate mixed samples from existing samples by picking
pairs that are close together in time, and making a new sample
where the different bases are N."""

import argh
from dataclasses import dataclass
import datetime
import random


@dataclass
class Sample:
    header: str
    sequence: str
    epochtime: int
    name: str


def parse_header(sample):
    try:
        elems = sample.header[1:].split(":")
        d = datetime.datetime.strptime(elems[1], "%d-%m-%Y")
        sample.epochtime = int(d.strftime("%s"))
        sample.name = elems[0]
    except Exception:
        pass


def parse_mfsl(mfsl_filepath):
    mfsl = open(mfsl_filepath).readlines()
    samples = list()
    for linenum, _ in enumerate(mfsl):
        if linenum % 2 == 1:
            s = Sample(
                sequence=mfsl[linenum].strip().upper(),
                header=mfsl[linenum - 1].strip(),
                name="",
                epochtime=0,
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


def mix_pair(pair):
    sam1, sam2 = pair
    new_seq = list()
    for i, (base1, base2) in enumerate(zip(sam1.sequence, sam2.sequence)):
        if base1 in "ACGT" and base2 in "ACGT" and base1 != base2:
            new_seq.append("N")
        else:
            if i % 2 == 0:
                new_seq.append(base1)
            else:
                new_seq.append(base2)
    s = Sample(
        header=f"{sam1.name}+{sam2.name}",
        sequence="".join(new_seq),
        epochtime=(sam1.epochtime + sam2.epochtime) / 2,
        name=f"{sam1.sequence}+{sam2.sequence}",
    )
    return s


def print_sample(sample):
    print(">" + sample.header)
    print(sample.sequence)


def print_mixed_samples(mfsl_filepath, number_of_mixed=10, days=7.0):
    """Add mixed samples generated by taking two random samples
    and turning their differences into ns.

    Output samples to stdout."""

    all_samples = parse_mfsl(mfsl_filepath)
    all_ok_sample_pairs = generate_all_sample_pairs(all_samples, days)
    picked_pairs = random.sample(all_ok_sample_pairs, k=number_of_mixed)

    for picked_pair in picked_pairs:
        m = mix_pair(picked_pair)
        # print(m)
        print_sample(m)


if __name__ == "__main__":
    argh.dispatch_command(print_mixed_samples)
