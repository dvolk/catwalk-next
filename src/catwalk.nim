import tables
import intsets
import strutils
import times
import jsony
import algorithm
import system

import symdiff

type
  Sequence = string

  CompressedSequence* = array[4, seq[int]]

  SampleStatus* = enum
    Unknown
    InvalidLength
    TooManyNs
    Ok

  Sample* = tuple
    status: SampleStatus
    diffsets: CompressedSequence
    n_positions: IntSet

  Mask* = tuple
    name: string
    positions: IntSet

  Settings* = tuple
    max_distance: int
    max_n_positions: int

  CatWalk* = tuple
    settings: Settings
    name: string
    reference_name: string
    reference_sequence: string
    mask: Mask
    active_samples: TableRef[int, Sample]
    all_sample_indexes: TableRef[string, int]
    all_sample_names: TableRef[int, string]

#
# CompressedSequence
#

proc empty_compressed_sequence(cs: var CompressedSequence) =
  for i in 0..3: cs[i] = @[]

proc compressed_sequence_counts*(cs: CompressedSequence) : seq[int] =
  result = @[]
  for i in 0..3: result.add(len(cs[i]))

proc count_diff2(cs1: CompressedSequence, cs2: CompressedSequence, sample1_n_positions: IntSet, sample2_n_positions: IntSet, max_distance: int) : int =
  return sum_sym_diff1(cs1[0], cs2[0],
                       cs1[1], cs2[1],
                       cs1[2], cs2[2],
                       cs1[3], cs2[3],
                       sample1_n_positions, sample2_n_positions,
                       max_distance)

proc ref_snp_distance(cs: CompressedSequence) : int =
  result = 0
  for i in 0..3:
    result += len(cs[i])

proc add_position(cs: var CompressedSequence, base: char, position: int) {.inline.} =
  let index = case base:
    of 'A': 0
    of 'C': 1
    of 'G': 2
    of 'T': 3
    else: 4
  cs[index].add(position)

#
# Sample
#

proc new_Sample*(): Sample =
  empty_compressed_sequence(result.diffsets)
  result.n_positions = initIntSet()
  result.status = Unknown

proc is_n_position(c: char): bool {.inline.} =
  c != 'A' and c != 'C' and c != 'G' and c != 'T'

proc reference_compress*(sample_sequence: string, ref_sequence: string, mask: Mask, max_n_positions: int): Sample =
  var
    sample = new_Sample()

  if sample_sequence.len != ref_sequence.len:
    sample.status = InvalidLength
    return

  for i in 0..ref_sequence.high:
    if sample_sequence[i] != ref_sequence[i] and not mask.positions.contains(i):
      if is_n_position(sample_sequence[i]):
        sample.n_positions.incl(i)
      else:
        sample.diffsets.add_position(sample_sequence[i], i)

  if sample.n_positions.len > max_n_positions:
    sample.status = TooManyNs
    empty_compressed_sequence(sample.diffsets)
    sample.n_positions = initIntSet()
  else:
    sample.status = Ok
  return sample

#
# Mask
#

proc new_Mask*(mask_name: string, mask_str: string): Mask =
  result.name = mask_name
  for line in mask_str.splitLines():
    try:
      result.positions.incl(parseInt(line))
    except ValueError:
      echo "Not an integer: '" & line & "'"

#
# CatWalk
#

proc new_CatWalk*(name: string, reference_name: string, reference_sequence: string, mask: Mask) : CatWalk =
  result.name = name
  result.reference_name = reference_name
  result.reference_sequence = reference_sequence
  result.mask = mask
  result.settings.max_distance = 20
  result.settings.max_n_positions = 130000
  result.active_samples = newTable[int, Sample]()
  result.all_sample_indexes = newTable[string, int]()
  result.all_sample_names = newTable[int, string]()

proc process_neighbours(c: var CatWalk, sample1: Sample, sample1_index: int, distance: int): seq[(int, int)] =
  if sample1.status != Ok:
    return
  for sample2_index in c.active_samples.keys:
    if sample2_index == sample1_index:
      continue
    let
      sample2 = c.active_samples[sample2_index]
    if sample2.status != Ok:
      continue
    let
      d = count_diff2(sample1.diffsets, sample2.diffsets, sample1.n_positions, sample2.n_positions, distance)
    if d <= distance:
      result.add((sample2_index, d))

proc get_neighbours*(c: var CatWalk, sample_name: string, distance: int) : seq[(string, int)] =
  let time1 = cpuTime()
  let
    sample_index = c.all_sample_indexes[sample_name]
    sample = c.active_samples[sample_index]
    neighbours = c.process_neighbours(sample, sample_index, distance)
  let dt = cpuTime() - time1
  let sam_num = c.active_samples.len
  echo "Performed " & $sam_num & " distance " & $distance & " comparisons on sample \"" & sample_name & "\" in " & $dt & " seconds (~" & $((1.0 / (dt.float32 / sam_num.float32)) / 1000).int & "k per second)"

  result = @[]
  for (neighbour_index, distance) in neighbours:
    result.add((c.all_sample_names[neighbour_index], distance))

let measure = false

proc add_sample*(c: var CatWalk, name: string, sequence: string, keep: bool) =
  let time1 = cpuTime()
  var
    sample = reference_compress(sequence, c.reference_sequence, c.mask, c.settings.max_n_positions)

  let sample_index = len(c.all_sample_indexes)
  c.all_sample_indexes[name] = sample_index
  c.all_sample_names[sample_index] = name

  if keep:
    c.active_samples[sample_index] = sample

  if measure:
    let l = len(c.all_sample_indexes)
    if (l < 1000 and l %% 100 == 0) or l %% 1000 == 0:
      let
        dt1 = cpuTime() - time1
        time2 = cpuTime()
        n = len(c.get_neighbours(c.all_sample_names[0], 20))
        dt2 = cpuTime() - time2
        mem = getOccupiedMem()
      echo $l & " " & $dt1 & " " & $dt2 & " " & $mem & " " & $n

proc add_sample_from_refcomp*(c: var CatWalk, name: string, refcomp_json: string, keep: bool) =
  let
    tbl = refcomp_json.fromJson(Table[string, seq[int]])
  var
    sample = new_Sample()
  sample.status = Ok
  sample.n_positions = toIntSet(tbl["N"])
  sample.diffsets[0] = tbl["A"]
  sample.diffsets[1] = tbl["C"]
  sample.diffsets[2] = tbl["G"]
  sample.diffsets[3] = tbl["T"]
  sample.diffsets[0].sort()
  sample.diffsets[1].sort()
  sample.diffsets[2].sort()
  sample.diffsets[3].sort()

  let sample_index = len(c.all_sample_indexes)
  c.all_sample_indexes[name] = sample_index
  c.all_sample_names[sample_index] = name

  if keep:
    c.active_samples[sample_index] = sample

proc add_samples_from_refcomp_array*(c: var CatWalk, names: string, refcomps: string) =
  let
    namesjs = names.fromJson(seq[string])
    refcompjs = refcomps.fromJson(seq[Table[string, seq[int]]])
  var i = 0
  for refcompj in refcompjs:
    var sample = new_Sample()
    let name = namesjs[i]
    sample.status = Ok
    sample.n_positions = toIntSet(refcompj["N"])
    sample.diffsets[0] = refcompj["A"]
    sample.diffsets[1] = refcompj["C"]
    sample.diffsets[2] = refcompj["G"]
    sample.diffsets[3] = refcompj["T"]
    sample.diffsets[0].sort()
    sample.diffsets[1].sort()
    sample.diffsets[2].sort()
    sample.diffsets[3].sort()
    let sample_index = len(c.all_sample_indexes)
    c.all_sample_indexes[name] = sample_index
    c.all_sample_names[sample_index] = name
    c.active_samples[sample_index] = sample
    i = i + 1

#
# test
#
when isMainModule:
  let
    mask = new_Mask("test", "0")
    rs = "AAACGT"
  var
    c = new_CatWalk("testcw", "testref", rs, mask)

  c.add_sample("s0", "AAACGT", true)
  c.add_sample("s1", "AAACGT", true)
  c.add_sample("s2", "AAACGC", true)

  assert c.get_neighbours("s0", -1) == []
  assert c.get_neighbours("s0", 0) == [("s1", 0)]
  assert c.get_neighbours("s0", 10) == [("s1", 0), ("s2", 1)]

  c.add_sample_from_refcomp("s3", """{"A": [], "C": [], "G": [], "T": [], "N": []}""", true)

  assert c.get_neighbours("s3", 10) == [("s1", 0),
                                        ("s2", 1),
                                        ("s0", 0)]

  echo "Tests passed."
