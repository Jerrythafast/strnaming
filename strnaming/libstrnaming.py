#!/usr/bin/env python3

#
# Copyright (C) 2020 Jerry Hoogenboom
#
# This file is part of STRNaming, an algorithm for generating simple,
# informative names for sequenced STR alleles in a standardised and
# automated manner.
#
# STRNaming is free software: you can redistribute it and/or modify it
# under the terms of the GNU Lesser General Public License as published
# by the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# STRNaming is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with STRNaming.  If not, see <http://www.gnu.org/licenses/>.
#

import math, re, sys, time

from functools import reduce

import json #FIXME, temp

# Version 2020.12.15.01

TRIM_FOR_SINGLETONS = True  # STRNaming was trained with 'True', but JS version always had 'False'. Not sure what's better.

# Maximum number of seconds that get_best_path() may spend.
MAX_SECONDS = 300

# Number of times a repeat must be repeated to make a 'significant' repeat stretch.
MANY_TIMES = 4

# STRNaming settings.
# TODO: Split these up into separate globals.
# TODO: Put alignment settings here?
NAMING_OPTIONS = {
    "min_repeats": 2,
    "min_repeat_length": 8,
    "min_structure_length": 20,
    "max_unit_length": 6,
    "max_gap": 8,
    "max_long_gap": 20,
    "bases_covered_factor": 0.15874379,
    "units_used_factor": -10.1786473,
    "units_used_multiplier": 1.09052798,
    "repeats_factor": 3.49237881,
    "preferred_repeats_factor": 0.65040629,
    "interruptions_factor": -9.56645361,
    "interruptions_multiplier": 1.41646677,
    "nice_interruptions_factor": 7.27483601,
    "interruption_bases_factor": -0.56939138,
    "prefix_suffix_factor": -1.785957
}

# Pattern that matches a CE allele number.
PAT_CE = re.compile("^(\d+)(?:[\.,](\d+))?$")

# Complementary bases.
COMPL = {"A": "T", "T": "A", "G": "C", "C": "G"}


class OutOfTimeException(Exception):
    pass
class OscillationException(Exception):
    pass


def print_options():
    print("%r" % (NAMING_OPTIONS,))
#print_options


def reverse_complement(sequence):
    """Return the reverse complement of the given DNA sequence."""
    return "".join(COMPL.get(x, x) for x in reversed(sequence))
#reverse_complement


def align(haystack, needle, prefix=False, debug=False):
    """
    Find the closest match of needle in haystack.
    If prefix is true, the match will start at the start of haystack.
    """
    try:
        return align.cache[haystack, needle, prefix]
    except KeyError:
        pass

    row_offset = len(needle) + 1
    matrix_match = [0] * row_offset * (len(haystack) + 1)
    matrix_direction = [0] * row_offset * (len(haystack) + 1)

    # Arrow enum constants.
    A_MATCH = 0b00001
    A_HORZ  = 0b00010
    A_VERT  = 0b00100

    # Settings.
    MATCH_SCORE = 1
    MISMATCH_SCORE = -1
    GAP_SCORE = -1

    max_score = len(needle) * GAP_SCORE
    max_pos = 0

    for i in range(1, len(matrix_match)):
        x = i % row_offset
        y = i // row_offset

        # Initialisation of first row and column.
        if x == 0:
            matrix_direction[i] |= A_VERT  # Can go vertically.
            if prefix:
                matrix_match[i] = y * GAP_SCORE
            continue
        if y == 0:
            matrix_direction[i] |= A_HORZ  # Can go horizontally.
            matrix_match[i] = x * GAP_SCORE
            continue

        if needle[x - 1] == haystack[y - 1]:
            match = MATCH_SCORE
        else:
            match = MISMATCH_SCORE

        options = [
            matrix_match[i - 1 - row_offset] + match,
            matrix_match[i - 1] + GAP_SCORE,
            matrix_match[i - row_offset] + GAP_SCORE]
        matrix_match[i] = max(options)
        if options[0] == matrix_match[i]:
            matrix_direction[i] |= A_MATCH  # Can go diagonally.
        if options[1] == matrix_match[i]:
            matrix_direction[i] |= A_HORZ  # Can go horizontally.
        if options[2] == matrix_match[i]:
            matrix_direction[i] |= A_VERT  # Can go vertically.
        if x == len(needle) and matrix_match[i] > max_score:
            max_score = matrix_match[i]
            max_pos = i

    if debug:
        print("MATCH")
        for i in range(0, len(matrix_match), row_offset):
            print(("%5i" * row_offset) % tuple(matrix_match[i : i + row_offset]))
        print("FLAGS")
        for i in range(0, len(matrix_direction), row_offset):
            print(("%3s|" * row_offset) % tuple("".join((
                "H" if x & A_HORZ else " ",
                "D" if x & A_MATCH  else " ",
                "V" if x & A_VERT else " "
            )) for x in matrix_direction[i : i + row_offset]))
        print("Traceback")


    # Backtracking.
    i = max_pos
    while (prefix and i) or i % row_offset:
        x = i % row_offset
        y = i // row_offset
        if debug:
            print("%i = (%i, %i)" % (i, x, y))

        if matrix_direction[i] & A_HORZ:
            # Go horizontally.  Deletion.
            i -= 1
        elif matrix_direction[i] & A_VERT:
            # Go vertically.  Insertion.
            i -= row_offset
        else:
            # Go diagonally.  Either match or mismatch.
            i -= 1 + row_offset

    # TODO: I could return the traceback to call variants from it in the naming phase.
    result = (i // row_offset, max_pos // row_offset, max_score)

    # Store the result in the cache.
    align.cache[haystack, needle, prefix] = result
    return result
align.cache = {}
#align


def detect_repeats(seq):
    """
    Return a matrix of max_unit_length x len(seq) elements.
    The value at position (j,i) gives the length of a repeat of units of
    length j (in number of units), ending in position i in the sequence.
    This function will only accept repeats of min_repeats or more units;
    it will report a length of 1 for any repeat shorter than this.
    """
    len_seq = len(seq)
    matrix = [[0]*len_seq for j in range(NAMING_OPTIONS["max_unit_length"])]
    for i in range(len_seq):
        for j in range(NAMING_OPTIONS["max_unit_length"]):
            if i <= j:
                matrix[j][i] = i + 1  # No full unit yet.
            elif seq[i] != seq[i - j - 1]:
                matrix[j][i] = j + 1  # Restart counting.
            else:
                matrix[j][i] = matrix[j][i - 1] + 1
    for i in range(len_seq):
        for j in range(NAMING_OPTIONS["max_unit_length"]):
            matrix[j][i] //= (j + 1)
            if 0 < matrix[j][i] < NAMING_OPTIONS["min_repeats"]:
                matrix[j][i] = 1;  # Ignore, too short.
    return matrix
#detect_repeats


def find_longest_repeat(matrix):
    """
    Get start and end position and unit length of the longest repeat in
    matrix.  Break ties by choosing the shortest unit, the earliest in
    the sequence (lowest indices into matrix).
    """
    maxes = [(((max(row)-1) or -1) + 1) * (j + 1) for j, row in enumerate(matrix)]
    max_len = max(maxes)
    max_unit = maxes.index(max_len)
    max_end = matrix[max_unit].index(max(matrix[max_unit])) + 1
    max_start = max_end - max_len
    return max_start, max_end, max_unit + 1
#find_longest_repeat


def recurse_find_longest_repeat(seq, matrix=None):
    """Generate tuples of (repeat length, unit length)."""
    # Find the longest repeat in the sequence.
    if matrix is None:
        matrix = detect_repeats(seq)
    start, end, unit = find_longest_repeat(matrix)

    # Degenerate case: no repeat found.  Short-circuit.
    length = end - start
    if length < NAMING_OPTIONS["min_repeat_length"] or length < NAMING_OPTIONS["min_repeats"] * unit:
        yield (len(seq), len(seq))
        return

    # First, recursively yield any parts prior to the found repeat.
    if start:
        for prefix in recurse_find_longest_repeat(seq[:start], [row[:start] for row in matrix]):
            yield prefix

    # Then, yield the repeat that was found in this call.
    yield (length, unit)

    # Finally, recursively yield any parts after the found repeat.
    if end < len(seq):
        for suffix in recurse_find_longest_repeat(seq[end:], detect_repeats(seq[end:])):
            yield suffix
#recurse_find_longest_repeat


def find_repeated_units(seq, cycle_units):
    """Return a dict of repeat units found in seq."""
    units = {}
    pos = 0
    for r_len, u_len in recurse_find_longest_repeat(seq):
        if r_len > u_len:
            # Repeated element.
            unit = seq[pos : pos + u_len]
            repeat = (pos, r_len)
            try:
                units[unit].append(repeat)
            except KeyError:
                units[unit] = [repeat]
        pos += r_len

    if not cycle_units:
        return units

    # Cycle units.
    newunits = {unit: repeats.copy() for unit, repeats in units.items()}
    for unit, repeats in units.items():
        for i in range(1, len(unit)):
            newunit = unit[i:] + unit[:i]
            regex = re.compile("^(" + newunit + "){2,}")
            for pos, r_len in repeats:
                start = pos + i
                newmatch = regex.match(seq[start : start + r_len])
                if newmatch:
                    newrepeat = (start, newmatch.end(0))
                    try:
                        newunits[newunit].append(newrepeat)
                    except KeyError:
                        newunits[newunit] = [newrepeat]
    return newunits
#find_repeated_units


def find_overlong_gap(units):
    """Return the longest interruption, if that interruption is longer than max_gap."""
    longest_gap = NAMING_OPTIONS["max_gap"]
    longest_start_pos = 0
    for i in range(1, len(units)):
        start_pos = units[i - 1][1]
        gap_length = units[i][0] - start_pos
        if gap_length > longest_gap:
            longest_gap = gap_length
            longest_start_pos = start_pos
    return (longest_start_pos, longest_start_pos + longest_gap) if longest_start_pos else None
#find_overlong_gap


def powersum(multiplier, count):
    """
    Compute sum(multiplier ** i for i in range(count)) efficiently.
    """
    power = 1
    powersum = 0
    for i in range(count):
        powersum += power
        power *= multiplier
    return powersum
#powersum


def calc_path_score(len_prefix, len_suffix, len_seq, len_large_gap, path, block_length, preferred_units=None,
                    is_refseq=False):
    """Calculate the score of the given path.  The higher the better."""
    if preferred_units is None:
        preferred_units = set()

    # These are the things we measure.
    bases_covered = 0
    units_used = set() if is_refseq else set(preferred_units)
    repeats = 0
    preferred_repeats = 0
    interruptions = []
    prefix_delta = len_prefix
    suffix_delta = len_suffix

    # Measure things about the path.
    pos = 0
    for start, end, unit, (repeat_length, unit_length, anchor, preferred) in path:
        repeat_count = repeat_length // unit_length
        bases_covered += repeat_length
        repeats += repeat_count
        if preferred:
            preferred_repeats += repeat_count
        units_used.add(unit)
        if start and not pos:
            prefix_delta -= start
        elif start > pos:
            interruptions.append(start - pos)
        pos = end
    if pos < len_seq:
        suffix_delta -= len_seq - pos

    # Count interruptions of the most common unit length.
    nice_interruptions = sum(1 for x in interruptions if x == block_length)
    large_gap_delta = len_large_gap - max(interruptions + [0]) if len_large_gap else 0

    # Calculate score.
    return sum((
        bases_covered * NAMING_OPTIONS["bases_covered_factor"],
        powersum(NAMING_OPTIONS["units_used_multiplier"], len(units_used)) *
            NAMING_OPTIONS["units_used_factor"],
        repeats * NAMING_OPTIONS["repeats_factor"],
        preferred_repeats * NAMING_OPTIONS["preferred_repeats_factor"],
        powersum(NAMING_OPTIONS["interruptions_multiplier"], len(interruptions)) *
            NAMING_OPTIONS["interruptions_factor"],
        nice_interruptions * NAMING_OPTIONS["nice_interruptions_factor"],
        sum(interruptions) * NAMING_OPTIONS["interruption_bases_factor"],
        0 if is_refseq else (abs(prefix_delta) + abs(suffix_delta) + abs(large_gap_delta)) *
            NAMING_OPTIONS["prefix_suffix_factor"]
    ))
#calc_path_score


def generate_scaffolds(repeats, scaffold=None, infos=None, anchored=None, orphans=None, start=0):
    if scaffold is None:
        for i, repeat in enumerate(repeats[start:], start):
            unit_len = len(repeat[2])
            # max_gap_len is equal to the unit length for non-preferred singletons, N/A otherwise.
            max_gap_len = unit_len if not repeat[3][3] and unit_len == repeat[1] - repeat[0] else 0
            new_scaffold = [repeat]
            new_infos = repeat[3][:2]
            new_anchored = {repeat[2]} if repeat[3][2] else set()
            new_orphans = set() if repeat[3][2] else {repeat[2]}
            yield (new_scaffold, new_infos, new_anchored, new_orphans, max_gap_len)
            yield from generate_scaffolds(repeats, new_scaffold, new_infos, new_anchored, new_orphans, i)
        return

    scaffold_end_pos = scaffold[-1][1]
    scaffold_end_unit = scaffold[-1][2]
    for i, repeat in enumerate(repeats[start:], start):
        # We're not putting repeats of the same unit adjacently. (Those
        # have been split up by upstream code to allow more junctions.)
        unit = repeat[2]
        if repeat[0] == scaffold_end_pos and unit != scaffold_end_unit:
            new_scaffold = scaffold + [repeat]
            new_infos = max(infos, repeat[3][:2])
            new_anchored = anchored
            new_orphans = orphans
            if repeat[3][2]:
                if unit not in anchored:
                    new_anchored = new_anchored | {unit}
                    if unit in orphans:
                        new_orphans = new_orphans - {unit}
            elif unit not in anchored and unit not in orphans:
                new_orphans = new_orphans | {unit}
            yield (new_scaffold, new_infos, new_anchored, new_orphans, 0)
            yield from generate_scaffolds(repeats, new_scaffold, new_infos, new_anchored, new_orphans, i)
#generate_scaffolds


def get_scaffolds(repeats):
    scaffolds = {}
    for scaffold in generate_scaffolds(repeats):
        start = scaffold[0][0][0]
        end = scaffold[0][-1][1]
        if start not in scaffolds:
            scaffolds[start] = {end: [scaffold]}
        elif end not in scaffolds[start]:
            scaffolds[start][end] = [scaffold]
        else:
            scaffolds[start][end].append(scaffold)
    return scaffolds
#get_scaffolds


def get_gaps(repeats, scaffolds):
    short_gaps = {}
    all_gaps = {}
    for start in set(repeat[1] for repeat in repeats):
        for end in set(repeat[0] for repeat in repeats):
            size = end - start
            if size <= 0 or (start in scaffolds and end in scaffolds[start]):
                continue  # Illegal gap: shoots hole in scaffold.
            if size <= NAMING_OPTIONS["max_long_gap"]:
                try:
                    all_gaps[start].add(end)
                except KeyError:
                    all_gaps[start] = {end}
            if size <= NAMING_OPTIONS["max_gap"]:
                try:
                    short_gaps[start].add(end)
                except KeyError:
                    short_gaps[start] = {end}
    return short_gaps, all_gaps
#get_gaps


def extend_ranges_with_gaps(ranges, gap_ranges):
    extended = {start: ends.copy() for start, ends in ranges.items()}
    for start, ends in ranges.items():
        subtree = extended[start]
        for gap_start, anchorables in ends.items():
            if gap_start in gap_ranges:
                for gap_end in gap_ranges[gap_start]:
                    for final_end, more_anchorables in ranges[gap_end].items():
                        subtree[final_end] = subtree.get(final_end, set()) | anchorables | more_anchorables
    return extended
#extend_ranges_with_gaps


def get_ranges(scaffolds, short_gaps, all_gaps, is_refseq):
    # TODO: This function is now hardcoded for MAX_GAPS=5. Maybe generalise?
    # NOTE: this function may report ranges that are not acutally possible due to
    # rules regarding e.g. singletons.
    # Lookup tables for starting positions and associated ending positions.
    # The positions are linked with the units that are anchorable along the way.
    start_end00 = {start:
            {end: reduce(set.union, (scaffold[2] for scaffold in scaffolds_here))
            for end, scaffolds_here in ends.items()}
        for start, ends in scaffolds.items()}
    start_end10 = extend_ranges_with_gaps(start_end00, short_gaps)
    start_end20 = extend_ranges_with_gaps(start_end10, short_gaps)
    start_end30 = extend_ranges_with_gaps(start_end20, short_gaps)
    start_end40 = extend_ranges_with_gaps(start_end30, short_gaps)
    if is_refseq:
        start_end50 = extend_ranges_with_gaps(start_end40, short_gaps)
    else:
        start_end11 = extend_ranges_with_gaps(start_end00, all_gaps)
        start_end21 = extend_ranges_with_gaps(start_end11, short_gaps)
        start_end31 = extend_ranges_with_gaps(start_end21, short_gaps)
        start_end41 = extend_ranges_with_gaps(start_end31, short_gaps)
        start_end51 = extend_ranges_with_gaps(start_end41, short_gaps)

    # Now, build a tree from these.
    node00 = (start_end00, {}, None, None)
    node10 = (start_end10, short_gaps, node00, None)
    node20 = (start_end20, short_gaps, node10, None)
    node30 = (start_end30, short_gaps, node20, None)
    node40 = (start_end40, short_gaps, node30, None)
    if is_refseq:
        return (start_end50, short_gaps, node40, None)
    node11 = (start_end11, all_gaps, node00, node00)
    node21 = (start_end21, all_gaps, node11, node10)
    node31 = (start_end31, all_gaps, node21, node20)
    node41 = (start_end41, all_gaps, node31, node30)
    node51 = (start_end51, all_gaps, node41, node40)
    return node51
#get_ranges


def recurse_gen_path(start_pos, end_pos, scaffolds, ranges,
                     anchored=set(), orphans=set(), prev_gap_len=-1, prev_unit=""):
    scaffolds_here = scaffolds[start_pos]
    try:
        for scaffold, longest_stretch, anchored2, orphans2, max_gap_len in scaffolds_here[end_pos]:
            if (not prev_unit and not scaffold[0][3][3]) or not scaffold[-1][3][3]:
                continue  # Can't start or end with a non-ref repeat unit.
            if orphans - anchored2 or orphans2 - anchored:
                continue  # Would end up with orphans.
            if not max_gap_len or max_gap_len >= prev_gap_len or scaffold[-1][2] == prev_unit:
                yield (scaffold, longest_stretch)
    except KeyError:
        pass  # No scaffolds reach exactly from start_pos to end_pos.

    # Insert a scaffold that reaches not quite to end_pos just yet.
    gaps = ranges[1]
    for gap_pos, scaffold_list in scaffolds_here.items():

        if gap_pos not in gaps:
            continue

        # Choose the next gap and continue the structure on the other side.
        for next_pos in gaps[gap_pos]:
            gap_len = next_pos - gap_pos
            min_singleton_len = gap_len if prev_gap_len == -1 or gap_len < prev_gap_len else prev_gap_len
            remaining_ranges = ranges[2 if gap_len <= NAMING_OPTIONS["max_gap"] else 3]

            try:
                anchorable_after_gap = remaining_ranges[0][next_pos][end_pos]
            except KeyError:
                continue  # Can't reach from next_pos to end_pos.

            for scaffold, longest_stretch2, anchored2, orphans2, max_gap_len in scaffold_list:
                last_unit = scaffold[-1][2]
                if not prev_unit and not scaffold[0][3][3]:
                    continue  # Can't start with this non-ref repeat unit.
                if max_gap_len and max_gap_len < min_singleton_len and not last_unit == prev_unit:
                    continue  # This singleton would get too isolated.
                now_anchored = anchored | anchored2
                now_orphaned = (orphans | orphans2) - now_anchored
                if not now_orphaned - anchorable_after_gap:
                    for substructure, longest_stretch1 in recurse_gen_path(next_pos, end_pos,
                            scaffolds, remaining_ranges, now_anchored, now_orphaned, gap_len,
                            last_unit):
                        yield (scaffold + substructure, max(longest_stretch1, longest_stretch2))
#recurse_gen_path


def gen_valid_paths(start_pos, end_pos, scaffolds, ranges, is_refseq):
    # Most validity rules are efficiently implemented in recurse_gen_path.
    for result in recurse_gen_path(start_pos, end_pos, scaffolds, ranges):
        path, (longest_stretch, block_length) = result

        # The path is invalid if it uses a repeat unit that is longer
        # than the dominant block length, or it has only singletons of a
        # non-preferred unit.
        invalid_unit_length = False
        repeated_units = set()
        singletons = set()
        for start, end, unit, (repeat_length, unit_length, anchor, preferred) in path:
            if unit_length > block_length:
                invalid_unit_length = True
                break
            if is_refseq or not preferred:
                if repeat_length > unit_length:
                    repeated_units.add(unit)
                else:
                    singletons.add(unit)
        if not invalid_unit_length and not singletons - repeated_units:
            if not is_refseq or longest_stretch >= block_length * MANY_TIMES:
                # For reference sequences, the longest stretch must be significant.
                yield result
#gen_valid_paths


def gen_all_paths(prefix, suffix, seq, repeats, is_refseq):
    scaffolds = get_scaffolds(repeats)
    short_gaps, all_gaps = get_gaps(repeats, scaffolds)
    ranges = get_ranges(scaffolds, short_gaps, all_gaps, is_refseq)

    # A *preferred* repeat unit is required as the starter/finisher.
    start_positions = set(x[0] for x in repeats if x[3][3])
    end_positions = set(x[1] for x in repeats if x[3][3])

    if not is_refseq:
        # Try anchoring the prefix and suffix.
        success = False
        if seq.startswith(prefix):
            len_prefix = len(prefix)
            has_prefix = len_prefix in start_positions
        else:
            zero, len_prefix, score = align(seq, prefix, prefix=True)
            has_prefix = score > 0 and len_prefix in start_positions
        if seq.endswith(suffix):
            pos_suffix = len(seq) - len(suffix)
            has_suffix = pos_suffix in end_positions
        else:
            zero, len_suffix, score = align(seq[-1::-1], suffix[-1::-1], prefix=True)
            pos_suffix = len(seq) - len_suffix
            has_suffix = score > 0 and pos_suffix in end_positions

        # First, try to reach exactly from prefix to suffix.
        if has_prefix and has_suffix and pos_suffix in ranges[0][len_prefix]:
            for result in gen_valid_paths(len_prefix, pos_suffix, scaffolds, ranges, is_refseq):
                success = True
                yield result
            if success:
                return  # Yielded path with anchored prefix AND suffix.

        if has_prefix:
            start_positions.remove(len_prefix)
        if has_suffix:
            end_positions.remove(pos_suffix)

        # If that did not work, force starting at the end of prefix
        # or finishing at the start of suffix separately.
        if has_prefix:
            for end_pos in end_positions & ranges[0][len_prefix].keys():
                for result in gen_valid_paths(len_prefix, end_pos, scaffolds, ranges, is_refseq):
                    success = True
                    yield result
        if has_suffix:
            for start_pos in start_positions:
                if pos_suffix in ranges[0][start_pos]:
                    for result in gen_valid_paths(start_pos, pos_suffix, scaffolds, ranges, is_refseq):
                        success = True
                        yield result
        if success:
            return  # Yielded path with anchored prefix OR suffix.

    for start_pos in start_positions:
        for end_pos in end_positions & ranges[0][start_pos].keys():
            if not is_refseq or end_pos - start_pos >= NAMING_OPTIONS["min_structure_length"]:
                yield from gen_valid_paths(start_pos, end_pos, scaffolds, ranges, is_refseq)
#gen_all_paths


def get_best_path(prefix, suffix, seq, repeats, preferred_units, ref_len_large_gap=-1):
    """Return combination of repeats that maximises calc_path_score."""
    starttime = time.monotonic()

    len_prefix = len(prefix) if prefix is not None else -1
    len_suffix = len(suffix) if suffix is not None else -1
    len_seq = len(seq)
    is_refseq = ref_len_large_gap == -1
    best_score = -sys.maxsize
    best_path = []

    for path, (_, block_length) in gen_all_paths(prefix, suffix, seq, repeats, is_refseq):
        # Timekeeping.
        if time.monotonic() - starttime > MAX_SECONDS:
            raise OutOfTimeException()

        # Calculate the score of the current path.
        score = calc_path_score(len_prefix, len_suffix, len_seq, ref_len_large_gap, path,
            block_length, preferred_units, is_refseq)

        if score >= best_score:
            if score == best_score:
                # Only update path if it is more 5'.
                better_path = False
                for a, b in zip(best_path, path):
                    if a[0] == b[0]:
                        continue
                    better_path = a[0] > b[0]
                    break
                if not better_path:
                    continue  # best_path is 5' of path.
            best_score = score
            best_path = path

    return best_score, [repeat[:3] for repeat in best_path]
#get_best_path


def find_repeat_stretches(seq, units, allow_bridges, allow_one, repeats=None):
    """Return a list of repeat stretches. The input units must be sorted by decreasing length."""
    # Find all stretches of each repeated unit.
    if repeats is None:
        repeats = []
    for unit in units:
        matches = []

        # Compile a list of all places where this unit is found.
        # len(stretch) must be >= 4.
        # NOTE: Hardcoded minimum length.
        regex = re.compile("(" + unit + "){" + str(4 if len(unit) < 2 else 2 if len(unit) < 4 or not allow_one else 1) + ",}")

        pos = 0
        while True:
            match = regex.search(seq[pos:])
            if match is None:
                break
            matches.append([pos + match.start(), len(match.group(0))])
            pos += match.start() + len(match.group(0)) - len(unit) + 1

        # Try bridging gaps to a single additional unit.
        if allow_bridges and len(unit) > 2 and (not allow_one or len(unit) < 4):
            for start, length in tuple(matches):
                end = start + length
                if start >= len(unit)*2 and seq[start - len(unit)*2 : start - len(unit)] == unit and not seq[:start - len(unit)*2].endswith(unit):
                    # Bridge one unit before stretch.
                    matches.append([start - len(unit)*2, len(unit)])
                if end + len(unit)*2 <= len(seq) and seq[end + len(unit) : end + len(unit)*2] == unit and seq[end + len(unit)*2 : end + len(unit)*3] != unit:
                    # Bridge one unit after stretch.
                    matches.append([end + len(unit), len(unit)])

        for match_start, match_len in matches:
            # Ignore sub-repeats within significant stretches of a longer unit.
            # Using those is not expected to give a high score, but their
            # presence may greatly increase calculation time.
            def is_overshadowed_by(longer_unit):
                if len(longer_unit[2]) <= len(unit):
                    return False  # Not a longer repeat unit.
                if match_len // len(unit) > MANY_TIMES:
                    return False  # This repeat is itself significant.
                if max(match_start, longer_unit[0]) >= min(match_start + match_len, longer_unit[1]):
                    return False  # The longer unit is not overlapping this repeat.
                if unit not in longer_unit[2] + longer_unit[2]:
                    return False  # The longer repeat unit does not contain the current one.
                repeated_times = (longer_unit[1] - longer_unit[0]) // len(longer_unit[2])
                if repeated_times < MANY_TIMES or (repeated_times == MANY_TIMES and unit not in longer_unit[2]):
                    return False  # The longer unit is not repeated a lot.
                # The longer unit is repeated a lot; ignore the current one.
                return True
            if not any(is_overshadowed_by(longer_unit) for longer_unit in repeats):
                # There is no significat repeat of a longer unit that overlaps with this match.
                repeats.append([match_start, match_start + match_len, unit])

    # Sort repeats.
    repeats.sort()

    # Remove singletons that fall completely within a significant repeat stretch
    # of a shorter unit that is embedded in the singleton.
    # For example, remove singletons of 'ACTA' within a significant repeat stretch of 'ACT'.
    i = 0
    while i < len(repeats):
        if (len(repeats[i][2]) == repeats[i][1] - repeats[i][0] and  # Current is singleton
                any(len(o_unit) < len(repeats[i][2]) and  # Other is shorter unit
                    repeats[i][0] >= o_start and repeats[i][1] <= o_end and  # Overlap
                    (o_end - o_start) // len(o_unit) > MANY_TIMES  # Other is significant repeat
                    for j, (o_start, o_end, o_unit) in enumerate(repeats))):
            # The current (singleton) unit overlaps completely with a significant repeat.
            del repeats[i]
        else:
            i += 1

    return repeats
#find_repeat_stretches


def find_everything(seq, unit_locations):  # FIXME, awful name.
    """
    The provided unit_locations are taken to represent the preferred units.
    Unit search will start around these; locations of newly found units
    will be added to the object passed in.
    """
    # Sort preferred units, longest come first.
    preferred_units = sorted(unit_locations, key=len, reverse=True)

    # Find stretches of the preferred units first,
    # then find more units outside those stretches.
    repeats = find_repeat_stretches(seq, preferred_units, True, True)
    for i, (start, end, preferred_unit) in enumerate(repeats):
        if start:
            # Find units before repeat i.
            locations = find_repeated_units(seq[:start], False)
            for unit in locations:
                try:
                    unit_locations[unit] += locations[unit]
                except KeyError:
                    unit_locations[unit] = locations[unit]
        if len(seq) - end:
            # Find units after repeat i.
            locations = find_repeated_units(seq[end:], False)
            for unit in locations:
                new_locations = [(x[0] + end, x[1]) for x in locations[unit]]
                try:
                    unit_locations[unit] += new_locations
                except KeyError:
                    unit_locations[unit] = new_locations
        if i:
            # Find units between repeat i-1 and repeat i.
            prev_end = repeats[i-1][1]
            if prev_end < start:
                locations = find_repeated_units(seq[prev_end : start], False)
                for unit in locations:
                    new_locations = [(x[0] + prev_end, x[1]) for x in locations[unit]]
                    try:
                        unit_locations[unit] += new_locations
                    except KeyError:
                        unit_locations[unit] = new_locations

    # If no stretches of preferred units were found, just find any stretches.
    if not repeats:
        locations = find_repeated_units(seq, False)
        for unit in locations:
            try:
                unit_locations[unit] += locations[unit]
            except:
                unit_locations[unit] = locations[unit]

    # Deduplicate unit_locations.
    for unit in unit_locations:
        unit_locations[unit].sort()
        i = 1
        while i < len(unit_locations[unit]):
            if tuple(unit_locations[unit][i-1]) == tuple(unit_locations[unit][i]):
                del unit_locations[unit][i]
            else:
                i += 1

    # Find all stretches of each newly-discovered repeated unit and add those to the list.
    units = sorted(unit_locations, key=len, reverse=True)
    find_repeat_stretches(
        seq,
        [unit for unit in units if unit not in preferred_units],
        False, False, repeats)

    # Where repeat stretches overlap, also include trimmed copies that don't overlap.
    repeat_trims = tuple(({0}, {0}) for x in range(len(repeats)))
    for i in range(len(repeats)):
        for j in range(i + 1, len(repeats)):
            overlap_len = repeats[i][1] - repeats[j][0]
            if overlap_len > 0:
                # Stretch i extends past the start of stretch j.
                # We can trim the end of stretch i, or the start of stretch j, to compensate.
                if repeats[i][2] in preferred_units and repeats[j][2] in preferred_units and (TRIM_FOR_SINGLETONS or (repeats[i][1] - repeats[i][0] != len(repeats[i][2]) and repeats[j][1] - repeats[j][0] != len(repeats[j][2]))):
                    # Always trimming if both repeats involve preferred units.
                    repeat_trims[i][1].add(math.ceil(float(overlap_len)/len(repeats[i][2]))*len(repeats[i][2]))
                    repeat_trims[j][0].add(math.ceil(float(overlap_len)/len(repeats[j][2]))*len(repeats[j][2]))
                else:
                    # Otherwise, only trimming if doing so leaves no gap.
                    if overlap_len % len(repeats[i][2]) == 0:
                        repeat_trims[i][1].add(overlap_len)
                    if overlap_len % len(repeats[j][2]) == 0:
                        repeat_trims[j][0].add(overlap_len)
    for i in range(len(repeat_trims)):
        for trim_start in repeat_trims[i][0]:
            for trim_end in repeat_trims[i][1]:
                # The -4 in the next line enforces the minimum repeat stretch length of 4nt here.
                # NOTE: hardcoded number.
                if (trim_start or trim_end) and repeats[i][1]-repeats[i][0]-4 > trim_start + trim_end:
                    repeats.append([repeats[i][0] + trim_start, repeats[i][1] - trim_end, repeats[i][2]])
    repeats.sort()

    # Deduplicate repeats.
    i = 1
    while i < len(repeats):
        if tuple(repeats[i-1]) == tuple(repeats[i]):
            del repeats[i]
        else:
            i += 1

    # Add (repeat_length, unit_length, anchor, preferred) to the end.
    for repeat in repeats:
        start, end, unit = repeat
        anchor = any(start < start_l + length and end > start_l
            for start_l, length in unit_locations[unit])
        repeat.append((end - start, len(unit), anchor, unit in preferred_units))

    # Return what we have found.
    return repeats
#find_everything


### Entry point for target analysis ###


def collapse_repeat_units(seq, prefix, suffix, preferred_units, overlong_gap):
    """
    Combine repeat units in the given sequence.
    """
    # Short-circuit if we're not told to collapse anything.
    if not preferred_units:
        return [[seq, 1]]

    # Deduplicate preferred units.
    preferred_units = [unit for i, unit in enumerate(preferred_units)
                       if unit and preferred_units.index(unit) == i]

    # Find all stretches of any unit that stands out.
    unit_locations = {preferred_unit: [(0, len(seq))] for preferred_unit in preferred_units}
    repeats = find_everything(seq, unit_locations)

    # Find best combination of repeat unit usage.
    score, path = get_best_path(prefix, suffix, seq, repeats, preferred_units, len(overlong_gap))
    return path
#collapse_repeat_units


collapse_repeat_units_refseq_three = set()  # List of sequences requiring a third round.
def collapse_repeat_units_refseq(seq):
    """Combine repeat units in the given sequence."""

    # Find some units to start with. These will be treated as 'preferred' in the first round.
    unit_locations = find_repeated_units(seq, True)

    previous_paths = set()  # FIXME, temp
    last_path = None
    large_gap_length = -1
    while True:
        # Make sure the longest repeat units come first.
        preferred_units = sorted(unit_locations, key=len, reverse=True)

        # Short-circuit if we can't collapse anything.
        if not preferred_units:
            return 0, []

        # Find all stretches of each repeated unit.
        # Note: unit_locations will be amended with newly-found units.
        repeats = find_everything(seq, unit_locations)

        # Find the best combination of repeat unit usage.
        prefix = None if not last_path else seq[:last_path[0][0]]
        suffix = None if not last_path else seq[last_path[-1][1]:]
        score, path = get_best_path(prefix, suffix, seq, repeats, preferred_units, large_gap_length)

        # TODO: Log a warning if repeats occur very close to either end of the sequence.
        # The name may not be stable if this is the case.
        """js:
        for ([start,end,unit] of path)
          if (start < 30 || seq.length-end < 30)
            console.log("Too close!", [start,seq.length-end,unit + "(" + ((end-start)/unit.length) + ")"]);
        """

        if last_path != path:
            # FIXME, temp:
            this_path = json.dumps(path)
            if this_path in previous_paths:
                print("last_path=%r" % (json.dumps(last_path),))
                print("this_path=%r" % (this_path,))
                print("previous_paths=%r" % (previous_paths,))
                raise OscillationException()
            previous_paths.add(this_path)
            if len(previous_paths) == 2:
                collapse_repeat_units_refseq_three.add(seq)

            # Repeat find_everything() and get_best_path(), now treating the used units as preferred.
            used_units = set(unit for start, end, unit in path)
            for unit in tuple(unit_locations.keys()):
                if unit not in used_units:
                    del unit_locations[unit]
                else:
                    # The preferred units can go anywhere in the next round.
                    unit_locations[unit] = [(0, len(seq))]
            last_path = path
            large_gap = find_overlong_gap(path)
            large_gap_length = large_gap[1] - large_gap[0] if large_gap else 0
        else:
            return score, path
#collapse_repeat_units_refseq


### Entry point for refseq analysis ###


def recurse_collapse_repeat_units_refseq(seq, *, offset=0):
    score, path = collapse_repeat_units_refseq(seq)
    if path and path[-1][1] - path[0][0] >= NAMING_OPTIONS["min_structure_length"]:

        # Maybe we can make structures before this one.
        pos = path[0][0] - NAMING_OPTIONS["max_long_gap"]
        if pos >= NAMING_OPTIONS["min_structure_length"]:
            yield from recurse_collapse_repeat_units_refseq(seq[:pos], offset=offset)

        yield [[repeat[0] + offset, repeat[1] + offset, repeat[2]] for repeat in path]

        # Maybe we can make structures after this one.
        pos = path[-1][1] + NAMING_OPTIONS["max_long_gap"]
        if len(seq) - pos >= NAMING_OPTIONS["min_structure_length"]:
            yield from recurse_collapse_repeat_units_refseq(seq[pos:], offset=pos+offset)
#recurse_collapse_repeat_units_refseq


# FIXME, get_genome_pos is copypasted from FDSTools but not all its features are used by STRNaming.

def get_genome_pos(location, x, *, invert=False):
    """Get the genome position of the x-th base in a sequence."""
    if invert:
        offset = 0
        for i in range(1, len(location)):
            if i % 2:
                # Starting position.
                pos = location[i]
            elif pos <= x <= location[i]:
                # x is in the current range
                break
            else:
                offset += location[i] - pos + 1
        else:
            if len(location) % 2:
                raise ValueError("Position %i is outside sequence range" % x)
        return offset + x - pos
    else:
        for i in range(1, len(location)):
            if i % 2:
                # Starting position.
                pos = location[i]
            elif location[i]-pos < x:
                # x is after this ending position
                x -= location[i] - pos + 1
            else:
                # x is before this ending position
                break
        return pos + x
#get_genome_pos


# TODO: merge align() into this
def align_affine(template, sequence, *, match_score=1, mismatch_score=-3, gap_open_score=-7, gap_extend_score=-2, cache=True, debug=False):
    """
    Perform a global alignment of sequence to template and return a
    list of variations detected.

    By default, the results of this function are cached.  Set cache to
    False to suppress caching the result and reduce memory usage.

    Setting debug to True will cause the alignment matrices to be
    printed to sys.stdout.  Be aware that they can be quite large.
    """
    try:
        return align_affine.cache[template, sequence, match_score, mismatch_score, gap_open_score, gap_extend_score]
    except KeyError:
        pass

    row_offset = len(template) + 1
    matrix_match = [0] * row_offset * (len(sequence) + 1)
    matrix_gap1 = [-sys.maxsize - 1] * row_offset * (len(sequence) + 1)
    matrix_gap2 = [-sys.maxsize - 1] * row_offset * (len(sequence) + 1)
    matrix_direction = [0] * row_offset * (len(sequence) + 1)

    # Matrix and arrow enum constants.
    M_MATCH = 0
    M_GAP1 = 1
    M_GAP2 = 2
    A_MATCH  = 0b00001
    A_HORZ_O = 0b00010
    A_HORZ_E = 0b00100
    A_VERT_O = 0b01000
    A_VERT_E = 0b10000

    for i in range(len(matrix_match)):
        x = i % row_offset
        y = i // row_offset

        # Initialisation of first row and column.
        if x == 0 or y == 0:
            if x != 0:
                # Top row.
                matrix_gap1[i] = gap_open_score + gap_extend_score * (x - 1)
                matrix_match[i] = matrix_gap1[i]
                matrix_direction[i] = A_HORZ_E | (A_HORZ_O if x == 1 else 0)
            elif y != 0:
                # Left column.
                matrix_gap2[i] = gap_open_score + gap_extend_score * (y - 1)
                matrix_match[i] = matrix_gap2[i]
                matrix_direction[i] = A_VERT_E | (A_VERT_O if y == 1 else 0)
            else:
                # Top left corner.
                matrix_direction[i] = A_MATCH
            continue


        match = match_score if template[x - 1] == sequence[y - 1] else mismatch_score

        options_gap1 = (
            matrix_match[i - 1] + gap_open_score,
            matrix_gap1[i - 1] + gap_extend_score)
        matrix_gap1[i] = max(options_gap1)
        if options_gap1[0] > options_gap1[1]:
            matrix_direction[i] |= A_HORZ_O  # Must exit M_GAP1 here.

        options_gap2 = (
            matrix_match[i - row_offset] + gap_open_score,
            matrix_gap2[i - row_offset] + gap_extend_score)
        matrix_gap2[i] = max(options_gap2)
        if options_gap2[0] > options_gap2[1]:
            matrix_direction[i] |= A_VERT_O  # Must exit M_GAP2 here.

        options = (
            matrix_match[i - 1 - row_offset] + match,
            matrix_gap1[i],
            matrix_gap2[i])
        matrix_match[i] = max(options)
        if options[0] == matrix_match[i]:
            matrix_direction[i] |= A_MATCH  # Can stay in M_MATCH here.
        if options[1] == matrix_match[i]:
            matrix_direction[i] |= A_HORZ_E  # Can enter M_GAP1 here.
        if options[2] == matrix_match[i]:
            matrix_direction[i] |= A_VERT_E  # Can enter M_GAP2 here.

    if debug:
        print("GAP1")
        for i in range(0, len(matrix_gap1), row_offset):
            print(("%5i" * row_offset) % tuple(matrix_gap1[i : i + row_offset]))
        print("GAP2")
        for i in range(0, len(matrix_gap2), row_offset):
            print(("%5i" * row_offset) % tuple(matrix_gap2[i : i + row_offset]))
        print("Match")
        for i in range(0, len(matrix_match), row_offset):
            print(("%5i" * row_offset) % tuple(matrix_match[i : i + row_offset]))
        print("FLAGS")
        for i in range(0, len(matrix_direction), row_offset):
            print(("%5s|" * row_offset) % tuple("".join((
                "h" if x & A_HORZ_O else " ",
                "H" if x & A_HORZ_E else " ",
                "D" if x & A_MATCH  else " ",
                "V" if x & A_VERT_E else " ",
                "v" if x & A_VERT_O else " "
            )) for x in matrix_direction[i : i + row_offset]))
        print("Traceback")


    # Backtracking.
    variants = []
    variant_template = 0
    variant_sequence = 0
    i = len(matrix_match) - 1
    in_matrix = M_MATCH  # May change before first step.
    while i >= 0:
        x = i % row_offset
        y = i // row_offset
        if debug:
            print("(%i, %i)" % (x, y))

        if in_matrix == M_MATCH:
            # Make gaps as soon as possible (pushed right).
            if matrix_direction[i] & A_HORZ_E:
                in_matrix = M_GAP1
            elif matrix_direction[i] & A_VERT_E:
                in_matrix = M_GAP2
            elif not (matrix_direction[i] & A_MATCH):
                raise ValueError(
                    "Alignment error: Dead route! (This is a bug.) [%s,%s]" % (template, sequence))

        if in_matrix == M_GAP1:
            # Go horizontally.  Deletion.
            variant_template += 1
            if matrix_direction[i] & A_HORZ_O:
                # End of gap, go diagonally after this.
                in_matrix = M_MATCH
            i -= 1
            continue

        if in_matrix == M_GAP2:
            # Go vertically.  Insertion.
            variant_sequence += 1
            if matrix_direction[i] & A_VERT_O:
                # End of gap, go diagonally after this.
                in_matrix = M_MATCH
            i -= row_offset
            continue

        # Go diagonally.  Either match or mismatch.
        if i != 0 and template[x - 1] != sequence[y - 1]:
            # Start/extend mismatch.
            variant_template += 1
            variant_sequence += 1

        else:
            # Match.  Flush variants.
            if variant_template or variant_sequence:
                variants.append(((x, variant_template), (y, variant_sequence)))
                variant_template = 0
                variant_sequence = 0
        i -= 1 + row_offset

    # Store the result in the cache.
    if cache:
        align_affine.cache[template, sequence, match_score, mismatch_score, gap_open_score, gap_extend_score] = variants
    return variants
#align_affine
align_affine.cache = {}


def name_variants(template, sequence, variants, *, location=(None, 1)):
    """
    Convert the list of variants (obtained through align_affine) to
    variant calls.  The format (nomenclature) of the
    returned variants depends on the location argument.

    If location is a tuple ("chromosome name", position), with any
    integer for the position, all variants are given as substitutions in
    the form posX>Y.  Insertions and deletions are written as pos.1->Y
    and posX>-, respectively.  The given position is that of the first
    base in the template.  With the location set to "suffix", a plus
    sign is prepended to position numbers and the first base in the
    template is pos=1.  With location set to "prefix", a minus sign is
    prepended and bases are counted from right to left instead.

    If location is a tuple ("M", position) with any integer for the
    position, variants are written following the mtDNA nomenclature
    guidelines.  The given position is that of the first base in the
    template.
    """
    variant_format = "%s%s>%s"
    if location == "prefix":
        location = ("prefix", -len(template))
    elif location == "suffix":
        # Include plus signs for position numbers.
        variant_format = "%+i%s>%s"
        location = ("suffix", 1)
    elif not isinstance(location, tuple) or len(location) < 2:
        raise ValueError("Unknown location %r. It should be 'prefix', "
            "'suffix', or a tuple (chromosome, position [, endpos])" % location)

    # Backtracking.
    named_variants = []
    for (pos_template, len_template), (pos_sequence, len_sequence) in sorted(variants):
        if location[0] == "M":
            # MtDNA variants are one-base-at-a-time.
            for j in range(max(len_template, len_sequence)-1, -1, -1):
                named_variants.append("%s%i%s%s" % (
                    template[pos_template + j] if j < len_template else "",#"-",
                    get_genome_pos(location, x + min(j, len_template - 1)),
                    ".%i" % (j - len_template + 1) if j >= len_template else "",
                    sequence[pos_sequence + j] if j < len_sequence else "DEL"))
        elif len_template == 0:
            # Insertions: "-130.1->C" instead of "-130->C".
            pos = get_genome_pos(location, pos_template - 1)
            if location[0] == "prefix":
                pos = "-%s" % abs(pos + 1)
            named_variants.append(variant_format % (pos, ".1-",
                sequence[pos_sequence : pos_sequence + len_sequence]))
        else:
            named_variants.append(variant_format % (
                get_genome_pos(location, pos_template),
                template[pos_template : pos_template + len_template],
                sequence[pos_sequence : pos_sequence + len_sequence] or "-"))

    if location[0] == "prefix":
        # Reverse variant order for prefix variants.
        named_variants.reverse()

    return named_variants
#name_variants


def call_variants(template, sequence, *, location=(None, 1), match_score=1, mismatch_score=-3, gap_open_score=-7, gap_extend_score=-2, cache=True, debug=False):
    variants = align_affine(template, sequence,
        match_score=match_score, mismatch_score=mismatch_score, gap_open_score=gap_open_score,
        gap_extend_score=gap_extend_score, cache=cache, debug=debug)
    return name_variants(template, sequence, variants, location=location)
#call_variants
