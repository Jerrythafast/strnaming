#
# Copyright (C) 2023 Jerry Hoogenboom
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
import itertools
import math
import re
import sys
import time

from functools import reduce

from .libsequence import align

# Maximum number of seconds that get_best_path() may spend.
MAX_SECONDS = 30
MAX_SECONDS_REFSEQ = 300

# Number of times a repeat must be repeated to make a significant enough repeat stretch for refseq analysis to kick in.
REFSEQ_MINIMUM_REPEATS = 4

# STRNaming settings.
# TODO: Split these up into separate globals.
# TODO: Put alignment settings here?
NAMING_OPTIONS = {
    "min_repeat_length": 8,
    "min_structure_length": 20,
    "max_unit_length": 6,
    "max_gap": 8,
    "max_long_gap": 20,
    "bases_covered_factor": 0.38033472,
    "units_used_factor": -8.8979406,
    "units_used_multiplier": 1.28993249,
    "repeats_factor": 2.12807513,
    "preferred_repeats_factor": 1.62639518,
    "interruptions_factor": -10.61080182,
    "interruptions_multiplier": 1.09212922,
    "nice_interruptions_factor": 7.23604507,
    "interruption_bases_factor": -0.46646098,
    "prefix_suffix_factor": -1.37392421
}


class OutOfTimeException(Exception):
    pass


class ComplexityException(Exception):
    pass


def detect_repeats(seq):
    """
    Return a matrix of max_unit_length x len(seq) elements.
    The value at position (j,i) gives the length of a repeat of units of
    length j (in number of units), ending in position i in the sequence.
    """
    len_seq = len(seq)
    matrix = [[0]*len_seq for j in range(NAMING_OPTIONS["max_unit_length"])]
    for j in range(NAMING_OPTIONS["max_unit_length"]):
        row = matrix[j]
        for i in range(len_seq):
            if i <= j:
                row[i] = i + 1  # No full unit yet.
            elif seq[i] != seq[i - j - 1]:
                row[i] = j + 1  # Restart counting.
            else:
                row[i] = row[i - 1] + 1
    for j in range(NAMING_OPTIONS["max_unit_length"]):
        row = matrix[j]
        for i in range(len_seq):
            row[i] //= (j + 1)
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
    if length < NAMING_OPTIONS["min_repeat_length"]:
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


def find_repeated_units(seq):
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
    return units
#find_repeated_units


def find_block_length(path, len_or_int=int):
    return max((end - start, len_or_int(unit)) for start, end, unit, *_ in path)[1]
#find_block_length


def find_overlong_gap(path):
    """Return the longest interruption, if that interruption is longer than max_gap."""
    longest_gap = NAMING_OPTIONS["max_gap"]
    longest_start_pos = 0
    for i in range(1, len(path)):
        start_pos = path[i - 1][1]
        gap_length = path[i][0] - start_pos
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


def calc_path_score(len_prefix, len_suffix, len_seq, ref_len_large_gap, path, block_length, preferred_units=None,
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
    large_gap_delta = ref_len_large_gap - max(interruptions, default=0) if ref_len_large_gap else 0
    #if not ref_len_large_gap and max(interruptions, default=0) > 8:
    #    # Idea for future version: Introducing a large gap while the refseq had none counts as an additional gap.
    #    # This greatly reduces (>90%) the number of reference structures with large gaps.
    #    raise ValueError("Did not want to introduce this yet")
    #    interruptions.append(0)

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


def unique_set(all_sets, this_set: frozenset):
    """
    Memory-saving helper function. Return an equivalent set from all_sets
    for this_set, or add this_set to all_sets and return this_set.
    Note that this_set must be a frozenset.
    """
    return all_sets.setdefault(this_set, this_set)
#unique_set


def generate_scaffolds(repeat_list, scaffold=None, anchored=None, orphans=None, start=0, sets=None):
    if sets is None:
        sets = {}
    if scaffold is None:
        empty_set = unique_set(sets, frozenset())
        for next_i, repeat in enumerate(itertools.islice(repeat_list, start, None), start + 1):
            unit_set = unique_set(sets, frozenset({repeat[2]}))
            anchor = repeat[3][2]
            new_scaffold = [repeat]
            new_anchored = unit_set if anchor else empty_set
            new_orphans = empty_set if anchor else unit_set
            yield (new_scaffold, new_anchored, new_orphans)
            yield from generate_scaffolds(repeat_list, new_scaffold, new_anchored, new_orphans, next_i, sets)
        return

    scaffold_end_pos = scaffold[-1][1]
    scaffold_end_unit = scaffold[-1][2]
    units_limit_reached = len(anchored | orphans) >= 6
    for next_i, repeat in enumerate(itertools.islice(repeat_list, start, None), start + 1):
        repeat_start_pos = repeat[0]
        if repeat_start_pos < scaffold_end_pos:
            continue
        if repeat_start_pos > scaffold_end_pos:
            break  # The repeat_list is sorted, so we're done here.

        unit = repeat[2]
        if unit == scaffold_end_unit:
            # We're not putting repeats of the same unit adjacently. (Those
            # have been split up by upstream code to allow more junctions.)
            continue

        unit_is_anchored = unit in anchored
        unit_is_orphan = unit in orphans
        if units_limit_reached and not unit_is_anchored and not unit_is_orphan:
            continue  # Limit the number of different units in one scaffold.

        new_scaffold = scaffold + [repeat]
        new_anchored = anchored
        new_orphans = orphans
        if repeat[3][2]:
            if not unit_is_anchored:
                new_anchored = unique_set(sets, new_anchored | {unit})
                if unit_is_orphan:
                    new_orphans = unique_set(sets, new_orphans - {unit})
        elif not unit_is_anchored and not unit_is_orphan:
            new_orphans = unique_set(sets, new_orphans | {unit})
        yield (new_scaffold, new_anchored, new_orphans)
        yield from generate_scaffolds(repeat_list, new_scaffold, new_anchored, new_orphans, next_i, sets)
#generate_scaffolds


def get_scaffolds(repeat_list):
    scaffolds = {}
    for n, scaffold in enumerate(generate_scaffolds(repeat_list)):
        if n == 5000000:
            raise ComplexityException("%r" % repeat_list)  # More than 5 million scaffolds?!
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


def get_gaps(repeat_list, scaffolds, seq):
    short_gaps = {}
    all_gaps = {}
    for _, start, unit1, _ in repeat_list:
        for end, _, unit2, _ in repeat_list:
            size = end - start
            if size <= 0 or (start in scaffolds and end in scaffolds[start]
                    and any(not scaffold[2] for scaffold in scaffolds[start][end])):
                continue  # Illegal gap: can be filled with scaffold without orphans.
            if size > NAMING_OPTIONS["max_long_gap"]:
                continue  # Gap is too large.
            gapseq = seq[start : end]
            if size == len(unit1) and gapseq == unit1:
                continue  # Illegal gap: same as previous unit.
            if size == len(unit2) and gapseq == unit2:
                continue  # Illegal gap: same as next unit.
            try:
                all_gaps[start][end] = gapseq
            except KeyError:
                all_gaps[start] = {end: gapseq}
            if size <= NAMING_OPTIONS["max_gap"]:
                try:
                    short_gaps[start][end] = gapseq
                except KeyError:
                    short_gaps[start] = {end: gapseq}
    return short_gaps, all_gaps
#get_gaps


def extend_with_gaps(extended, gap_ranges, ranges_before, ranges_after):
    for start, ends in ranges_before.items():
        subtree = extended[start]
        for gap_start, (anchorables, orphans) in ends.items():
            if gap_start in gap_ranges:
                for gap_end in gap_ranges[gap_start]:
                    for final_end, (more_anchorables, more_orphans) in ranges_after[gap_end].items():
                        must_use_units = orphans | more_orphans
                        if len(must_use_units) > 6:
                            # This combination of before-range and after-range is not valid,
                            # because it uses too many different units. Note that this filter is
                            # not perfect, as we only know the number of must-orphan units.
                            continue
                        if final_end in subtree:
                            existing_anchorables, existing_orphans = subtree[final_end]
                            new_anchorables = existing_anchorables | anchorables | more_anchorables
                            new_orphans = (must_use_units & existing_orphans) - new_anchorables
                            subtree[final_end] = (new_anchorables, new_orphans)
                        else:
                            new_anchorables = anchorables | more_anchorables
                            new_orphans = must_use_units - new_anchorables
                            subtree[final_end] = (new_anchorables, new_orphans)
#extend_with_gaps


def extend_ranges_with_gaps(gap_ranges, *prior_ranges):
    # The prior_ranges should be sorted by the number of gaps allowed, ascending.
    extended = {start: ends.copy() for start, ends in prior_ranges[-1].items()}
    for ranges_before, ranges_after in zip(prior_ranges, prior_ranges[::-1]):
        # [0] -> [(0,0)]    adding first gap: allow no gaps on either side
        # [0,1] -> [(0,1),(1,0)]    adding second gap: allow 1 gap on 1 side only
        # [0,1,2] -> [(0,2),(1,1),(2,0)]    adding third gap: allow 2 gaps on one side only, or 1 gap on both sides
        # [0,1,2,3] -> [(0,3),(1,2),(2,1),(3,0)]    adding fourth gap: allow 3 gaps on one side only, or 2 gaps on one side with 1 on the other
        # [0,1,2,3,4] -> [(0,4),(1,3),(2,2),(3,1),(4,0)]    adding fifth gap: allow 4 gaps on one side only, 3 gaps on one side with 1 on the other, or 2 gaps on both sides
        extend_with_gaps(extended, gap_ranges, ranges_before, ranges_after)
    return extended
#extend_ranges_with_gaps


def get_ranges(scaffolds, short_gaps, all_gaps, is_refseq):
    # TODO: This function is now hardcoded for MAX_GAPS=5. Maybe generalise?
    # NOTE: this function may report ranges that are not acutally possible due to
    # rules regarding e.g. singletons.
    # Lookup tables for starting positions and associated ending positions.
    # The positions are linked with the units that are anchorable along the way.
    start_end00 = {start:
            {end: (
                reduce(frozenset.union, (scaffold[1] for scaffold in scaffolds_here)),
                reduce(frozenset.intersection, (scaffold[2] for scaffold in scaffolds_here)))
            for end, scaffolds_here in ends.items()}
        for start, ends in scaffolds.items()}
    start_end10 = extend_ranges_with_gaps(short_gaps, start_end00)
    start_end20 = extend_ranges_with_gaps(short_gaps, start_end00, start_end10)
    start_end30 = extend_ranges_with_gaps(short_gaps, start_end00, start_end10, start_end20)
    start_end40 = extend_ranges_with_gaps(short_gaps, start_end00, start_end10, start_end20, start_end30)
    if is_refseq:
        start_end50 = extend_ranges_with_gaps(short_gaps, start_end00, start_end10, start_end20, start_end30, start_end40)
    else:
        start_end11 = extend_ranges_with_gaps(all_gaps, start_end00)
        start_end21 = extend_ranges_with_gaps(all_gaps, start_end00, start_end10)
        start_end31 = extend_ranges_with_gaps(all_gaps, start_end00, start_end10, start_end20)
        start_end41 = extend_ranges_with_gaps(all_gaps, start_end00, start_end10, start_end20, start_end30)
        start_end51 = extend_ranges_with_gaps(all_gaps, start_end00, start_end10, start_end20, start_end30, start_end40)

    # Now, build a tree from these.
    node00 = (start_end00, {}, None, None, start_end00)
    node10 = (start_end10, short_gaps, node00, None, start_end00)
    node20 = (start_end20, short_gaps, node10, None, start_end00)
    node30 = (start_end30, short_gaps, node20, None, start_end00)
    node40 = (start_end40, short_gaps, node30, None, start_end00)
    if is_refseq:
        return (start_end50, short_gaps, node40, None, start_end00)
    node11 = (start_end11, all_gaps, node00, node00, start_end00)
    node21 = (start_end21, all_gaps, node11, node10, start_end00)
    node31 = (start_end31, all_gaps, node21, node20, start_end00)
    node41 = (start_end41, all_gaps, node31, node30, start_end00)
    node51 = (start_end51, all_gaps, node41, node40, start_end00)
    return node51
#get_ranges


def recurse_gen_path(start_pos, end_pos, scaffolds, ranges, endtime,
                     anchored=frozenset(), orphans=frozenset(), prev_gap=""):
    # Timekeeping.
    if time.monotonic() > endtime:
        raise OutOfTimeException()
    scaffolds_here = scaffolds[start_pos]
    try:
        units_so_far = anchored | orphans
        for scaffold, anchored2, orphans2 in scaffolds_here[end_pos]:
            if prev_gap == scaffold[0][2]:
                continue  # Scaffold must not start with unit equal to previous gap.
            if (not prev_gap and not scaffold[0][3][3]) or not scaffold[-1][3][3]:
                continue  # Can't start or end with a non-ref repeat unit.
            if orphans - anchored2 or orphans2 - anchored:
                continue  # Would end up with orphans.
            if len(units_so_far | anchored2 | orphans2) > 6:
                continue  # Would end up with too many different units.
            yield scaffold
    except KeyError:
        pass  # No scaffolds reach exactly from start_pos to end_pos.

    # Insert a scaffold that reaches not quite to end_pos just yet.
    gaps = ranges[1]
    ranges_without_gaps_here = ranges[-1][start_pos]
    for gap_pos, scaffold_list in scaffolds_here.items():

        if gap_pos not in gaps:
            continue

        anchorable_before_gap, orphaned_before_gap = ranges_without_gaps_here[gap_pos]
        anchorable_before = anchored | anchorable_before_gap  # Possibly.
        orphaned_before = orphans | orphaned_before_gap  # Certainly.

        # Choose the next gap and continue the structure on the other side.
        for next_pos, gapseq in gaps[gap_pos].items():
            gap_len = next_pos - gap_pos
            remaining_ranges = ranges[2 if gap_len <= NAMING_OPTIONS["max_gap"] else 3]

            try:
                anchorable_after_gap, orphaned_after_gap = remaining_ranges[0][next_pos][end_pos]
            except KeyError:
                continue  # Can't reach from next_pos to end_pos.

            if (orphaned_before | orphaned_after_gap) - (anchorable_before | anchorable_after_gap):
                continue  # Using this gap implies we will end up with orphans.

            # The following check is meant to improve performance, but it might take more time than it saves.
            #if len(anchored | orphaned_before | orphaned_after_gap) > 6:
            #    continue  # Using this gap implies we will be using too many different units.

            for scaffold, anchored2, orphans2 in scaffold_list:
                if prev_gap == scaffold[0][2]:
                    continue  # Scaffold must not start with unit equal to previous gap.
                if gapseq == scaffold[-1][2]:
                    continue  # Scaffold must not end with unit equal to the new gap.
                if not prev_gap and not scaffold[0][3][3]:
                    continue  # Can't start with this non-ref repeat unit.
                now_anchored = anchored | anchored2
                if orphaned_after_gap - now_anchored:
                    continue  # Unavoidable orphans at the other side of the gap.
                now_orphaned = (orphans | orphans2) - now_anchored
                if now_orphaned - anchorable_after_gap:
                    continue  # Can't anchor all our orphans at the other side of the gap.
                if len(now_anchored | now_orphaned) > 6:
                    continue  # Would end up with too many different units.
                for substructure in recurse_gen_path(next_pos, end_pos,
                        scaffolds, remaining_ranges, endtime,
                        now_anchored, now_orphaned, gapseq):
                    yield scaffold + substructure
#recurse_gen_path


def gen_all_paths(prefix, suffix, seq, repeat_list, is_refseq, endtime):
    if not is_refseq:
        # Filter out repeat stretches that fall completely in the prefix or suffix.
        if seq.startswith(prefix):
            len_prefix = len(prefix)
        else:
            zero, len_prefix, score = align(seq, prefix, prefix=True)
            len_prefix = len_prefix if score > 0 else None
        if len_prefix:
            repeat_list = [repeat for repeat in repeat_list if repeat[1] > len_prefix]
        if seq.endswith(suffix):
            pos_suffix = len(seq) - len(suffix)
        else:
            zero, len_suffix, score = align(seq[-1::-1], suffix[-1::-1], prefix=True)
            pos_suffix = len(seq) - len_suffix if score > 0 else None
        if pos_suffix:
            repeat_list = [repeat for repeat in repeat_list if repeat[0] < pos_suffix]

    scaffolds = get_scaffolds(repeat_list)
    short_gaps, all_gaps = get_gaps(repeat_list, scaffolds, seq)
    ranges = get_ranges(scaffolds, short_gaps, all_gaps, is_refseq)

    # A *preferred* repeat unit is required as the starter/finisher.
    start_positions = set(x[0] for x in repeat_list if x[3][3])
    end_positions = set(x[1] for x in repeat_list if x[3][3])

    if not is_refseq:
        # Try anchoring the prefix and suffix.
        success = False
        has_prefix = len_prefix in start_positions
        has_suffix = pos_suffix in end_positions

        # First, try to reach exactly from prefix to suffix.
        if has_prefix and has_suffix and pos_suffix in ranges[0][len_prefix] and not ranges[0][len_prefix][pos_suffix][1]:
            for result in recurse_gen_path(len_prefix, pos_suffix, scaffolds, ranges, endtime):
                success = True
                yield result
            if success:
                return  # Yielded path with anchored prefix AND suffix.
            # Turns out that pos_suffix is not reachable from len_prefix.
            del ranges[0][len_prefix][pos_suffix]

        # If that did not work, force starting at the end of prefix
        # or finishing at the start of suffix separately.
        if has_prefix:
            for end_pos in end_positions & ranges[0][len_prefix].keys():
                if not ranges[0][len_prefix][end_pos][1]:
                    for result in recurse_gen_path(len_prefix, end_pos, scaffolds, ranges, endtime):
                        success = True
                        yield result
            if not success:
                # Turns out we can't start from len_prefix.
                start_positions.remove(len_prefix)
        if has_suffix:
            for start_pos in start_positions:
                if pos_suffix in ranges[0][start_pos] and not ranges[0][start_pos][pos_suffix][1]:
                    for result in recurse_gen_path(start_pos, pos_suffix, scaffolds, ranges, endtime):
                        success = True
                        yield result
            if not success:
                # Turns out we can't end at pos_suffix.
                end_positions.remove(pos_suffix)
        if not success:
            # No path with anchored prefix OR suffix found; try everything else.
            for start_pos in start_positions:
                for end_pos in end_positions & ranges[0][start_pos].keys():
                    if not ranges[0][start_pos][end_pos][1]:
                        yield from recurse_gen_path(start_pos, end_pos, scaffolds, ranges, endtime)

    else:
        # Refseq structures must include a significant repeat of at least 8 nt.
        significant_repeats = [r for r in repeat_list if r[3][0] >= max(NAMING_OPTIONS["min_repeat_length"], r[3][1] * REFSEQ_MINIMUM_REPEATS)]
        minimal_end_pos_table = {}
        minimal_end_pos = len(seq)
        for repeat in reversed(significant_repeats):
            minimal_end_pos = min(minimal_end_pos, repeat[1])
            minimal_end_pos_table[repeat[0]] = minimal_end_pos
        significant_repeats = iter(significant_repeats)
        significant_repeat = next(significant_repeats, None)
        for start_pos in sorted(start_positions):
            while significant_repeat is not None and significant_repeat[0] < start_pos:
                significant_repeat = next(significant_repeats, None)
            if significant_repeat is None:
                # Start position comes after start of last significant repeat.
                return
            minimal_end_pos = max(
                start_pos + NAMING_OPTIONS["min_structure_length"],
                minimal_end_pos_table[significant_repeat[0]])
            for end_pos in end_positions & ranges[0][start_pos].keys():
                if end_pos >= minimal_end_pos and not ranges[0][start_pos][end_pos][1]:
                    yield from recurse_gen_path(start_pos, end_pos, scaffolds, ranges, endtime)
#gen_all_paths


def get_best_path(prefix, suffix, seq, repeat_list, preferred_units, ref_block_length, ref_len_large_gap, offset=0):
    """Return combination of repeats that maximises calc_path_score."""
    len_prefix = len(prefix) if prefix is not None else -1
    len_suffix = len(suffix) if suffix is not None else -1
    len_seq = len(seq)
    is_refseq = ref_block_length == 0
    best_score = -sys.maxsize
    best_path = []
    endtime = time.monotonic() + (MAX_SECONDS_REFSEQ if is_refseq else MAX_SECONDS)

    for path in gen_all_paths(prefix, suffix, seq, repeat_list, is_refseq, endtime):
        block_length = ref_block_length or find_block_length(path, len_or_int=len)

        # Calculate the score of the current path.
        score = calc_path_score(len_prefix, len_suffix, len_seq, ref_len_large_gap, path,
            block_length, preferred_units, is_refseq)

        if score >= best_score:
            if score == best_score:
                # Tiebreak #1: prefer fewer stretches.
                len_diff = len(best_path) - len(path)
                better_path = len_diff > 0
                if not better_path and not len_diff:
                    # Tiebreak #2: Prefer most 5' starting position.
                    for a, b in zip(best_path, path):
                        if a[0] == b[0]:
                            continue
                        better_path = a[0] > b[0]
                        break
                if not better_path:
                    continue
            best_score = score
            best_path = path
    return best_score, [[repeat[0] + offset, repeat[1] + offset, repeat[2]] for repeat in best_path]
#get_best_path


def find_repeat_stretches(seq, units, allow_one, preferred_units, repeat_list=None):
    """Return a list of repeat stretches. The input units must be sorted by decreasing length."""
    # Find all stretches of each repeated unit.
    if repeat_list is None:
        repeat_list = []
    for unit in units:

        # Compile a list of all places where this unit is found.
        # len(stretch) must be >= 3.
        # NOTE: Hardcoded minimum length; singletons of 2nt or shorter are not included.
        regex = re.compile("(" + unit + "){" + str(3 if len(unit) == 1 else 2 if len(unit) == 2 or not allow_one else 1) + ",}")

        # Ignore repeats that completely overlap a repeat of a longer unit.
        # Using those is not expected to give a high score, but their
        # presence may greatly increase calculation time.
        # Due to the way this function is called, preferred repeats will not be ignored
        # if they overlap with unpreferred repeats. This is the desired behaviour.
        # Example markers where this filter has an effect: DYF387S1, PentaD, DYS448.
        ignore_ranges = [(o_start, o_end) for o_start, o_end, o_unit in repeat_list
            if len(unit) < len(o_unit)  # Other is longer unit (will be 4+)
            and o_end - o_start >= 8]  # Other is repeated

        pos = 0
        while True:
            match = regex.search(seq[pos:])
            if match is None:
                break
            match_start = pos + match.start()
            match_end = match_start + len(match.group(0))
            pos = match_end - len(unit) + 1
            if not any(match_start >= o_start and match_end <= o_end for o_start, o_end in ignore_ranges):
                repeat_list.append([match_start, match_end, unit])

    # Sort repeats.
    repeat_list.sort()

    # Remove singletons that completely overlap an 8nt+ repeat of a preferred unit.
    significant_ranges = [(o_start, o_end)
        for o_start, o_end, o_unit in repeat_list
        if o_unit in preferred_units
        and o_end - o_start >= 8]
    i = 0
    while i < len(repeat_list):
        start, end, unit = repeat_list[i]
        if (len(unit) == end - start and  # Current is singleton
                any(start >= o_start and end <= o_end  # Overlap
                    for o_start, o_end in significant_ranges)):
            # The current (singleton) unit overlaps completely with a significant repeat.
            del repeat_list[i]
        else:
            i += 1

    return repeat_list
#find_repeat_stretches


def trim_overlapping_repeats(repeat_list, preferred_units):
    """
    Where repeat stretches overlap, also include trimmed copies that don't overlap.
    """
    repeat_trims = tuple(({0}, {0}) for x in range(len(repeat_list)))
    for i in range(len(repeat_list)):
        for j in range(i + 1, len(repeat_list)):
            overlap_len = repeat_list[i][1] - repeat_list[j][0]
            if overlap_len <= 0:
                break  # The repeats are sorted, so subsequent j will not overlap.
            # Stretch i extends past the start of stretch j.
            # We can trim the end of stretch i, or the start of stretch j, to compensate.
            if repeat_list[i][2] in preferred_units and repeat_list[j][2] in preferred_units:
                # Always trimming if both repeats involve preferred units.
                repeat_trims[i][1].add(math.ceil(float(overlap_len)/len(repeat_list[i][2]))*len(repeat_list[i][2]))
                repeat_trims[j][0].add(math.ceil(float(overlap_len)/len(repeat_list[j][2]))*len(repeat_list[j][2]))
            else:
                # Otherwise, only trimming if doing so leaves no gap.
                if overlap_len % len(repeat_list[i][2]) == 0:
                    repeat_trims[i][1].add(overlap_len)
                if overlap_len % len(repeat_list[j][2]) == 0:
                    repeat_trims[j][0].add(overlap_len)
    for i in range(len(repeat_trims)):
        start, end, unit = repeat_list[i]
        # Don't trim unpreferred units such that they violate min_repeat_length.
        # For ref units, a lower threshold of 3 nt is applied.
        # NOTE: hardcoded number.
        min_length = 3 if unit in preferred_units else NAMING_OPTIONS["min_repeat_length"]
        for trim_start in repeat_trims[i][0]:
            for trim_end in repeat_trims[i][1]:
                trim_length = trim_start + trim_end
                if trim_length and end - start - trim_length >= min_length:
                    repeat_list.append([start + trim_start, end - trim_end, unit])
    repeat_list.sort()

    # Deduplicate repeats.
    i = 1
    while i < len(repeat_list):
        if tuple(repeat_list[i-1]) == tuple(repeat_list[i]):
            del repeat_list[i]
        else:
            i += 1
#trim_overlapping_repeats


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
    repeat_list = find_repeat_stretches(seq, preferred_units, True, preferred_units)
    for i, (start, end, preferred_unit) in enumerate(repeat_list):
        if start:
            # Find units before repeat i.
            locations = find_repeated_units(seq[:start])
            for unit in locations:
                try:
                    unit_locations[unit] += locations[unit]
                except KeyError:
                    unit_locations[unit] = locations[unit]
        if len(seq) - end:
            # Find units after repeat i.
            locations = find_repeated_units(seq[end:])
            for unit in locations:
                new_locations = [(x[0] + end, x[1]) for x in locations[unit]]
                try:
                    unit_locations[unit] += new_locations
                except KeyError:
                    unit_locations[unit] = new_locations
        if i:
            # Find units between repeat i-1 and repeat i.
            prev_end = repeat_list[i-1][1]
            if prev_end < start:
                locations = find_repeated_units(seq[prev_end : start])
                for unit in locations:
                    new_locations = [(x[0] + prev_end, x[1]) for x in locations[unit]]
                    try:
                        unit_locations[unit] += new_locations
                    except KeyError:
                        unit_locations[unit] = new_locations

    # If no stretches of preferred units were found, just find any stretches.
    if not repeat_list:
        locations = find_repeated_units(seq)
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
        False, preferred_units, repeat_list)

    # Where repeat stretches overlap, also include trimmed copies that don't overlap.
    trim_overlapping_repeats(repeat_list, preferred_units)

    # Add (repeat_length, unit_length, anchor, preferred) to the end.
    for repeat in repeat_list:
        start, end, unit = repeat
        anchor = any(start < start_l + length and end > start_l
            for start_l, length in unit_locations[unit])
        repeat.append((end - start, len(unit), anchor, unit in preferred_units))

    # Return what we have found.
    return repeat_list
#find_everything


### Entry point for target analysis ###


def collapse_repeat_units(seq, prefix, suffix, preferred_units, ref_block_length, ref_len_large_gap):
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
    repeat_list = find_everything(seq, unit_locations)

    # Find best combination of repeat unit usage.
    score, path = get_best_path(prefix, suffix, seq, repeat_list, preferred_units, ref_block_length, ref_len_large_gap)
    return path
#collapse_repeat_units


"""
Differences Naming after Refseq:
* Can start to introduce one gap of more than 8 nt (up to 20 nt) without inflicting the prefix_suffix_factor penalty
* Block length no longer varies per path, so nice gap points may change for previously-not-optimal paths
* Preferred-unit repeats < 8 nt become anchoring, so can now choose to use only < 8 nt repeats of a unit
* Some (or many) units may become unpreferred, making repeat trimming and filtering much less aggressive


Refseq Rules:

collapse_repeat_units_refseq:
  * 1. Start search around the longest repeat, cancel if not 8+ nt AND 4+ repeats
  * 2. Extend region to include nearby repeats of 8+ nt, if they reach within 20 nt of current region
  * 3. Determine set of repeat units that have a repeat of 8+ nt within the region
  * 4. Extend region to include all nearby instances of those units; stop before the first gap strictly larger than 8 nt
  * 5. Loop back to 2 if more 8+ nt repeats come within 20 nt of the current region
  * 6. Call find_repeat_stretches: all units are preferred and singletons are allowed
  * 7. Call trim_overlapping_repeats: all overlapping repeats are trimmed, as long as the trimmed repeat is at least 3 nt long
  * 8. All repeats are tagged as preferred; all repeats of 8+ nt are tagged as anchoring
  * 9. After get_best_path, only units that are used in non-singleton positions are tagged as preferred
  * 10. Perform regular allele naming (find_everything + get_best_path), looping back to 9 until the same name is obtained consecutively


Naming Rules:

collapse_repeat_units:
  * Units used in the refseq name are the preferred units
  * Just call find_everything and then get_best_path to obtain the best path


Repeat Finding Rules:

find_everything:
  * Call find_repeat_stretches with the preferred units, allowing singletons
  * If any repeat stretches of preferred units are found, call find_repeated_units on the follwing ranges:
    * All sequence before the start of a repeat stretch of a preferred unit (for each repeat stretch)
    * All sequence after the end of a repeat stretch of a preferred unit (for each repeat stretch)
    * All sequence between two consecutive repeat stretches of preferred units (for each pair of consecutive repeat stretches)
  * Otherwise, call find_repeated_units on the entire sequence once
  * Call find_repeat_stretches with the non-preferred units, disallowing singletons, adding to the list of preferred repeats
  * Call trim_overlapping_repeats
  * Mark as anchoring:
    * All repeats (including singletons) of preferred units
    * Non-preferred repeats that overlap a location (even by 1 nt) for that unit as returned by any find_repeated_units call

find_repeat_stretches:
  * Minimum stretch length is 3 nt
  * Singletons are only allowed for preferred repeat units
  * Ignore repeats that completely overlap a 8nt+ repeat of a longer unit (this requires sorted input)
  * Remove singletons that completely overlap a 8+ nt stretch of a preferred unit

trim_overlapping_repeats:
  * Add trimmed copies of repeat stretches with a minimum length of 3 nt for preferred units and 8 nt for non-preferred units
  * Where two preferred repeats overlap, add trimmed copies of both
  * If one (or both) of the repeats is non-preferred, only add trimmed copies if doing so leaves no gap

find_repeated_units:
  * Recursively find the longest non-overlapping repeat and return its starting position and length
  * The cycle_units option is no longer used  --> removed

recurse_find_longest_repeat:
  * Report only repeat stretches of 8+ nt

find_longest_repeat:
  * Break ties by choosing the shortest unit, the earliest in the sequence


Best Path Rules:

calc_path_score
  * For non-refseq, preferred units always count as used

get_best_path
  * Choose the path with the highest score
  * In case of a tie, prefer fewer stretches if possible
  * Prefer the most 5' starting position(s) otherwise
  * For non-refseq, block length is fixed to refseq's block length

gen_all_paths
  * Start and end with a repeat marked as *preferred*
  * For non-refseq:
    * Ignore repeat stretches completely in prefix or suffix
    * Only yield paths that reach exactly from prefix to suffix if possible
    * Else, yield all paths that start at prefix or end at suffix
    * If none of that is possible, just yield all possible paths
  * For refseq:
    * Yield only paths that span a stretch of REFSEQ_MINIMUM_REPEATS=4 repeats and 8nt or longer
    * Yield only paths that are at least min_structure_length=20 nt

recurse_gen_path
  * Require that paths start and end with a repeat marked as *preferred*
  * Require that paths leave no orphans
  * Don't make gaps with sequence equal to an adjacent repeat unit

get_ranges
  * Make at most 5 gaps
  * Don't make large gaps (>8 nt) in reference sequences

get_gaps
  * Don't suggest gaps that can be filled with a scaffold without orphans
  * Don't suggest gaps that have the same sequence as adjacent units

generate_scaffolds
  * Don't put repeats of the same unit adjacently
"""