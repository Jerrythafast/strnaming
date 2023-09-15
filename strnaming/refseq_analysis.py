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
import re

from .libstrnaming import MANY_TIMES, NAMING_OPTIONS, ComplexityException, OutOfTimeException, \
                          get_best_path, find_block_length, find_everything, find_overlong_gap, \
                          find_repeat_stretches, trim_overlapping_repeats

#import numpy as np  # Only imported when actually running this code.

import json #FIXME, temp


class OscillationException(Exception):
    pass


def detect_repeats_np(seq):
    """
    Return a matrix of max_unit_length x len(seq) elements.
    The value at position (j,i) gives the length of a repeat of units of
    length j (in number of units), ending in position i in the sequence.
    """
    sa = np.fromiter(map(ord, seq), count=len(seq), dtype=np.uint8)
    matrix = np.empty((NAMING_OPTIONS["max_unit_length"], len(seq)), dtype=np.int64)
    for j in range(NAMING_OPTIONS["max_unit_length"]):
        # Inspired by https://stackoverflow.com/a/42129610
        jp1 = j + 1
        a_ext = np.concatenate(([0], [0]*jp1, sa[:-jp1] == sa[jp1:], [0]))
        idx = np.flatnonzero(a_ext[1:] != a_ext[:-1])
        a_ext[1:][idx[1::2]] = idx[::2] - idx[1::2]
        a_ext[jp1] = jp1
        a_ext[1:-1].cumsum(out=matrix[j])
        matrix[j] //= jp1
    return matrix
#detect_repeats_np


def find_longest_repeat_np(matrix):
    """
    Get start and end position and unit length of the longest repeat in
    matrix.  Break ties by choosing the shortest unit, the earliest in
    the sequence (lowest indices into matrix).
    """
    max_repeats = matrix.max(1)
    max_lengths = np.where(max_repeats == 1, 0, max_repeats * np.arange(1, 7))
    max_unit = max_lengths.argmax()
    max_len = max_lengths[max_unit]
    max_end = matrix[max_unit,:].argmax() + 1
    max_start = max_end - max_len
    return max_start, max_end, max_unit + 1
#find_longest_repeat_np


def grow_significant_repeats(matrix, locations, start, end, max_distance=20):
    """Reduce start and increase end by including significant repeats within max_distance nt."""
    while True:
        # Try to extend the scope to the left.
        pos_offset = max(0, start - max_distance)
        significant_end_positions = np.flatnonzero(locations[pos_offset : start])
        if not significant_end_positions.size:
            break
        # Found at least one significant repeat that ends within max_distance nt of current start.
        # Adjust start position to the start of that stretch.
        pos = significant_end_positions[0] + pos_offset  # Inclusive end position.
        ulen = locations[pos]
        start = pos + 1 - matrix[ulen - 1, pos] * ulen
    while True:
        # Try to extend the scope to the right.
        # We have to look further to the right, because long repeats are marked only on their ending position.
        sublocations = locations[end : end + max_distance + 47]  # 8*6 - 1; if sublocations[-1]==6, a 8x6nt repeat starts at the max_distance base after current scope end.
        # Mask out repeats that are potentially too far to the right: only want repeats that certainly start within max_distance nt of the scope NOW.
        # If they should be included (for being close to a previous significant repeat), they will get included in the next loop.
        significant_end_positions = np.flatnonzero(sublocations > (np.maximum(np.arange(len(sublocations)) - (max_distance-1), 0) // 8))
        if not significant_end_positions.size:
            break
        # Found at least one significant repeat that starts within max_distance nt of current end.
        # Adjust end position to the end of that stretch.
        end += significant_end_positions[-1] + 1  # +1 to convert from inclusive to exclusive position.
    return int(start), int(end)  # Convert from np.int64 to native Python int
#grow_significant_repeats


def find_repetitive_end(seq, units, pos=0):
    """Return furthest position in seq that can be reached with repeats when leaving no gap > 8 nt."""
    regexes = {unit: re.compile("(" + unit + "){" + str(3 if len(unit) == 1 else 2 if len(unit) == 2 else 1) + ",}") for unit in units}
    prev_pos = pos - 1
    while pos > prev_pos:
        prev_pos = pos
        for unit, regex in regexes.items():
            # Start looking slightly before current pos, to allow moving forward by less than one repeat.
            search_start_pos = max(0, pos - (2 if len(unit) == 1 else 3 if len(unit) == 2 else len(unit) - 1))
            match = regex.search(seq, search_start_pos, search_start_pos + 50)  # Trade-off loops vs excessive regex scan.
            if match and match.start() - pos <= 8:
                pos = match.end()
    return pos
#find_repetitive_end


def grow_scope(seq, matrix, locations, start, end):
    """
    Look for repeats in seq around seq[start:end] and return the found start, end, units.

    New repeat units must have an 8nt+ repeat within 20 nt of any units found before.
    The scope is extended to include all occurences of found units, leaving no gaps of >8 nt.
    These two steps are repeated until the algorithm converges.
    """

    # Repeatedly find any repeats of 8+ nt that are not far away from repeats found earlier.
    start, end = grow_significant_repeats(matrix, locations, start, end)

    while True:

        # Determine the set of repeat units for which we have found a repeat of 8+ nt.
        units = set(seq[pos-ulen+1 : pos+1] for pos, ulen in enumerate(locations[start:end], start) if ulen)

        # Find out where the furthest repeats of these units end, leaving no gaps of >8 nt.
        # Look in the 5' direction first, by reversing all sequences (and indexing).
        start = len(seq) - find_repetitive_end(seq[::-1], {unit[::-1] for unit in units}, len(seq) - start)
        end = find_repetitive_end(seq, units, end)

        # Double-check that there are no additional significant repeat stretches within the scope now.
        new_start, new_end = grow_significant_repeats(matrix, locations, start, end)
        if new_start < start or new_end > end:
            # New significant repeat stretches entered the scope, loop again.
            start = new_start
            end = new_end
        else:
            return start, end, units
#grow_scope


def collapse_repeat_units_refseq(seq, *, offset=0):
    if len(seq) < NAMING_OPTIONS["min_structure_length"]:
        return 0, []

    # Import numpy now.
    global np
    import numpy as np

    # First, find the longest repeat.
    matrix = detect_repeats_np(seq)
    start, end, u_len = find_longest_repeat_np(matrix)

    # If the longest repeat is not significant (8+ nt and 4+ repeats), stop.
    if end - start < max(NAMING_OPTIONS["min_repeat_length"], u_len * MANY_TIMES):
        return 0, []

    # Find all stretches of at least 8 nt.
    # Implementation: Find any elements in matrix where matrix[j, i] >= 8/(j+1).
    locations = np.zeros(len(seq), dtype=np.int64)
    for j in range(6):
        u_len = j + 1
        repeated_locations = np.logical_not(locations) & (matrix[j,:] >= NAMING_OPTIONS["min_repeat_length"] / u_len)
        locations[repeated_locations] = u_len

    # locations[i] now contains the u_len of a repeat ending at position i (inclusive!), if that repeat is 8+ nt, else 0.
    # Note that repeats of 8+ nt using different-length units cannot end in the same place.

    start, end, units = grow_scope(seq, matrix, locations, start, end)

    # We now have a well-defined scope!
    range_seq = seq[start : end]
    offset += start

    # Find all repeats of the significant units within the scope range.
    # NOTE: We will not get longer-unit singletons overlapping a 4+ repeat, e.g., ACTA[1]'s within ACT[4+].
    # NOTE: In strnaming-1.1.4 we also wouldn't get repeats that completely overlap a longer unit, e.g., A[3] within AAAG[x] or AGAA[2+].
    # That latter filter has been changed to apply only if there is overlap with a 4+ repeat, e.g., AAAG[4+] or AGAA[4+].
    repeats = find_repeat_stretches(range_seq, sorted(units, key=len, reverse=True), True, units)

    # TODO: Now maybe we can do a possible performance trick:
    # Drop all repeats that can't be reached from a significant repeat (8+ nt and 4+ repeats) with at most 5 gaps of at most 8 nt.
    # Subsequently also drop all repeats that can't be anchored anymore. Loop to drop more.

    # Perform repeat stretch trimming such that partially overlapping repeats can be combined (maintain minimum 3+ nt).
    trim_overlapping_repeats(repeats, units)

    # Add (repeat_length, unit_length, anchor, preferred) to the end of all repeats.
    # We will set preferred=True on everything and anchor=True if repeat length >= 8.
    # This means that the structure may start and end with any stretch of any length,
    # but requires that there is at least one repeat of length >= 8 for each unit used.
    # We will consider all possible locations of all units that have at least one
    # stretch of 8+ nt.
    for repeat in repeats:
        repeat_start, repeat_end, repeat_unit = repeat
        repeat_length = repeat_end - repeat_start
        repeat.append((repeat_length, len(repeat_unit), repeat_length >= 8, True))

    # Find the best combination of repeat stretches.
    # Doing this in a loop to make sure that the ref allele gets the same
    # repeat structure as the foregoing reference sequence analysis suggested.
    # In some cases, the ref allele name has a higher score than the initial outcome,
    # as allele naming has more freedom (e.g., large gap, repeat anchoring rules).
    previous_path_jsons = []  # FIXME, temp
    previous_paths = []  # FIXME, temp
    previous_path = None
    while True:
        # Short-circuit if we can't collapse anything.
        if not repeats:
            return 0, []

        if not previous_path:
            # Bootstrap values.
            prefix = None
            suffix = None
            block_length = 0
            large_gap_length = 0
        else:
            prefix = range_seq[:previous_path[0][0]-offset]
            suffix = range_seq[previous_path[-1][1]-offset:]
            block_length = find_block_length(previous_path, len_or_int=len)
            large_gap = find_overlong_gap(previous_path)
            large_gap_length = large_gap[1] - large_gap[0] if large_gap else 0

        # Find the best combination of repeat unit usage.
        score, path = get_best_path(prefix, suffix, range_seq, repeats, units, block_length, large_gap_length, offset=offset)

        if path and previous_path != path:
            this_path = json.dumps(path)
            if this_path in previous_path_jsons:
                # FIXME, temp: error out if we get stuck in a circle.
                osci_str = "\n".join((
                    "prev_path=%r" % (json.dumps(previous_path),),
                    "this_path=%r" % (this_path,),
                    "previous_paths=%r" % (previous_path_jsons,),
                    "involves_block_length=%r" % (len({find_block_length(p, len_or_int=len) for p in previous_paths}) > 1,),
                    "involves_short_only=%r" % (any((set(u for s,e,u in p) - set(u for s,e,u in p if e-s>=8)) for p in previous_paths),),
                    "involves_large_gap=%r" % (any(find_overlong_gap(p) for p in previous_paths),)))
                raise OscillationException(osci_str)
            previous_path_jsons.append(this_path)
            previous_paths.append(path)
            previous_path = path

            # Repeat find_everything() and get_best_path(), now treating the used units as preferred.
            # Singletons in the path are not counted as used units.
            # This is doing 'real' allele naming, not applying specific refseq rules.
            units = set(unit for s, e, unit in path if e - s > len(unit))
            repeats = find_everything(range_seq, {unit: [(0, len(range_seq))] for unit in units})
        else:
            return score, path
#collapse_repeat_units_refseq


def wrap_collapse_repeat_units_refseq(seq, *, offset=0):
    try:
        score, path = collapse_repeat_units_refseq(seq, offset=offset)
    except OscillationException as ex:
        return None, ex
    except OutOfTimeException as ex:
        return None, ex
    except ComplexityException as ex:
        return None, ex
    if not path or path[-1][1] - path[0][0] < NAMING_OPTIONS["min_structure_length"]:
        return None, None
    else:
        return path, None
#wrap_collapse_repeat_units_refseq


def recurse_collapse_repeat_units_refseq(seq, *, offset=0, status_callback=None):
    ranges = [(seq, offset)]
    while ranges:
        sub_seq, start = ranges.pop()
        structure, exception = wrap_collapse_repeat_units_refseq(sub_seq, offset=start)
        end = start + len(sub_seq)
        if structure:
            yield structure

            # Maybe we can make structures before this one.
            range_before = structure[0][0] - start
            if range_before >= NAMING_OPTIONS["min_structure_length"]:
                sub_seq = seq[start-offset : structure[0][0]-offset]
                ranges.append((sub_seq, start))

            # Maybe we can make structures after this one.
            range_after = end - structure[-1][1]
            if range_after >= NAMING_OPTIONS["min_structure_length"]:
                sub_seq = seq[structure[-1][1]-offset : end-offset]
                ranges.append((sub_seq, structure[-1][1]))

        if status_callback is not None:
            status_callback(start, end, structure, exception)
#recurse_collapse_repeat_units_refseq
