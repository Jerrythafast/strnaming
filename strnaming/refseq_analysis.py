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
import re
import sys

from .libstrnaming import REFSEQ_MINIMUM_REPEATS, NAMING_OPTIONS, \
                          ComplexityException, OutOfTimeException, \
                          get_best_path, find_block_length, find_everything, find_overlong_gap, \
                          find_repeat_stretches, trim_overlapping_repeats
from .refseq_cache import get_refseq

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
    sa = np.fromiter(map(ord, seq), count=len(seq), dtype=np.int32)
    matrix = np.empty((NAMING_OPTIONS["max_unit_length"], len(seq)), dtype=np.int32)
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
    return int(start), int(end)  # Convert from np.int32 to native Python int
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
    """
    Return the optimal score and repeat structure containing the longest repeat in seq.
    """
    if len(seq) < NAMING_OPTIONS["min_structure_length"]:
        return 0, []

    # Import numpy now.
    global np
    import numpy as np

    # First, find the longest repeat.
    matrix = detect_repeats_np(seq)
    start, end, u_len = find_longest_repeat_np(matrix)

    # If the longest repeat is not significant (8+ nt and 4+ repeats), stop.
    if end - start < max(NAMING_OPTIONS["min_repeat_length"], u_len * REFSEQ_MINIMUM_REPEATS):
        return 0, []

    # Find all stretches of at least 8 nt.
    # Implementation: Find any elements in matrix where matrix[j, i] >= 8/(j+1).
    locations = np.zeros(len(seq), dtype=np.int32)
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
    repeat_list = find_repeat_stretches(range_seq, sorted(units, key=len, reverse=True), True, units)

    # TODO: Now maybe we can do a possible performance trick:
    # Drop all repeats that can't be reached from a significant repeat (8+ nt and 4+ repeats) with at most 5 gaps of at most 8 nt.
    # Subsequently also drop all repeats that can't be anchored anymore. Loop to drop more.

    # Perform repeat stretch trimming such that partially overlapping repeats can be combined (maintain minimum 3+ nt).
    trim_overlapping_repeats(repeat_list, units)

    # Add (repeat_length, unit_length, anchor, preferred) to the end of all repeats.
    # We will set preferred=True on everything and anchor=True if repeat length >= 8.
    # This means that the structure may start and end with any stretch of any length,
    # but requires that there is at least one repeat of length >= 8 for each unit used.
    # We will consider all possible locations of all units that have at least one
    # stretch of 8+ nt.
    for repeat in repeat_list:
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
        if not repeat_list:
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
        score, path = get_best_path(prefix, suffix, range_seq, repeat_list, units, block_length, large_gap_length, offset=offset)

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
            repeat_list = find_everything(range_seq, {unit: [(0, len(range_seq))] for unit in units})
        else:
            return score, path
#collapse_repeat_units_refseq


def wrap_collapse_repeat_units_refseq(seq, *, offset=0):
    """
    Like collapse_repeat_units_refseq, but return a 2-tuple (structure, exception)
    instead. If OscillationException, OutOfTimeException or ComplexityException was
    raised, return (None, exception). Else, return (structure, None).
    If no repeat structure was found, (None, None) is returned.
    """
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


def recurse_collapse_repeat_units_refseq(seq, *, offset=0, status_callback=None, generate=True):
    """
    Repeatedly call wrap_collapse_repeat_units_refseq to identify all repeat structures in
    the given sequence.

    If provided, the status_callback function is called with four arguments
    (start, exclusive_end, structure, exception) to communicate results on sub-ranges.
    Subsequent calls do not necessarily describe the results in a particular order.

    If generate=True (the default), this function will yield (as a generator) the
    identified repeat structures.
    """
    ranges = [(seq, offset)]
    while ranges:
        sub_seq, start = ranges.pop()
        structure, exception = wrap_collapse_repeat_units_refseq(sub_seq, offset=start)
        end = start + len(sub_seq)
        if structure:
            if generate:
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


def get_scopes_from_seq(seq):
    """Return a sorted list of unique, possibly-overlapping scopes on seq."""
    # Build an overview of all 4+ repeat locations of 8+ nt.
    import numpy as np
    matrix = detect_repeats_np(seq)
    locations_4r_8nt = np.zeros(len(seq), dtype=np.int32)
    locations_8nt = np.zeros(len(seq), dtype=np.int32)
    for j in range(6):
        u_len = j + 1
        row = matrix[j,:]
        repeated_locations = np.logical_not(locations_4r_8nt) & (row >= max(8 / u_len, 4))
        locations_4r_8nt[repeated_locations] = u_len
        repeated_locations = np.logical_not(locations_8nt) & (row >= 8 / u_len)
        locations_8nt[repeated_locations] = u_len

    # locations_4r_8nt[i] now contains the u_len of a repeat ending at position i (inclusive!), if that repeat is 8+ nt and at least 4 repeats, else 0.
    # locations_8nt[i] now contains the u_len of a repeat ending at position i (inclusive!), if that repeat is 8+ nt, else 0.

    # Iterate over *ALL* 4+ repeats of 8+ nt to collect all unique scopes.
    scopes = set()
    for end in np.flatnonzero(locations_4r_8nt):
        ulen = locations_4r_8nt[end]
        start = end + 1 - matrix[ulen - 1, end] * ulen
        start, end, units = grow_scope(seq, matrix, locations_8nt, start, end + 1)
        if end - start >= NAMING_OPTIONS["min_structure_length"]:
            scopes.add((int(start), int(end)))
    return sorted(scopes)
#get_scopes_from_seq


def find_scopes(seq, offset):
    """Yield unique, non-overlapping scopes on seq, counting from the given offset."""
    scopes = get_scopes_from_seq(seq)
    if not scopes:
        return

    # Merge overlapping scopes while generating.
    start, end = scopes[0]
    for next_start, next_end in itertools.islice(scopes, 1, None):
        if next_start >= end:
            yield start + offset, end + offset
            start = next_start
            end = next_end
        elif next_end > end:
            end = next_end
    yield start + offset, end + offset
#find_scopes


def analyse_range(threads, chr, start, end, status_callback):
    """
    Analyse the given range of genomic sequence.
    The end position is exlusive.

    The provided status_callback function is called with four
    arguments (start, exclusive_end, structure, exception) to
    communicate results on sub-ranges. This function is called,
    possibly concurrently, from within subprocesses.
    """
    from concurrent.futures import ProcessPoolExecutor
    with ProcessPoolExecutor(threads) as pool:
        seq = get_refseq(chr, start, end - 1)
        pos = start
        for scope_start, scope_end in find_scopes(seq, start):
            if scope_start > pos:
                status_callback(pos, scope_start, None, None)
            pool.submit(recurse_collapse_repeat_units_refseq,
                seq[scope_start-start:scope_end-start],
                offset=scope_start, status_callback=status_callback, generate=False)
            if scope_end < end:
                status_callback(scope_end, end, None, None)
            pos = scope_end
#analyse_range


def print_callback(chr, start, end, structure, exception):
    """Print structures to stdout and other statuses to stderr."""
    if structure:
        print("%s\t%s" % (chr, "\t".join(str(x) for x in structure for x in x)))
        sys.stdout.flush()
    else:
        status = "DEFINITE"
        exception_type = type(exception)
        if exception_type is OscillationException:
            status = "FAILOSCI"
        elif exception_type is OutOfTimeException:
            status = "FAILTIME"
        elif exception_type is ComplexityException:
            status = "FAILSIZE"
        sys.stderr.write("%s chr%s:%i..%i\n" % (status, chr, start, end - 1))
        sys.stderr.flush()
#print_callback


if __name__ == "__main__":
    import argparse
    import functools
    import libstrnaming
    parser = argparse.ArgumentParser()
    parser.add_argument("--chr", required=True)
    parser.add_argument("--start", type=int, required=True)
    parser.add_argument("--end", type=int, required=True)
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--maxtime", type=int, default=300)
    args = parser.parse_args()
    libstrnaming.MAX_SECONDS = args.maxtime
    libstrnaming.MAX_SECONDS_REFSEQ = args.maxtime
    chr = args.chr.upper().lstrip("CHR")
    callback = functools.partial(print_callback, chr)
    analyse_range(args.threads, chr, args.start, args.end, callback)
