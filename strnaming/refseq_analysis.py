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

import sys

from .libstrnaming import NAMING_OPTIONS, find_repeated_units, find_block_length, \
                          find_overlong_gap, get_best_path, find_everything

import json #FIXME, temp


class OscillationException(Exception):
    pass


def ranges_of(repeats, len_seq):
    if not repeats:
        return []
    ranges = []
    range_start = None
    range_end = None
    for start, end, unit, infos in repeats:
        if range_start is None:
            # Start new range.
            range_start = max(0, start - NAMING_OPTIONS["max_long_gap"])
            range_end = min(len_seq, end + NAMING_OPTIONS["max_long_gap"])
        elif start <= range_end:
            # Extend current range.
            range_end = min(len_seq, max(range_end, end + NAMING_OPTIONS["max_long_gap"]))
        else:
            # Illegally large gap between repeats.
            if range_end - range_start >= NAMING_OPTIONS["min_structure_length"]:
                ranges.append((range_start, range_end))
            range_start = max(0, start - NAMING_OPTIONS["max_long_gap"])
            range_end = min(len_seq, end + NAMING_OPTIONS["max_long_gap"])
    ranges.append((range_start, range_end))
    return ranges
#ranges_of


def collapse_repeat_units_refseq(seq, *, allow_close=0):
    """Combine repeat units in the given sequence."""

    # Find some units to start with. These will be treated as 'preferred' in the first round.
    unit_locations = find_repeated_units(seq, True)

    previous_paths = set()  # FIXME, temp
    previous_path = None
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

        # Initially, repeats can have multiple disjoint ranges.
        ranges = [] if previous_path else ranges_of(repeats, len(seq))
        if len(ranges) > 1:
            # There is an illegally large gap between the repeats;
            # re-analyse this sequence in separate parts instead.
            score = -sys.maxsize
            path = []
            for range_start, range_end in ranges:
                range_seq = seq[range_start : range_end]
                unit_locations = find_repeated_units(range_seq, True)
                preferred_units = sorted(unit_locations, key=len, reverse=True)
                repeats = find_everything(range_seq, unit_locations)
                new_score, new_path = get_best_path(None, None, range_seq, repeats, preferred_units, 0, 0)
                if new_score > score:
                    score = new_score
                    path = [[start + range_start, end + range_start, unit]
                        for start, end, unit in new_path] if range_start else new_path
        else:
            # Find the best combination of repeat unit usage.
            if not previous_path:
                # Bootstrap values.
                prefix = None
                suffix = None
                block_length = 0
                large_gap_length = 0
            else:
                prefix = seq[:previous_path[0][0]]
                suffix = seq[previous_path[-1][1]:]
                block_length = find_block_length(previous_path, len_or_int=len)
                large_gap = find_overlong_gap(previous_path)
                large_gap_length = large_gap[1] - large_gap[0] if large_gap else 0
            score, path = get_best_path(prefix, suffix, seq, repeats, preferred_units, block_length, large_gap_length)

        # Log a warning if repeats occur very close to either end of the sequence.
        # The name may not be stable if this is the case.
        if path and path[0][0] < 30 and not allow_close & 1:
            sys.stderr.write("Path #%i starts at position %i in analysed sequence\n" % (len(previous_paths) + 1, path[0][0]))
        if path and len(seq) - path[-1][1] < 30 and not allow_close & 2:
            sys.stderr.write("Path #%i ends at %i positions before end of sequence\n" % (len(previous_paths) + 1, len(seq) - path[-1][1]))

        if path and previous_path != path:
            # FIXME, temp:
            this_path = json.dumps(path)
            if this_path in previous_paths:
                sys.stderr.write("prev_path=%r\n" % (json.dumps(previous_path),))
                sys.stderr.write("this_path=%r\n" % (this_path,))
                sys.stderr.write("previous_paths=%r\n" % (previous_paths,))
                raise OscillationException()
            previous_paths.add(this_path)

            # Repeat find_everything() and get_best_path(), now treating the used units as preferred.
            # Singletons in the path are not counted as used units.
            used_units = set(unit for start, end, unit in path if end - start > len(unit))
            unit_locations = {unit: [(0, len(seq))] for unit in used_units}
            previous_path = path
        else:
            return score, path
#collapse_repeat_units_refseq


def recurse_collapse_repeat_units_refseq(seq, *, offset=0, allow_close=0):
    score, path = collapse_repeat_units_refseq(seq, allow_close=allow_close)
    if path and path[-1][1] - path[0][0] >= NAMING_OPTIONS["min_structure_length"]:

        # Maybe we can make structures before this one.
        pos = path[0][0]
        if pos >= NAMING_OPTIONS["min_structure_length"]:
            yield from recurse_collapse_repeat_units_refseq(seq[:pos], offset=offset, allow_close=allow_close|2)

        yield [[repeat[0] + offset, repeat[1] + offset, repeat[2]] for repeat in path]

        # Maybe we can make structures after this one.
        pos = path[-1][1]
        if len(seq) - pos >= NAMING_OPTIONS["min_structure_length"]:
            yield from recurse_collapse_repeat_units_refseq(seq[pos:], offset=pos+offset, allow_close=allow_close|1)
#recurse_collapse_repeat_units_refseq
