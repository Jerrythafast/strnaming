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
import sys

# Patterns that match blocks of TSSV-style sequences and allele names.
PAT_TSSV_BLOCK = re.compile("([ACGT]+)\((\d+)\)")
PAT_ALLELENAME_BLOCK = re.compile("([ACGT]+)\[(\d+)\]")
PAT_ALLELENAME_GAP = re.compile("\[((?: ?\d+(?:\.1)?[ACGT-]+>[ACGT-]+)*)\]")

# Patterns that match a single variant.
PAT_VARIANT_STR = re.compile(
    "^(?P<pos>[-+]\d+)(?:\.(?P<ins>1))?"
    "(?P<old>(?:(?<=\.1)-)|(?<!\.1)[ACGT]+)>"
    "(?!(?P=old))(?P<new>[ACGT]+|-)$")
PAT_VARIANT_SNP = re.compile(
    "^(?P<pos>\d+)(?:\.(?P<ins>1))?"
    "(?P<old>(?:(?<=\.1)-)|(?<!\.1)[ACGTN]+)>"
    "(?!(?P=old))(?P<new>[ACGT]+|-)$")
PAT_VARIANT_MT = re.compile(
    "^(?P<old>(?P<a>[ACGT])|-?)"
    "(?P<pos>\d+)(?(a)|(?:\.(?P<ins>\d+))?)"
    "(?P<new>[ACGT-]|DEL)$")

# Complementary bases.
COMPL = str.maketrans("ACGTUMKRYBDHVacgtumkrybdhv", "TGCAAKMYRVHDBtgcaakmyrvhdb")


def reverse_complement(sequence):
    """Return the reverse complement of the given DNA sequence."""
    return sequence[::-1].translate(COMPL)


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

    result = (i // row_offset, max_pos // row_offset, max_score)

    # Store the result in the cache.
    align.cache[haystack, needle, prefix] = result
    return result
align.cache = {}
#align


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
def align_affine(template, sequence, *, match_score=1, mismatch_score=-3, gap_open_score=-7,
                 gap_extend_score=-2, cache=True, debug=False):
    """
    Perform a global alignment of sequence to template and return a
    list of variations detected.

    By default, the results of this function are cached.  Set cache to
    False to suppress caching the result and reduce memory usage.

    Setting debug to True will cause the alignment matrices to be
    printed to sys.stdout.  Be aware that they can be quite large.
    """
    try:
        return align_affine.cache[
            template, sequence, match_score, mismatch_score, gap_open_score, gap_extend_score]
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

        gap_open = matrix_match[i - 1] + gap_open_score
        gap_extend = matrix_gap1[i - 1] + gap_extend_score
        if gap_open > gap_extend:
            matrix_direction[i] |= A_HORZ_O  # Must exit M_GAP1 here.
            matrix_gap1[i] = gap_open
        else:
            matrix_gap1[i] = gap_extend

        gap_open = matrix_match[i - row_offset] + gap_open_score
        gap_extend = matrix_gap2[i - row_offset] + gap_extend_score
        if gap_open > gap_extend:
            matrix_direction[i] |= A_VERT_O  # Must exit M_GAP2 here.
            matrix_gap2[i] = gap_open
        else:
            matrix_gap2[i] = gap_extend

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
        align_affine.cache[template, sequence, match_score, mismatch_score, gap_open_score,
            gap_extend_score] = variants
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
            for j in range(max(len_template, len_sequence)):
                named_variants.append("%s%i%s%s" % (
                    template[pos_template + j] if j < len_template else "",
                    get_genome_pos(location, pos_template + min(j, len_template - 1)),
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


def call_variants(template, sequence, *, location=(None, 1), match_score=None, mismatch_score=None,
                  gap_open_score=None, gap_extend_score=None, cache=True, debug=False):
    if location[0] == "M":
        match_score = 1 if match_score is None else match_score
        mismatch_score = -1 if mismatch_score is None else mismatch_score
        gap_open_score = -2 if gap_open_score is None else gap_open_score
        gap_extend_score = -1 if gap_extend_score is None else gap_extend_score
    else:
        match_score = 1 if match_score is None else match_score
        mismatch_score = -3 if mismatch_score is None else mismatch_score
        gap_open_score = -7 if gap_open_score is None else gap_open_score
        gap_extend_score = -2 if gap_extend_score is None else gap_extend_score
    variants = align_affine(template, sequence,
        match_score=match_score, mismatch_score=mismatch_score, gap_open_score=gap_open_score,
        gap_extend_score=gap_extend_score, cache=cache, debug=debug)
    return name_variants(template, sequence, variants, location=location)
#call_variants


def mutate_sequence(seq, variants, *, location=None):
    """Apply the given variants to the given sequence."""
    if not isinstance(location, tuple) or len(location) < 2:
        pattern = PAT_VARIANT_STR
        location = (None, 0)
    elif location[0] == "M":
        pattern = PAT_VARIANT_MT
        location = (location[0], location[1] - 1) + tuple(location[2:])
    else:
        pattern = PAT_VARIANT_SNP
        location = (location[0], location[1] - 1) + tuple(location[2:])

    seq = [[]] + [[base] for base in seq]
    for variant in variants:
        vm = pattern.match(variant)
        if vm is None:
            raise ValueError("Unrecognised variant '%s'" % variant)
        pos = int(vm.group("pos"))
        ins = int(vm.group("ins") or 0)
        old = vm.group("old")
        new = vm.group("new")
        if old == "-":
            old = ""
        if new == "-" or new.lower() == "del":
            new = ""
        if pos < 0:
            pos += len(seq)
            if ins:
                pos -= 1
        pos = get_genome_pos(location, pos, invert=True)
        if pos < 0 or (pos == 0 and not ins) or pos >= len(seq):
            raise ValueError("Position of variant '%s' is outside sequence range" % variant)
        if not ins and old and old != "".join("".join(x[:1]) for x in seq[pos : pos + len(old)]):
            raise ValueError(
                "Incorrect original sequence in variant '%s'; should be '%s'!"
                % (variant, "".join("".join(x[:1]) for x in seq[pos : pos + len(old)])))
        elif not ins and not old:
            # MtDNA substitution with reference base omitted.
            old = "".join("".join(x[:1]) for x in seq[pos : pos + len(new)])
        if not ins:
            # Remove old bases, retaining those inserted between/after.
            seq[pos : pos + len(old)] = [[""] + x[1:] for x in seq[pos : pos + len(old)]]
            # Place new entirely in the position of the first old base.
            seq[pos][0] = new
        else:
            # Insert new exactly ins positions after pos.
            while len(seq[pos]) <= ins:
                seq[pos].append("")
            seq[pos][ins] = new
    return "".join("".join(x) for x in seq)
#mutate_sequence
