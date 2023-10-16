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
import gzip

from pathlib import Path

# Constants for file location.
BUILTINDIR = Path(__file__).parent / "data"
REFFILENAME = "reference-structures-{chromosome}.bin.gz"


# HEADER 48 BITS = 6 BYTES:
# Start position requires 28 bits
# Length fits within 12 bits (0..4095, max=651)
# Size fits within 8 bits (0..255, max=57)

# BLOCKS ARE 1 OR 2 BYTES:     (unit length 0 indicates interruption)
# Short: uuur rrrr             (u=unit length 0..6; r=repeat length 1..32)
# Long:  111u uurr  rrrr rrrr  (u=unit length 0..6; r=repeat length 1..1024)


def _parse_structures_with_separation(instream, min_separation):
    """
    Parse and yield structures from instream, requiring at least
    min_separation nucleotides between any two structures.
    Print (to stdout) the lines containing structures that are
    too close to each other.
    The input must be sorted.
    """
    prev = None
    for line in instream:
        line = line.strip("\r\n")
        chr, *structure = line.split("\t")
        structure = list(structure)
        for i in range(0, len(structure), 3):
            structure[i] = int(structure[i])
            structure[i+1] = int(structure[i+1])
        too_close = False
        if prev:
            prev_chr, prev_structure, prev_line, prev_too_close = prev
            sep = min_separation if prev_chr != chr else structure[0] - prev_structure[-2]
            if sep < 0:
                raise ValueError(
                    "Bad structure order: Chr%s, %i-%i, %i-%i" %
                    (chr, prev_structure[0], prev_structure[-2], structure[0], structure[-2]))
            if sep < min_separation:
                # Too close to previous structure!
                if not prev_too_close:
                    print(prev_line)
                print(line)
                too_close = True
            elif not prev_too_close:
                # We can safely output the previous structure now.
                yield prev_chr, prev_structure
        prev = (chr, structure, line, too_close)
    if prev and not prev[-1]:
        # Last structure was not too close, output it!
        yield prev[0], prev[1]

def _filter_large_gaps(structures):
    """
    Yield structures that don't contain a large interruption.
    Print (to stdout) the structures that have one.
    """
    from .libstrnaming import find_overlong_gap
    for chr, structure in structures:
        if find_overlong_gap(list(zip(*[iter(structure)]*3))) is not None:
            print("%s\t%s" % (chr, "\t".join(map(str, structure))))
        else:
            yield chr, structure

def _txt2bin(structures):
    """
    Read structures from the given iterable and write them to the
    built-in reference structure data files.
    The input must be sorted.
    """
    outfiles = {}
    last_end = {}
    for chr, structure in structures:
        if chr not in outfiles:
            filename = REFFILENAME.format(chromosome=chr)
            reffile = BUILTINDIR / filename
            outfiles[chr] = gzip.open(str(reffile), "xb")
            last_end[chr] = 0
        start, *_, end, _ = structure
        if start < last_end[chr]:
            raise ValueError(
                "Bad structure order: chr%s, ?-%i, %i-%i" % (chr, last_end[chr], start, end))
        last_end[chr] = end
        length = end - start
        if start > 0x0FFFFFFF:
            raise ValueError("Unexpectedly high start position")
        if length > 0x0FFF:
            raise ValueError("Unexpectedly long structure")
        this_bytes = bytearray()
        pos = structure[0]
        for u_start, u_end, unit in zip(*[iter(structure)]*3):
            if u_start > pos:
                gap_len_minus_one = u_start - pos - 1
                if gap_len_minus_one > 0x1F:
                    # Could extend to two bytes, but:
                    raise ValueError("Unexpectedly long gap")
                this_bytes.append(gap_len_minus_one)
            unit_len = len(unit)
            repeats_minus_one = (u_end - u_start) // unit_len - 1
            if unit_len > 6:
                raise ValueError("Unexpectedly long repeat unit")
            if repeats_minus_one < 0x1F:
                # Pack unit_len in upper 3 bits, repeats in lower 5 bits.
                this_bytes.append((unit_len << 5) | repeats_minus_one)
            elif repeats_minus_one < 0x03FF:
                # First byte: upper 3 bits are 111, then 3 bits unit_len,
                # then upper 2 bits of repeats.
                # Second byte: lower 8 bits of repeats.
                this_bytes.append(0xE0 | (unit_len << 2) | (repeats_minus_one >> 8))
                this_bytes.append(repeats_minus_one & 0xFF)
            else:
                raise ValueError("Unexpectedly long repeat stretch")
            pos = u_end
        if len(this_bytes) > 0xFF:
            raise ValueError("Unexpectedly long structure description")
        header = bytearray()
        header.append(start & 0xFF)
        header.append((start >> 8) & 0xFF)
        header.append((start >> 16) & 0xFF)
        header.append(((start >> 24) & 0x0F) | ((length >> 4) & 0xF0))
        header.append(length & 0xFF)
        header.append(len(this_bytes))
        outfiles[chr].write(header)
        outfiles[chr].write(this_bytes)
    for outfile in outfiles.values():
        outfile.close()


class Reader:
    def __init__(self, chromosome):
        self.reffile = BUILTINDIR / REFFILENAME.format(chromosome=chromosome)

    def __enter__(self):
        self.s_start = self.s_end = self.num_bytes = None
        if self.reffile.is_file():
            self.instream = gzip.open(str(self.reffile), "rb")
            self.next_header()
        else:
            self.instream = None
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        if self.instream is not None:
            self.instream.close()
            self.s_start = self.s_end = self.num_bytes = None

    def next_header(self):
        if self.num_bytes:
            # Seek to next header.
            self.instream.seek(self.num_bytes, whence=1)
        header = self.instream.read(6)
        if header:
            self.s_start = (((header[3] & 0x0F) << 24) |
                (header[2] << 16) | (header[1] << 8) | header[0])
            self.s_end = self.s_start + (((header[3] & 0xF0) << 4) | header[4])
            self.num_bytes = header[5]
            return True
        else:  # EOF
            self.s_start = self.s_end = self.num_bytes = None
            return False

    def next_structure(self):
        this_bytes = self.instream.read(self.num_bytes)
        self.num_bytes = 0  # Have read them now.
        structure = []
        pos = self.s_start
        i = 0
        while i < len(this_bytes):
            byte = this_bytes[i]
            unit_len = byte >> 5
            if unit_len == 7:
                i += 1
                unit_len = (byte & 0x1C) >> 2
                repeats = ((byte & 0x03) << 8 | this_bytes[i]) + 1
            else:
                repeats = (byte & 0x1F) + 1
            if unit_len:
                stretch_end = pos + repeats * unit_len
                structure.append([pos, stretch_end, unit_len])
                pos = stretch_end
            else:
                pos += repeats  # Skip over interruption.
            i += 1
        if pos != self.s_end:
            raise ValueError("Possible corruption in file %s" % self.reffile)
        self.next_header()
        return structure

    def skip_to(self, start):
        if self.s_end is not None:
            while self.s_end <= start and self.next_header():
                pass

    def gen_until(self, end):
        while self.s_start is not None and self.s_start <= end:
            yield self.next_structure()


def gen_within_range(chromosome, start, end):
    """
    Generate all structures in the given range of the genome. The end position is inclusive.
    """
    with Reader(chromosome) as reader:
        reader.skip_to(start)
        yield from reader.gen_until(end)


def gen_within_ranges(chromosome, ranges):
    """
    Generate all structures in the given ranges of a chromosome.
    The ranges should be a sorted iterable of (start, end) pairs.
    The end positions are inclusive.
    """
    with Reader(chromosome) as reader:
        for start, end in ranges:
            reader.skip_to(start)
            yield from reader.gen_until(end)


if __name__ == "__main__":
    import sys
    if sys.argv[1] == "txt2bin":
        with open(sys.argv[2], "rt") as instream:
            _txt2bin(_filter_large_gaps(_parse_structures_with_separation(instream, 20)))
    elif sys.argv[1] == "bin2txt":
        with Reader(sys.argv[2]) as reader:
            for structure in reader.gen_until(sys.maxsize):
                print("%s\t%s" % (sys.argv[2], "\t".join(
                    "%i\t%i\t%s" % (stretch[0], stretch[1], "N" * stretch[2])
                    for stretch in structure)))
