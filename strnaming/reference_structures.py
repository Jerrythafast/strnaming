#
# Copyright (C) 2021 Jerry Hoogenboom
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


def _txt2bin(instream):
    """Included for reference."""
    outfiles = {}
    last_end = {}
    for line in instream:
        chr, *structure = line.strip("\r\n").split("\t")
        if chr not in outfiles:
            filename = REFFILENAME.format(chromosome=chr)
            reffile = BUILTINDIR / filename
            outfiles[chr] = gzip.open(str(reffile), "xb")
            last_end[chr] = 0
        structure = list(structure)
        for i in range(0, len(structure), 3):
            structure[i] = int(structure[i])
            structure[i+1] = int(structure[i+1])
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


def _parse_header(header):
    start = ((header[3] & 0x0F) << 24) | (header[2] << 16) | (header[1] << 8) | (header[0])
    length = ((header[3] & 0xF0) << 4) | header[4]
    num_bytes = header[5]
    return (start, start + length, num_bytes)


def get_within_range(chromosome, start, end):
    """
    Get all structures in the given range of the genome. The end position is inclusive.
    """
    end += 1
    structures = []
    filename = REFFILENAME.format(chromosome=chromosome)
    reffile = BUILTINDIR / filename
    if not reffile.is_file():
        return structures
    with gzip.open(str(reffile), "rb") as instream:
        while True:
            # Read header.
            header = instream.read(6)
            if not header:
                break  # EOF
            s_start, s_end, num_bytes = _parse_header(header)
            if s_start >= end:
                break  # Done (assuming input is sorted by s_start).
            if s_end <= start:
                instream.seek(num_bytes, whence=1)
            else:
                # Structure within range.
                this_bytes = instream.read(num_bytes)
                structure = []
                pos = s_start
                i = 0
                while i < num_bytes:
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
                if pos != s_end:
                    raise ValueError("Possible corruption in file %s" % reffile)
                structures.append(structure)
    return structures
