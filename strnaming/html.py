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

from . import classes


HTML_START = """<!DOCTYPE html>
<html>
<head>
<title>STRNaming Results</title>
<style type="text/css">
table {
    font-family: monospace;
    text-align: left;
    white-space: nowrap;
    margin-bottom: 2em;
}
h1 {
    font-family: sans-serif;
    margin: 0;
}
tr {
    background-color: rgba(0, 91, 142, 0.2);
}
tr:hover {
    font-weight: bold;
    background-color: rgba(0, 91, 142, 0.5);
    transition: background-color 0.2s ease;
}
.seqrepeat0 {
    background-color: rgba(176, 80, 48, 0.75);
}
.seqrepeat1 {
    background-color: rgba(176, 144, 80, 0.75);
}
.seqrepeat2 {
    background-color: rgba(160, 176, 48, 0.75);
}
.seqrepeat3 {
    background-color: rgba(48, 144, 80, 0.75);
}
.seqrepeat4 {
    background-color: rgba(144, 80, 144, 0.75);
}
.alignright td, .alignright th {
    text-align: right;
}
.notice {
    position: fixed;
    left: 0;
    bottom: 0;
    background-color: white;
    font-family: sans-serif;
}
</style>
</head>
<body ondblclick="this.classList.toggle('alignright')">
<div class="notice">Double-click to toggle text alignment.</div>\n"""
MARKER_START = "<h1>%s</h1>"
TABLE_START = "<table><tr><th>Marker</th><th>Allele Name</th><th>Sequence</th></tr>\n"
SPAN_TEMPLATE = "<span class='seqrepeat%i'>%%s</span>"
TABLE_ROW = "<tr><td>%s</td><td>%s</td><td>%s</td></tr>\n"
TABLE_END = "</table>\n"
MARKER_END = ""
HTML_END = "</body>\n</html>\n"
NUM_CLASSES = 5


def colorize_refseq(structure_store, marker_range, *, context=50):
    # Get extended range with additional context sequence.
    chromosome, start, end = marker_range.location
    if context:
        start -= context
        end += context + 1
        structure_store.load_within_range(chromosome, start, end)
        wider_range = classes.ReportedRange(structure_store, chromosome, start, end)
    else:
        end += 1
        wider_range = marker_range

    # Get reference sequence and a copy with context in lowercase.
    refseq_upper = refseq = structure_store.get_refseq_store().get_refseq(chromosome, start, end)
    if context:
        refseq = refseq[:context].lower() + refseq[context:-context] + refseq[-context:].lower()

    # Colorize reference sequence.
    normalized_seq = wider_range.normalize_sequence(refseq, lowercase_inserts=True)
    normalized_upper = normalized_seq.upper()

    stretches = wider_range.get_stretches(refseq_upper, normalized_seq=normalized_upper)
    colored_name = wider_range.get_name(refseq_upper, color=True, normalized_seq=normalized_upper,
        stretches=stretches)
    colored_sequence = wider_range.colorize_sequence(refseq_upper, normalized_seq=normalized_seq,
        stretches=stretches)
    return colored_name, colored_sequence


def colorize(marker_range, sequence):
    stretches = marker_range.get_stretches(sequence)
    colored_name = marker_range.get_name(sequence, color=True, stretches=stretches)
    colored_sequence = marker_range.colorize_sequence(sequence, stretches=stretches)
    return colored_name, colored_sequence


def make_page(ranges_store, instream, outstream, *,
        page_headers=True, marker_headers=True, marker_refseqs=True):
    if page_headers:
        outstream.write(HTML_START)
    if not marker_headers:
        outstream.write(TABLE_START)
    prev_marker = None
    for marker, sequence in (line.rstrip("\r\n").split() for line in instream):
        marker_range = ranges_store.get_range(marker)
        if marker_headers and not prev_marker == marker:
            if prev_marker is not None:
                outstream.write(TABLE_END)
                outstream.write(MARKER_END)
            outstream.write(MARKER_START % marker)
            outstream.write(TABLE_START)
        if marker_refseqs and not prev_marker == marker:
            colored_name, colored_sequence = colorize_refseq(
                ranges_store.get_structure_store(), marker_range)
            outstream.write(TABLE_ROW %
                ("REF" if marker_headers else marker, colored_name, colored_sequence))
        colored_name, colored_sequence = colorize(marker_range, sequence)
        outstream.write(TABLE_ROW % (marker, colored_name, colored_sequence))
        prev_marker = marker
    if not marker_headers or prev_marker is not None:
        outstream.write(TABLE_END)
    if marker_headers:
        outstream.write(MARKER_END)
    if page_headers:
        outstream.write(HTML_END)
