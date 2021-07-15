#!/usr/bin/env python3

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

import argparse
import io
import pkgutil
import re
import sys
import textwrap
#import cProfile  # Imported only if the -d/--debug option is specified


from . import usage, version, classes, refseq_cache, html


OPTIMIZE = False  # FIXME; temp. Turn this to 'True' to run optimize.py.
if OPTIMIZE:
    from . import optimize




class _HelpFormatter(argparse.HelpFormatter):
    _pat_paragraph_delim = re.compile("\n\n+")
    def _fill_text(self, text, width, indent):
        # Reflow (wrap) description text, but maintain paragraphs.
        return "\n\n".join(
            textwrap.fill(self._whitespace_matcher.sub(" ", p).strip(), width,
                          initial_indent=indent, subsequent_indent=indent)
            for p in self._pat_paragraph_delim.split(text))
    #_fill_text
#_HelpFormatter


def stream_package_data(filename):
    return io.StringIO(pkgutil.get_data("strnaming", "data/" + filename).decode())
#stream_package_data


def run_strnaming(arguments):
    """
    Run STRNaming command with the given arguments.
    """
    ranges_store = classes.ReportedRangeStore(
        stream_package_data("ranges_%s.txt" % arguments.ranges),
        load_structures=True)
    if arguments.html:
        html.make_page(ranges_store, arguments.instream, arguments.outstream,
            page_headers=arguments.html_page_headers,
            marker_headers=arguments.html_marker_headers,
            marker_refseqs=arguments.html_marker_refseqs)
        return
    for marker, sequence in (line.rstrip("\r\n").split() for line in arguments.instream):
        arguments.outstream.write("%s\t%s\n" % (marker, ranges_store.get_range(marker).get_name(sequence)))
        if arguments.unbuffered:
            arguments.outstream.flush()
#run_strnaming


def run_refseq_cache(arguments):
    """
    Cache a portion of reference sequence.
    """
    rng = re.match(r"(?:[Cc][Hh][Rr])?([1-9XYM]|1[0-9]|2[0-2]):(\d+)\.\.(\d+)", arguments.range)
    if rng is None:
        raise ValueError("Invalid range specification '%s'" % arguments.range)
    try:
        refseq_cache.get_refseq(rng.group(1), int(rng.group(2)), int(rng.group(3)))
    except Exception as e:
        print("Something went wrong while trying to cache refseq range %s" % arguments.range)
        raise
    print("The reference sequence was successfully cached. The cache location is: %s"
        % refseq_cache.CACHEDIR)
#run_refseq_cache


def main():
    """
    Main entry point.
    """
    parser = argparse.ArgumentParser(formatter_class=_HelpFormatter, prog="strnaming",
        description=usage[0])
    parser.add_argument("-v", "--version", action="version", version=version("strnaming"),
        help="show version number and exit")
    parser.add_argument("-d", "--debug", action="store_true",
        help="if specified, additional debug output is given")
    if OPTIMIZE:
        optimize.add_arguments(parser)
        parser.set_defaults(func=optimize.run)
    else:
        subparsers = parser.add_subparsers(title="available actions", metavar="ACTION")
        subparser = subparsers.add_parser("name-sequences", formatter_class=_HelpFormatter,
            help="convert sequences to allele names (default)",
            description="generate allele names for one or more DNA sequences")
        subparser.add_argument("-U", "--unbuffered", action="store_true",
            help="if specified, the output buffer is flushed after every line")
        subparser.add_argument("-r", "--ranges", metavar="uas-frr", choices=("uas-frr",),
            required=True,
            help="specify 'uas-frr' to use UAS flanking region report ranges (mandatory; future "
                 "versions of STRNaming will offer more options)")
        subparser.add_argument("instream", metavar="IN", nargs="?", type=argparse.FileType("tr"),
            default="-",
            help="an input file or stream containing a marker name and a raw sequence on each line, "
                 "separated by whitespace, which will be used to process in STRNaming (default: read "
                 "from stdin)")
        subparser.add_argument("outstream", metavar="OUT", nargs="?", type=argparse.FileType("tw"),
            default=sys.stdout,
            help="output file, always uses a tab character as separator (default: write to stdout)")
        html_group = subparser.add_argument_group("colorized output (HTML format)")
        html_group.add_argument("--html", action="store_true",
            help="write an HTML-formatted output file with colorized repeats")
        html_group.add_argument("--html-no-head", dest="html_page_headers", action="store_false",
            help="exclude <head> section from HTML output; only write visible page content")
        html_group.add_argument("--html-no-headings", dest="html_marker_headers",
            action="store_false", help="exclude per-marker headings in the HTML output")
        html_group.add_argument("--html-no-refseq", dest="html_marker_refseqs",
            action="store_false", help="exclude marker reference sequence in the HTML output")
        subparser.set_defaults(func=run_strnaming)
        subparser = subparsers.add_parser("refseq-cache", formatter_class=_HelpFormatter,
            help="download and cache reference sequence",
            description="Makes sure the given portion of reference sequence is available in "
                        "STRNaming's offline reference sequence cache, downloading the sequence "
                        "from the Ensembl REST API (http://rest.ensembl.org) if necessary.")
        subparser.add_argument("range", metavar="RANGE",
            help="range of sequence to download, for example 'chr2:1489653..1489689'")
        subparser.set_defaults(func=run_refseq_cache)
    try:
        # Assume the user wants help if they just type 'strnaming'.
        if len(sys.argv) == 1:
            sys.argv.append("-h")
        args, unknowns = parser.parse_known_args()
    except Exception as error:
        parser.error(error)
    try:
        if unknowns:
            # Politely inform the user about unknown arguments.
            parser.error(
                "The following arguments are not known. Please check spelling "
                "and argument order: '%s'." % "', '".join(unknowns))
        if args.debug and not OPTIMIZE:
            import cProfile
            cProfile.runctx("args.func(args)", globals(), locals(), sort="tottime")
        else:
            args.func(args)
    except Exception as error:
        if args.debug:
            raise
        parser.error(error)
#main


if __name__ == "__main__":
    main()
