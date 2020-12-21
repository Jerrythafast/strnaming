#!/usr/bin/env python

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

import argparse
import io
import pkgutil
import re
import sys
import textwrap
#import cProfile  # Imported only if the -d/--debug option is specified


from . import usage, version, classes


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
    refseq_store = classes.ReferenceSequenceStore(stream_package_data("refseqs.txt"))
    structure_store = classes.ReferenceStructureStore(
        stream_package_data("structures.txt"), refseq_store=refseq_store)
    ranges_store = classes.ReportedRangeStore(
        stream_package_data("ranges_%s.txt" % arguments.ranges), structure_store=structure_store)
    for marker, sequence in (line.rstrip("\r\n").split() for line in arguments.instream):
        arguments.outstream.write("%s\t%s\n" % (marker, ranges_store.get_range(marker).get_name(sequence)))
        if arguments.unbuffered:
            arguments.outstream.flush()
#run


def run(arguments):
    if OPTIMIZE:
        optimize.run(arguments)
    else:
        run_strnaming(arguments)
#run


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
    else:
        parser.add_argument("-U", "--unbuffered", action="store_true",
            help="if specified, the output buffer is flushed after every line")
        parser.add_argument("-r", "--ranges", metavar="uas-frr", choices=("uas-frr",),
            required=True,
            help="specify 'uas-frr' to use UAS flanking region report ranges (mandatory; future "
                 "versions of STRNaming will offer more options)")
        parser.add_argument("instream", metavar="IN", nargs="?", type=argparse.FileType("tr"),
            default="-",
            help="an input file or stream containing a marker name and a raw sequence on each line, "
                 "separated by whitespace, which will be used to process in STRNaming (default: read "
                 "from stdin)")
        parser.add_argument("outstream", metavar="OUT", nargs="?", type=argparse.FileType("tw"),
            default=sys.stdout,
            help="output file, always uses a tab character as separator (default: write to stdout)")

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
        if args.debug:
            import cProfile
            cProfile.runctx("run(args)", globals(), locals(), sort="tottime")
        else:
            run(args)
    except Exception as error:
        if args.debug:
            raise
        parser.error(error)
#main


if __name__ == "__main__":
    main()
