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

from . import length_adjustments, libsequence, libstrnaming, refseq_analysis, refseq_cache, reference_structures, html
# TODO: Proper argument/input checking for all classes/methods!


class ReferenceSequenceStore:
    def __init__(self, *, autoload=True):
        """
        Construct a new ReferenceSequenceStore.

        If autoload is set to True (the default), the get_refseq()
        method will load the sequence from disk cache or Ensembl by
        default if necessary.
        """
        self.autoload = autoload
        self.refseqs = {}

    def load_from_cache(self, chromosome, start, end):
        """
        Load a portion of reference sequence from the refseq cache.

        The end position is exclusive.
        If the requested sequence does not reside in the refseq cache on
        disk, it is automatically downloaded from the Ensembl REST API.
        """
        self.add_refseq(chromosome, start, refseq_cache.get_refseq(chromosome, start, end - 1))

    def add_refseq(self, chromosome, start, seq):
        """
        Store a portion of reference sequence.
        """
        end = start + len(seq)  # Exclusive.
        if chromosome not in self.refseqs:
            self.refseqs[chromosome] = [[start, end, seq]]
        else:
            parts = self.refseqs[chromosome]
            i = 0
            while i <= len(parts):
                if i == len(parts) or parts[i][0] > start:
                    # New part starts before part i.
                    if i and parts[i-1][1] >= start:
                        if parts[i-1][1] < end:
                            # New part extends part i-1.
                            parts[i-1][2] += seq[parts[i-1][1] - start:]
                            parts[i-1][1] = end
                    else:
                        # New part comes between part i-1 and part i.
                        parts.insert(i, [start, end, seq])
                        i += 1  # Otherwise, i suddenly points to the new part.

                    while i < len(parts) and parts[i][0] <= end:
                        if parts[i][1] > end:
                            # New part (now stored as part i-1) extends into part i.
                            parts[i-1][2] += parts[i][2][parts[i-1][1] - parts[i][0]:]
                            parts[i-1][1] = parts[i][1]
                        # Part i is now completely covered by the new part.
                        del parts[i]
                    break  # Successfully inserted new sequence.
                i += 1

    def get_refseq(self, chromosome, start, end, *, autoload=None):
        """
        Get a portion of reference sequence.

        The end position is exclusive.
        If autoload=True, the sequence is loaded from disk cache or
        downloaded from Ensembl if necessary.
        """
        try:
            try:
                for part_start, part_end, seq in self.refseqs[chromosome]:
                    if part_start <= start <= end <= part_end:
                        return seq[start - part_start : end - part_start]
                raise ValueError(
                    "Refseq range %i..%i of chromosome %s not available" % (start, end - 1, chromosome))
            except KeyError:
                raise ValueError("No refseq for chromosome %s available" %  chromosome)
        except ValueError:
            if autoload is not False and (self.autoload or autoload):
                self.load_from_cache(chromosome, start, end)
                return self.get_refseq(chromosome, start, end)
            raise


class ReferenceStructureStore:
    def __init__(self, struct_input=None, *, refseq_store=None, analyse_refseqs=False):
        """
        Construct a new ReferenceStructureStore.

        If refseq_store is provided, it is used as the backing store of
        reference sequene data.
        If analyse_refseqs=True, the find_in_refseq() method is called for
        each portion of reference sequence present in refseq_store.
        If struct_input is provided, it is passed to load_from_tsv().
        """
        self.structures = {}
        if refseq_store is None:
            self.refseq_store = ReferenceSequenceStore()
        else:
            self.refseq_store = refseq_store
        if analyse_refseqs:
            for chromosome, segments in self.refseq_store.refseqs.items():
                for start, end, seq in segments:
                    self.find_in_refseq(chromosome, start, end)
        if struct_input is not None:
            self.load_from_tsv(struct_input)

    def load_within_range(self, chromosome, start, end):
        """
        Load human genome reference STR structure data.

        The end position is exclusive.
        """
        for structure in reference_structures.gen_within_range(chromosome, start, end - 1):
            self.add_structure(chromosome, structure)

    def load_within_ranges(self, chromosome, ranges):
        """
        Load human genome reference STR structure data.

        The ranges should be a sorted iterable of (start, end) pairs.
        The end positions are exclusive.
        """
        for structure in reference_structures.gen_within_ranges(
                chromosome, ((start, end - 1) for start, end in ranges)):
            self.add_structure(chromosome, structure)

    def load_from_tsv(self, input):
        """
        Load custom reference STR structure data from an input stream.

        The input is a stream with one STR structure definition per line.
        An STR structure definition consists of a chromosome number,
        followed by one or more stretch definitions, separated by whitespace.
        A stretch definition consists of a start position, end position
        (exclusive), and repeat unit sequence, separated by whitespace.
        """
        for line in input:
            chromosome, *stretches = line.rstrip("\r\n").split()
            structure = []
            for i in range(0, len(stretches), 3):
                structure.append([int(stretches[i]), int(stretches[i+1]), len(stretches[i+2])])
            self.add_structure(chromosome, structure)

    def find_in_refseq(self, chromosome, start, end, autoload=False):
        """
        Automatically detect STR structures in a portion of reference sequence.

        The end position is exclusive.
        If autoload=True, the human genome reference sequence for this range is
        loaded from disk cache or downloaded from Ensembl if necessary.
        """
        seq = self.refseq_store.get_refseq(chromosome, start, end, autoload=autoload)
        for structure in refseq_analysis.recurse_collapse_repeat_units_refseq(seq, offset=start):
            self.add_structure(
                chromosome, [[s_start, s_end, len(unit)] for s_start, s_end, unit in structure])

    def add_structure(self, chromosome, stretches):
        """
        Store an STR structure.

        The stretches consist of lists (or tuples) containing the start
        position, end position (exclusive) and repeat unit length.
        """
        start = stretches[0][0]
        end = stretches[-1][1]
        if chromosome not in self.structures:
            self.structures[chromosome] = [tuple(stretches)]
        else:
            structures = self.structures[chromosome]
            for i, structure in enumerate(structures):
                if start == structure[0][0] and end == structure[-1][1]:
                    break  # Already have this structure.
                if end >= structure[0][0] and start <= structure[-1][1]:
                    raise ValueError("Reference structures cannot touch or overlap")
                if start < structure[0][0]:
                    # New structure goes before structure i.
                    structures.insert(i, tuple(stretches))
                    break
            else:
                # Structure goes at the end of the list.
                structures.append(tuple(stretches))

    def get_structures(self, chromosome, start, end):
        """
        Get a list of STR structures on the given portion of the genome.

        The end position is exclusive.
        The structures must have been detected by prior calls to
        load_from_tsv(), find_in_refseq(), or add_structure().
        If no known structures exist at the given portion of the genome,
        an empty list is returned.
        """
        structures = []
        try:
            for structure in self.structures[chromosome]:
                if end > structure[0][0] and start < structure[-1][1]:
                    structures.append(structure)
        except KeyError:
            pass  # No structures for this chromosome, that's fine.
        return structures

    def get_refseq_store(self):
        return self.refseq_store


class ReportedRangeStore:
    def __init__(self, ranges_input=None, *, structure_store=None, load_structures=False):
        """
        Construct a new ReportedRangeStore.

        If structure_store is provided, it is used as the backing store
        of reference STR structure and reference sequence data.
        If ranges_input is provided, it is passed to load_from_tsv(),
        together with the load_structures argument if given.
        """
        self.ranges = {}
        self.structure_store = structure_store if structure_store else ReferenceStructureStore()
        if ranges_input is not None:
            self.load_from_tsv(ranges_input, load_structures=load_structures)

    def load_from_tsv(self, input, *, load_structures=False):
        """
        Load reported ranges from a tab-separated input stream.

        The first line contains headers, which minimally includes: 'name',
        'chromosome', 'start', 'end', and 'options'.

        The end position is inclusive.  If sequences are expected in reverse
        complement orientation, the start and end positions should be swapped.

        If load_structures=True, ReferenceStructureStore.load_within_range() is
        called prior to configuring each ReportedRange.
        """
        headers = input.readline().rstrip("\r\n").split("\t")
        try:
            ci = {name: headers.index(name) for name in
                ("name", "chromosome", "start", "end", "options")}
        except ValueError:
            raise ValueError("The reported ranges file should contain tab-separated columns "
                "named 'name', 'chromosome', 'start', 'end', and 'options'.")
        PAT_OPTIONS = re.compile("([^=,]+)=([^=,]+)")
        for line in input:
            values = line.rstrip("\r\n").split("\t")
            options = dict(PAT_OPTIONS.findall(values[ci["options"]]))
            start = int(values[ci["start"]])
            end = int(values[ci["end"]])
            if end < start:
                start, end = end, start
                options["reverse_complement"] = True
            self.add_range(values[ci["name"]], values[ci["chromosome"]], start,
                end + 1, load_structures=load_structures, options=options)

    def add_complex_range(self, name, genome_position, *, options=None):
        """
        Store a new reported range and return its ReportedRange object.

        The end positions in genome_position are inclusive.
        """
        if name in self.ranges:
            raise ValueError("Duplicate fragment name '%s'" % name)
        self.ranges[name] = ComplexReportedRange(
            self.structure_store, genome_position, options=options)
        return self.ranges[name]

    def add_range(self, name, chromosome, start, end, *, load_structures=False, options=None):
        """
        Store a new reported range and return its ReportedRange object.

        The end position is exclusive.
        If load_structures=True, ReferenceStructureStore.load_within_range() is
        called prior to configuring the ReportedRange.
        """
        if name in self.ranges:
            raise ValueError("Duplicate fragment name '%s'" % name)
        if load_structures:
            self.structure_store.load_within_range(chromosome, start, end)
        self.ranges[name] = ReportedRange(
            self.structure_store, chromosome, start, end, options=options)
        return self.ranges[name]

    def get_range(self, name):
        """
        Get a description of the reported range with the given name.
        """
        try:
            return self.ranges[name]
        except KeyError:
            raise ValueError("Range '%s' not found" % name)

    def get_ranges(self):
        return self.ranges;

    def get_structure_store(self):
        return self.structure_store


class ComplexReportedRange:
    def __init__(self, structure_store, genome_position, *, options=None):
        if not len(genome_position) % 2 or len(genome_position) < 3:
            raise ValueError(
                "Genome position %r invalid; it must consist of a chromosome and one or more "
                "pairs of start and end positions." % (genome_position,))
        self.options = {} if options is None else {key: value for key, value in options.items()}
        self.location = tuple(genome_position)
        self.refseq = ""
        refseq_store = structure_store.get_refseq_store()
        chromosome = genome_position[0]
        for i in range(1, len(genome_position), 2):
            start, end = genome_position[i : i + 2]
            if structure_store.get_structures(chromosome, start, end + 1):
                raise ValueError(
                    "Complex range overlaps a known STR structure in range Chr%s:%i..%i" %
                    (chromosome, start, end))
            self.refseq += refseq_store.get_refseq(chromosome, start, end + 1)

    def get_option(self, option, default=None):
        return self.options.get(option, default)

    def has_option(self, option):
        return option in self.options

    def get_tssv(self, seq, *, as_string=True):
        return seq + "(1)" if as_string else [[seq, 1, None]]

    def from_tssv(self, tssv):
        return tssv[:-3]

    def get_name(self, seq):
        variants = libsequence.call_variants(self.refseq, seq, location=self.location)
        return " ".join(variants) or "REF"

    def from_name(self, name):
        """Convert allele name to a raw sequence."""
        if name == "REF":
            return self.refseq
        return libsequence.mutate_sequence(self.refseq, name.split(), location=self.location)


class ReportedRange:  # TODO: this could extend ComplexReportedRange to avoid code duplication.
    def __init__(self, structure_store, chromosome, start, end, *, options=None):
        refseq_store = structure_store.get_refseq_store()
        longest_stretch = 1

        self.options = {} if options is None else {key: value for key, value in options.items()}
        self.location = (chromosome, start, end - 1)
        self.limit = int(self.options.get("limit", 0))
        self.library = []
        self.block_length = 1
        self.length_adjust = 0
        self.preinsert = ""
        self.postinsert = ""
        self.reverse_complement = self.options.get("reverse_complement", False)

        # Analyse sequence to extract repeat units.
        first_structure_start = None
        pos = start
        for structure in structure_store.get_structures(chromosome, start, end):

            # Drop any repeat stretches that have more than an entire repeat
            # outside either end of the reported range.
            filtered_structure = list(filter(lambda stretch: (
                    start <= min(stretch[0] + stretch[2], stretch[1] - 1) and
                    end >= max(stretch[1] - stretch[2], stretch[0] + 1)),
                structure))

            # If we now end up with too short an STR structure, drop it.
            if not filtered_structure or filtered_structure[-1][1] - filtered_structure[0][0] < libstrnaming.NAMING_OPTIONS["min_structure_length"]:
                continue

            # We need to remember where the first actually included structure began, as
            # we don't want to apply the length_adjustments of those that we dropped.
            if first_structure_start is None:
                first_structure_start = structure[0][0]

            # Set length adjustment according to prefix/infix/preinsert.
            self.length_adjust += pos - structure[0][0]

            # Store end position of last structure for later suffix/postinsert adjustment.
            last_structure_end = structure[-1][1]

            # Determine the block_length of the full reference structure.
            block_length = libstrnaming.find_block_length(structure)

            # Find unique repeat units.
            units = [refseq_store.get_refseq(chromosome, stretch[0], stretch[0] + stretch[2])
                for stretch in sorted(structure,
                    key=lambda stretch: (stretch[1]-stretch[0], stretch[2]), reverse=True)]
            units = [unit for index, unit in enumerate(units) if units.index(unit) == index]

            # From here on, we are not interested in the removed stretches anymore.
            structure = filtered_structure

            # Construct a regex pattern.
            regex_blocks = []
            regex_pos = pos
            for stretch_start, stretch_end, unit_length in structure:
                if stretch_start > regex_pos:
                    regex_blocks.append(
                        "%s" % refseq_store.get_refseq(chromosome, regex_pos, stretch_start))
                regex_blocks.append("((%s){%i,})" % (refseq_store.get_refseq(
                        chromosome, stretch_start, stretch_start + unit_length),
                        3 if unit_length == 1 else 2 if unit_length == 2 else 1))
                regex_pos = stretch_end
            if regex_pos < end:
                regex_blocks.append("%s" % refseq_store.get_refseq(chromosome, regex_pos, end))
            regex = re.compile("^" + "".join(regex_blocks) + "$")

            # Store preinsert/postinsert if the structure extends slightly outside range.
            if start > structure[0][0]:
                self.preinsert = refseq_store.get_refseq(chromosome, structure[0][0], start)
            if end < structure[-1][1]:
                self.postinsert = refseq_store.get_refseq(chromosome, end, structure[-1][1])

            # Update block_length if this structure contains the longest stretch so far.
            # NOTE: The ReportedRange's block_length is determined by in-range stretches only!
            for stretch_start, stretch_end, unit_length in structure:
                stretch_length = stretch_end - stretch_start
                if stretch_length < longest_stretch:
                    continue
                if stretch_length == longest_stretch and unit_length < self.block_length:
                    continue
                longest_stretch = stretch_length
                self.block_length = unit_length

            prefix = refseq_store.get_refseq(
                chromosome, pos, structure[0][0]) if pos < structure[0][0] else ""
            if self.library:
                self.library[-1]["suffix"] = prefix
            pos = structure[-1][1]

            overlong_gap = libstrnaming.find_overlong_gap(structure)

            self.library.append({
                # Sequence included in the range before/after the STR structure.
                "prefix": prefix,
                "suffix": "",

                # Sequence of the largest interruption within the STR structure.
                "overlong_gap": refseq_store.get_refseq(
                    chromosome, overlong_gap[0], overlong_gap[1]) if overlong_gap else "",

                # The repeat units in the reference STR structure.
                "units": units,

                # The unit length of the longest repeat in the reference STR structure.
                "block_length": block_length,

                # A compiled regular expression that matches the reference structure.
                "regex": regex
            })

        if self.library:
            if pos < end:
                # Save suffix of last repeat structure in the library.
                self.library[-1]["suffix"] = refseq_store.get_refseq(chromosome, pos, end)
            self.length_adjust += last_structure_end - end + length_adjustments.get_adjustment(
                chromosome, first_structure_start, last_structure_end, len(self.library) > 1)
        else:
            self.refseq = refseq_store.get_refseq(chromosome, start, end)
            self.length_adjust = None

    def get_option(self, option, default=None):
        return self.options.get(option, default)

    def has_option(self, option):
        return option in self.options

    def get_ce(self, seq):
        if self.limit and len(seq) == self.limit:
            return "?"
        length = len(seq) + self.length_adjust
        sign = "-" if length < 0 else ""
        length = abs(length)
        return "%s%s%s" % (sign, length // self.block_length,
            "." + str(length % self.block_length) if length % self.block_length else "")

    def normalize_sequence(self, seq, *, lowercase_inserts=False):
        if self.reverse_complement:
            seq = libsequence.reverse_complement(seq)
        if lowercase_inserts:
            return self.preinsert.lower() + seq + self.postinsert.lower()
        return self.preinsert + seq + self.postinsert

    def denormalize_sequence(self, seq):
        if self.preinsert and not seq.startswith(self.preinsert):
            raise ValueError("Sequene should start with " + self.preinsert)
        if self.postinsert and not seq.endswith(self.postinsert):
            raise ValueError("Sequene should end with " + self.postinsert)
        seq = seq[len(self.preinsert) : -len(self.postinsert) if self.postinsert else len(seq)]
        return libsequence.reverse_complement(seq) if self.reverse_complement else seq

    def get_stretches(self, seq, *, normalized_seq=None):
        # Short-circuit if no STR structures are within range.
        if not self.library:
            return []
        if normalized_seq is None:
            normalized_seq = self.normalize_sequence(seq)
        stretches = []
        pos = 0
        for i, part in enumerate(self.library):
            suffix = part["suffix"]
            if i < len(self.library) - 1:
                # Find end position of this infix.
                # If suffix in part["prefix"], look beyond it!
                skip = part["prefix"].rfind(suffix) + 1
                pos += skip
                try:
                    end = normalized_seq.index(suffix, pos) + len(suffix)
                except ValueError:
                    suffix_start, suffix_end, score = libsequence.align(normalized_seq[pos:], suffix)
                    end = pos + suffix_end
                    suffix = normalized_seq[pos + suffix_start : end]
                pos -= skip
            else:
                end = len(normalized_seq)
            match = part["regex"].match(normalized_seq[pos:end])
            if match:
                path = [match.span(group) + (match.group(group + 1),)
                        for group in range(1, part["regex"].groups, 2)]
            else:
                try:
                    path = libstrnaming.collapse_repeat_units(normalized_seq[pos:end],
                        part["prefix"], suffix, part["units"], part["block_length"], len(part["overlong_gap"]))
                except libstrnaming.OutOfTimeException:
                    # TODO: Add range name to this message (and possibly others).
                    sys.stderr.write("STRNaming Ran out of time while analysing sequence %s\n"
                        % (normalized_seq[pos:end],))
                    path = []
            stretches += [[s_start + pos, s_end + pos, unit, i] for s_start, s_end, unit in path]
            pos = stretches[-1][1] if path else end
            if pos >= len(normalized_seq):
                # If we have repeats right up to the end,
                # we can't make any more structures...
                break
        return stretches

    def color_template_for(self, unit, library_part, *,
            num_classes=html.NUM_CLASSES, span_template=html.SPAN_TEMPLATE):
        try:
            return span_template % (self.library[library_part]["units"].index(unit) % num_classes)
        except ValueError:
            return "%s"  # Not a reference unit.

    def colorize_sequence(self, seq, *, normalized_seq=None, stretches=None):
        if normalized_seq is None:
            normalized_seq = self.normalize_sequence(seq, lowercase_inserts=True)
        if not self.library:
            return normalized_seq
        if stretches is None:
            stretches = self.get_stretches(seq, normalized_seq=normalized_seq.upper())
        pos = 0
        colored = []
        for start, end, unit, i in stretches:
            if start > pos:
                # Insert a gap or prefix.
                colored.append(normalized_seq[pos : start])
            colored.append(self.color_template_for(unit, i) % normalized_seq[start : end])
            pos = end
        if pos < len(normalized_seq):
            # Add the suffix.
            colored.append(normalized_seq[pos:])
        return "".join(colored)

    def get_tssv(self, seq, *, as_string=True, normalized_seq=None, stretches=None):
        if normalized_seq is None:
            normalized_seq = self.normalize_sequence(seq)
        if not self.library:
            return normalized_seq + "(1)" if as_string else [[normalized_seq, 1, None]]
        if stretches is None:
            stretches = self.get_stretches(seq, normalized_seq=normalized_seq)
        pos = 0
        tssv_seq = []
        for start, end, unit, i in stretches:
            if start > pos:
                # Insert a gap or prefix.
                tssv_seq.append([normalized_seq[pos : start], 1, i])
            tssv_seq.append([unit, (end - start) // len(unit), i])
            pos = end
        if pos < len(normalized_seq):
            # Add the suffix.
            tssv_seq.append([normalized_seq[pos:], 1, len(self.library)-1])
        if as_string:
            return "".join("%s(%i)" % (unit, count) for unit, count, i in tssv_seq)
        return tssv_seq

    def from_tssv(self, tssv):
        """Convert TSSV-style sequence to a raw sequence."""
        return self.denormalize_sequence("".join(block[0] * int(block[1])
            for block in libsequence.PAT_TSSV_BLOCK.findall(tssv)))

    def get_name(self, seq, *, color=False, normalized_seq=None, stretches=None):
        if self.limit and len(seq) == self.limit:
            # TODO: Give sequence-descriptive name nonetheless!
            return "CE?_TODO_UAS_INCOMPLETE_SEQUENCE"
        if normalized_seq is None:
            normalized_seq = self.normalize_sequence(seq)

        if not self.library:
            # Report variants if there are not STRs on this range.
            variants = libsequence.call_variants(
                self.refseq, normalized_seq, location=self.location)
            return " ".join(variants) or "REF"

        if stretches is None:
            stretches = self.get_stretches(seq, normalized_seq=normalized_seq)

        prefix = self.library[0]["prefix"]
        suffix = self.library[-1]["suffix"]
        blocks = []
        variants = []

        pos = stretches[0][0] if stretches else len(normalized_seq)
        before = normalized_seq[:pos]
        if prefix:
            if not before:
                # Deleted prefix completely.
                variants.append("-%s%s>-" % (len(prefix), prefix))
            elif before.startswith(prefix):
                if pos > len(prefix):
                    if stretches:
                        # Insertion at end of prefix only.
                        blocks.append(before[len(prefix):] + "[1]")
                    else:
                        # Remainder of sequence could be suffix.
                        pos = len(prefix)
            else:
                # Modified prefix.
                prefix_variants = libsequence.call_variants(prefix, before, location="prefix")
                if prefix_variants[0].startswith("-0."):
                    if stretches:
                        # Rename next-to-STR insertion to an unrepeated unit.
                        blocks.append(prefix_variants[0][6:] + "[1]")
                    else:
                        # Remainder of sequence could be suffix.
                        pos -= len(prefix_variants[0]) - 6
                    prefix_variants = prefix_variants[1:]
                variants += prefix_variants
        elif before:
            if stretches:
                # No prefix, but there is an insertion at the start of the fragment.
                blocks.append(before + "[1]")
            else:
                # All of the sequence could be suffix.
                pos = 0

        prev_i = 0
        for start, end, unit, i in stretches:
            if start > pos:
                gap = normalized_seq[pos:start]
                if start - pos > libstrnaming.NAMING_OPTIONS["max_gap"] and (i > prev_i or self.library[i]["overlong_gap"]):
                    # Handle overlong gap.
                    gap_ref = self.library[i]["prefix" if i > prev_i else "overlong_gap"]
                    if gap_ref == gap:
                        blocks.append("[]")
                    else:
                        this_variants = libsequence.call_variants(gap_ref, gap, location="suffix")
                        if this_variants[0].startswith("+0."):
                            # Rename initial insertion to an unrepeated unit.
                            blocks.append(this_variants[0][6:] + "[1]")
                            this_variants = this_variants[1:]
                        # TODO: Move insertions at end into blocks too.
                        this_variants = " ".join(x[1:] for x in this_variants)
                        if len(this_variants) < len(gap):
                            blocks.append("[" + this_variants + "]")
                        else:
                            blocks.append(gap + "[1]")
                else:
                    blocks.append(gap + "[1]")
            block = "%s[%i]" % (unit, (end - start) // len(unit))
            blocks.append(self.color_template_for(unit, i) % block if color else block)
            prev_i = i
            pos = end

        after = normalized_seq[pos:]
        if suffix:
            if not after:
                # Deleted suffix completely.
                variants.append("+1" + suffix + ">-")
            elif after.endswith(suffix):
                if len(after) > len(suffix):
                    # Insertion at start of suffix only.
                    blocks.append(after[:-len(suffix)] + "[1]")
            else:
                # Modified suffix.
                suffix_variants = libsequence.call_variants(suffix, after, location="suffix")
                if suffix_variants[0].startswith("+0."):
                    # Rename next-to-STR insertion to an unrepeated unit.
                    blocks.append(suffix_variants[0][6:] + "[1]")
                    suffix_variants = suffix_variants[1:]
                variants += suffix_variants
        elif after:
            # No suffix, remaining sequence is insertion at end of fragment.
            blocks.append(after + "[1]")

        # Construct allele name.
        parts = ["CE" + self.get_ce(seq), "".join(blocks)]
        if variants:
            parts.append("_".join(variants))
        return "_".join(parts)

    def from_name(self, name):
        """Convert allele name to a raw sequence."""
        if name == "CE?_TODO_UAS_INCOMPLETE_SEQUENCE":
            raise ValueError("Cannot reconstruct incomplete UAS sequence")
        if not self.library:
            return self.denormalize_sequence(self.refseq if name == "REF" else
                libsequence.mutate_sequence(self.refseq, name.split(), location=self.location))

        allele = name.split("_")
        if len(allele) < 2:
            raise ValueError("Invalid allele name '%s' encountered!" % name)

        # Get and mutate prefix and suffix.
        prefix = self.library[0]["prefix"]
        suffix = self.library[-1]["suffix"]
        variants = [[], []]
        for variant in allele[2:]:
            if variant[0] == "-":
                variants[0].append(variant)
            elif variant[0] == "+":
                variants[1].append(variant)
            else:
                raise ValueError("Unrecognised variant '%s'" % variant)
        if variants[0]:
            if not prefix:
                raise ValueError(
                    "Encountered prefix variants %r, but range has no prefix!" % (variants[0],))
            prefix = libsequence.mutate_sequence(prefix, variants[0])
        if variants[1]:
            if not suffix:
                raise ValueError(
                    "Encountered suffix variants %r, but range has no suffix!" % (variants[1],))
            suffix = libsequence.mutate_sequence(suffix, variants[1])

        # See how many overlong gaps we're dealing with.
        overlong_refs = [seq for part_i, part in enumerate(self.library) for seq in
                (part["prefix"] if part_i else "", part["overlong_gap"]) if seq]
        num_overlong_gaps = len(libsequence.PAT_ALLELENAME_GAP.findall(name))
        max_overlong_gaps = len(overlong_refs)
        if max_overlong_gaps < num_overlong_gaps:
            raise ValueError("Invalid allele name for this range; too many overlong gaps!")
        if num_overlong_gaps and max_overlong_gaps != num_overlong_gaps:
            raise ValueError("Cannot reconstruct sequence; skipping only some overlong gaps "
                "is not supported yet")

        # Reconstruct the sequence.
        seq = [prefix]
        overlong_i = 0
        for i, name_part in enumerate(libsequence.PAT_ALLELENAME_GAP.split(allele[1])):
            if i % 2:
                # Overlong gap.
                seq.append(overlong_refs[overlong_i] if not name_part else
                    libsequence.mutate_sequence(
                        overlong_refs[overlong_i], name_part.split(), location=(None, 1)))
                overlong_i += 1
            else:
                # Explicitly visible STR structure.
                for block in libsequence.PAT_ALLELENAME_BLOCK.finditer(name_part):
                    seq.append(block.group(1) * int(block.group(2)))
        seq.append(suffix)
        return self.denormalize_sequence("".join(seq))
