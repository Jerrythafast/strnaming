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
import re

from . import libstrnaming
# TODO: Proper argument/input checking for all classes/methods!


class ReferenceSequenceStore:
    def __init__(self, input=None, *, online=False):
        """
        Construct a new ReferenceSequenceStore.

        If input is given, it is passed to load_from_tsv().
        If online is set to True, the get_refseq() function will
        download the sequence from Ensembl by default if necessary.
        """
        self.online = online
        self.refseqs = {}
        if input is not None:
            self.load_from_tsv(input)

    def load_from_tsv(self, input):
        """
        Load reference sequences from a tab-separated input stream.

        The first line contains headers, which minimally includes:
        'chromosome', 'start', and  sequence'.
        """
        headers = input.readline().rstrip("\r\n").split("\t")
        try:
            ci = {name: headers.index(name) for name in ("chromosome", "start", "sequence")}
        except ValueError:
            raise ValueError("The genome reference file should contain tab-separated columns "
                "named 'chromosome', 'start', and 'sequence'.")
        for line in input:
            values = line.rstrip("\r\n").split("\t")
            self.add_refseq(
                values[ci["chromosome"]], int(values[ci["start"]]), values[ci["sequence"]])

    def load_from_ensembl(self, chromosome, start, end):
        """
        Load a portion of reference sequence from Ensembl.

        The end position is exclusive.
        """
        # TODO: Error handling.
        import urllib.request
        self.add_refseq(chromosome, start, urllib.request.urlopen(
            "http://rest.ensembl.org/sequence/region/human/%s:%i..%i?content-type=text/plain" % (
                chromosome, start, end - 1)).read().decode("UTF-8"))

    def add_refseq(self, chromosome, start, seq):
        """
        Store a portion of reference sequence.
        """
        end = start + len(seq)  # Exclusive.
        if chromosome not in self.refseqs:
            self.refseqs[chromosome] = [(start, end, seq)]
        else:
            parts = self.refseqs[chromosome]
            # TODO: handle overlap with / extension of existing parts.
            for i, part in enumerate(parts):
                if start < part[0]:
                    if i and start <= parts[i-1][1]:
                        raise ValueError("Range %i-%i of chromosome %s requires merger with range %i-%i, which is currently unsupported." % (start, end, chromosome, parts[i-1][0], parts[i-1][1]))
                    if end >= part[0]:
                        raise ValueError("Range %i-%i of chromosome %s requires merger with range %i-%i, which is currently unsupported." % (start, end, chromosome, part[0], part[1]))
                    parts.insert(i, (start, end, seq))
                    break
            else:
                if parts and start <= parts[-1][1]:
                    raise ValueError("Range %i-%i of chromosome %s requires merger with range %i-%i, which is currently unsupported." % (start, end, chromosome, parts[-1][0], parts[-1][1]))
                parts.append((start, end, seq))

    def get_refseq(self, chromosome, start, end, *, online=None):
        """
        Get a portion of reference sequence.

        The end position is exclusive.
        If online=True, the sequence is downloaded from Ensembl if necessary.
        """
        try:
            try:
                for part_start, part_end, seq in self.refseqs[chromosome]:
                    if part_start <= start <= end <= part_end:
                        return seq[start - part_start : end - part_start]
                raise ValueError(
                    "Refseq range %i-%i of chromosome %s not available" % (start, end, chromosome))
            except KeyError:
                raise ValueError("No refseq for chromosome %s available" %  chromosome)
        except ValueError:
            if online is not False and (self.online or online):
                self.load_from_ensembl(chromosome, start, end)
                return self.get_refseq(chromosome, start, end)
            raise


class ReferenceStructureStore:
    def __init__(self, struct_input=None, *, refseq_store=None, refseq_input=None, analyse_refseqs=False):
        """
        Construct a new ReferenceStructureStore.

        If refseq_store is provided, it is used as the backing store of
        reference sequene data.
        If refseq_input is provided, reference sequence data is loaded
        from it using refseq_store.load_from_tsv(refseq_input).
        If analyse_refseqs=True, the find_in_refseq() method is called for
        each portion of reference sequence present in refseq_store and refseq_input.
        If struct_input is provided, it is passed to load_from_tsv().
        """
        self.structures = {}
        if refseq_store is None:
            self.refseq_store = ReferenceSequenceStore(refseq_input)
        else:
            self.refseq_store = refseq_store
            if refseq_input is not None:
                refseq_store.load_from_tsv(refseq_input)
        if analyse_refseqs:
            for chromosome, segments in self.refseq_store.refseqs.items():
                for start, end, seq in segments:
                    self.find_in_refseq(chromosome, start, end)
        if struct_input is not None:
            self.load_from_tsv(struct_input)

    def load_from_tsv(self, input):
        """
        Load reference STR structure data.

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
                structure.append([int(stretches[i]), int(stretches[i+1]), stretches[i+2]])
            self.add_structure(chromosome, structure)

    def find_in_refseq(self, chromosome, start, end, online=False):
        """
        Automatically detect STR structures in a portion of reference sequence.

        The end position is exclusive.
        If online=True, the sequence is downloaded from Ensembl if necessary.
        """
        seq = self.refseq_store.get_refseq(chromosome, start, end, online=online)
        for structure in libstrnaming.recurse_collapse_repeat_units_refseq(seq, offset=start):
            self.add_structure(chromosome, structure)

    def add_structure(self, chromosome, stretches):
        """
        Store an STR structure.

        The stretches consist of lists (or tuples) containing the start
        position, end position (exclusive) and repeat unit sequence.
        """
        start = stretches[0][0]
        end = stretches[-1][1]
        if chromosome not in self.structures:
            self.structures[chromosome] = [tuple(stretches)]
        else:
            structures = self.structures[chromosome]
            for i, structure in enumerate(structures):
                if start < structure[0][0]:
                    structures.insert(i, tuple(stretches))
                    break
            else:
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
    def __init__(self, ranges_input=None, *, structure_store=None):
        """
        Construct a new ReportedRangeStore.

        If structure_store is provided, it is used as the backing store
        of reference STR structure and reference sequence data.
        If ranges_input is provided, it is passed to load_from_tsv().
        """
        self.ranges = {}
        self.structure_store = structure_store if structure_store else ReferenceStructureStore()
        if ranges_input is not None:
            self.load_from_tsv(ranges_input)

    def load_from_tsv(self, input):
        """
        Load reported ranges from a tab-separated input stream.

        The first line contains headers, which minimally includes: 'name',
        'chromosome', 'start', 'end', 'ref_ce', and 'options'.

        The end position is inclusive.  If sequences are expected in reverse
        complement orientation, the start and end positions should be swapped.
        """
        headers = input.readline().rstrip("\r\n").split("\t")
        try:
            ci = {name: headers.index(name) for name in
                ("name", "chromosome", "start", "end", "ref_ce", "options")}
        except ValueError:
            raise ValueError("The reported ranges file should contain tab-separated columns "
                "named 'name', 'chromosome', 'start', 'end', 'ref_ce', and 'options'.")
        PAT_OPTIONS = re.compile("([^=,]+)=([^=,]+)")
        for line in input:
            values = line.rstrip("\r\n").split("\t")
            options = dict(PAT_OPTIONS.findall(values[ci["options"]]))
            reverse_complement = False
            start = int(values[ci["start"]])
            end = int(values[ci["end"]])
            if end < start:
                start, end = end, start
                reverse_complement = True
            self.add_range(values[ci["name"]], values[ci["chromosome"]], start,
                end + 1, values[ci["ref_ce"]], reverse_complement, options)

    def add_range(self, name, chromosome, start, end, ref_ce, reverse_complement, options):
        """
        Store a reported range.

        The end position is exclusive.
        The reference CE number can be an empty string if unknown.
        The options are passed to the ReportedRange constructor. If the "aka"
        option is present, this reported range will be made available under
        that name as well.
        """
        if name in self.ranges:
            raise ValueError("Duplicate fragment name '%s'" % name)
        self.ranges[name] = ReportedRange(self.structure_store, chromosome, start, end, ref_ce, reverse_complement, options)
        if "aka" in options:
            if options["aka"] in self.ranges:
                raise ValueError("Duplicate fragment name '%s'" % options["aka"])
            self.ranges[options["aka"]] = self.ranges[name]

    def get_range(self, name):
        """
        Get a description of the reported range with the given name.
        """
        try:
            return self.ranges[name]
        except KeyError:
            raise ValueError("Range '%s' not found" % name)


class ReportedRange:
    def __init__(self, structure_store, chromosome, start, end, ref_ce, reverse_complement, options):
        match_ce = libstrnaming.PAT_CE.match(ref_ce)
        refseq_store = structure_store.get_refseq_store()
        longest_stretch = 1

        self.limit = int(options.get("limit", 0))
        self.library = []
        self.block_length = 1
        self.length_adjust = 0
        self.preinsert = ""
        self.postinsert = ""
        self.reverse_complement = reverse_complement

        # Analyse sequence to extract repeat units.
        pos = start
        for structure in structure_store.get_structures(chromosome, start, end):

            # Drop any repeat stretches that have more than an entire repeat
            # outside either end of the reported range.
            structure = list(filter(lambda stretch: (
                    start <= min(stretch[0] + len(stretch[2]), stretch[1] - 1) and
                    end >= max(stretch[1] - len(stretch[2]), stretch[0] + 1)),
                structure))

            # If we now end up with too short an STR structure, drop it.
            if not structure or structure[-1][1] - structure[0][0] < libstrnaming.NAMING_OPTIONS["min_structure_length"]:
                continue

            # Store preinsert/postinsert if the structure extends slightly outside range.
            if start > structure[0][0]:
                self.preinsert = refseq_store.get_refseq(chromosome, structure[0][0], start)
                self.length_adjust -= start - structure[0][0]
            if end < structure[-1][1]:
                self.postinsert = refseq_store.get_refseq(chromosome, end, structure[-1][1])
                self.length_adjust -= structure[-1][1] - end

            # Update block_length if this structure contains the longest stretch so far.
            for stretch_start, stretch_end, unit in structure:
                stretch_length = stretch_end - stretch_start
                if stretch_length < longest_stretch:
                    continue
                unit_length = len(unit)
                if stretch_length == longest_stretch and unit_length < self.block_length:
                    continue
                longest_stretch = stretch_length
                self.block_length = unit_length

            # Update length adjustment.
            if pos < structure[0][0]:
                self.length_adjust += structure[0][0] - pos
            if match_ce:
                self.length_adjust += structure[-1][1] - structure[0][0]

            prefix = refseq_store.get_refseq(chromosome, pos, structure[0][0]) if pos < structure[0][0] else ""
            if self.library:
                if not prefix:
                    # TODO: Block this from happening in StructureStore.
                    sys.stderr.write("BadLibrary=%r\nNextStructure=%r\n" % (library, structure))
                    raise ValueError("Unseparated structure encountered")
                self.library[-1]["suffix"] = prefix
            pos = structure[-1][1]

            overlong_gap = libstrnaming.find_overlong_gap(structure)

            self.library.append({
                # Sequence included in the range before/after the STR structure.
                "prefix": prefix,
                "suffix": "",

                # Sequence of the largest interruption within the STR structure.
                "overlong_gap": refseq_store.get_refseq(chromosome, overlong_gap[0], overlong_gap[1]) if overlong_gap else "",

                # The repeat units in the STR structure.
                "units":  [stretch[2] for stretch in sorted(structure,
                    key=lambda stretch: (stretch[1]-stretch[0], len(stretch[2])), reverse=True)]
            })

        if self.library:
            if pos < end:
                # Save suffix of last repeat structure in the library.
                self.library[-1]["suffix"] = refseq_store.get_refseq(chromosome, pos, end)
                self.length_adjust += end - pos
        else:
            self.refseq = refseq_store.get_refseq(chromosome, start, end)
            self.location = (chromosome, start)
            if match_ce:
                self.length_adjust += end - start
            else:
                self.length_adjust = None
        if match_ce:
            # Update length adjustment with user-provided CE allele number.
            self.length_adjust -= int(match_ce.group(1)) * self.block_length + int(match_ce.group(2) or 0)

    def get_ce(self, seq):
        if self.limit and len(seq) == self.limit:
            return "?"
        length = len(seq) - self.length_adjust
        sign = "-" if length < 0 else ""
        length = abs(length)
        return "%s%s%s" % (sign, length // self.block_length,
            "." + str(length % self.block_length) if length % self.block_length else "")

    def normalize_sequence(self, seq):
        if self.reverse_complement:
            seq = libstrnaming.reverse_complement(seq)
        return self.preinsert + seq + self.postinsert

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
                try:
                    end = normalized_seq.index(suffix, pos) + len(suffix)
                except ValueError:
                    suffix_start, suffix_end, score = libstrnaming.align(normalized_seq[pos:], suffix)
                    end = pos + suffix_end
                    suffix = normalized_seq[pos + suffix_start : end]
            else:
                end = len(normalized_seq)
            path = libstrnaming.collapse_repeat_units(normalized_seq[pos:end],
                part["prefix"], suffix, part["units"], part["overlong_gap"])
            stretches += [[s_start + pos, s_end + pos, unit, i] for s_start, s_end, unit in path]
            pos = stretches[-1][1] if path else end
            if pos >= len(normalized_seq):
                # If we have repeats right up to the end,
                # we can't make any more structures...
                break
        return stretches

    def get_tssv(self, seq, as_string=False, *, normalized_seq=None):
        if normalized_seq is None:
            normalized_seq = self.normalize_sequence(seq)
        if not self.library:
            return [[normalized_seq, 1, None]]
        pos = 0
        tssv_seq = []
        for start, end, unit, i in self.get_stretches(seq, normalized_seq=normalized_seq):
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

    def get_name(self, seq, *, normalized_seq=None):
        if self.limit and len(seq) == self.limit:
            # TODO: Give sequence-descriptive name nonetheless!
            return "CE?_TODO_UAS_INCOMPLETE_SEQUENCE"
        if normalized_seq is None:
            normalized_seq = self.normalize_sequence(seq)

        if not self.library:
            # Report variants if there are not STRs on this range.
            variants = libstrnaming.call_variants(
                self.refseq, normalized_seq, location=self.location)
            if self.length_adjust is None:
                return " ".join(variants) or "REF"
            return "CE" + self.get_ce(seq) + "_" + ("_".join(variants) or "REF")

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
                prefix_variants = libstrnaming.call_variants(prefix, before, location="prefix")
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
                        this_variants = libstrnaming.call_variants(gap_ref, gap, location="suffix")
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
            blocks.append("%s[%i]" % (unit, (end - start) // len(unit)))
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
                suffix_variants = libstrnaming.call_variants(suffix, after, location="suffix")
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
