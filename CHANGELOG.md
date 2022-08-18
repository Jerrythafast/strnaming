STRNaming Changelog
===================

### Version 1.1.3 (18 August 2022)
* Fixed an issue that caused STRNaming to sometimes favour a longer name with
  the same score.

### Version 1.1.2 (10 May 2022)
* Updated CE allele numbering of DYS612 (+6 units).
* Added reference structure for SE33.
* Updated hardcoded reference length adjustment table to suppress second
  structure 5' of DYS522.
* Added double-click-to-toggle-text-alignment feature to HTML output.
* Improved reference sequence analysis speed.
* Fixed mtDNA reference sequence download URL.

### Version 1.1.1 (19 July 2021)
* Fixed an issue with CE allele numbering that occurred for reporting ranges
  that started or ended halfway into a structure with a hardcoded reference
  length adjustment.
* Updated table of hardcoded reference length adjustments to include more loci.

### Version 1.1.0 (15 July 2021)
* Updated scoring criteria to minimize unintended side-effects of other changes.
* Fixed HPRTB allele numbering: previously, the CE allele number calculated
  for a given sequence was one higher than it should be.
* Allele names are now permitted to contain repeats of a unit that exceeds the
  dominant unit length of a locus.
* Short repeat stretches that only partially overlap with a significant repeat
  of a longer unit are no longer discarded.
* Fixed bug that disallowed making interruptions which could be filled exactly
  with an 'orphan' repeat, thereby forcing the use of a compatible 'anchor'.
* Reference sequence analysis now guarantees that all repeat units in the
  final result are actually repeated.
* Reference repeat units only found outside the reported range are now included
  in the list of preferred units when generating allele names.
* The minimum repeat stretch length has been reduced from 4 nt to 3 nt.
* Stretches of reference repeat units may now be trimmed down to 3 nt (was 5) to
  combine them with other overlapping stretches; stretches of non-reference
  repeat units now have a higher minimum of 8 nt (up from 5).
* Mononucleotide stretches of 5nt are no longer considered a valid repeat when
  overlapped completely by a hexanucleotide repeat unit.
* When multiple structures with differing numbers of repeat stretches result in
  the same score, always select the one with the fewest repeat stretches.
* Don't consider names that include an interruption of which the sequence is equal
  to an adjacent repeat unit (e.g., CCTA[2]CCTA[1]TCTA[2]).
* New built-in reference sequence cache and bundled reference structure data for
  nearly the entire human genome.
* New mandatory ACTION command-line argument to choose between allele naming
  and reference sequence downloading.
* New options for generating HTML-formatted output with repeat unit coloring.
* The ReferenceSequenceStore.add_refseq() function can now handle overlapping
  ranges of reference sequence.
* Improved support for naming non-STR ranges, including mtDNA.
* Added support for calculating the sequence for a given name.
* Greatly improved performance of reverse_complement().
* New ComplexReportedRange class that supports combining parts of reference
  sequence (but does not support STR structures as of now). Useful for reporting
  mtDNA and mRNA data.
* The OutOfTimeException is now handled gracefully (by writing a message to the
  standard error stream and reporting no repeat stretches).
* Reference sequence analysis will now raise a ComplexityException if the
  provided reference sequence results in generating over 5 million scaffolds.


### Version 1.0.0 (21 December 2020)
Initial release.
