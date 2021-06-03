STRNaming Changelog
===================

### Version 1.1.0 (in progress)
* Fixed HPRTB allele numbering: previously, the CE allele number calculated
  for a given sequence was one higher than it should be.
* Reference sequence analysis now guarantees that all repeat units in the
  final result are actually repeated.
* Reference repeat units only found outside the reported range are now included
  in the list of preferred units when generating allele names.
* Don't consider names that include an interruption of which the sequence is equal
  to an adjacent repeat unit (e.g., CCTA[2]CCTA[1]TCTA[2]).
* New built-in reference sequence cache and reference structure plumbing.
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


### Version 1.0.0 (21 december 2020)
Initial release.
