STRNaming
=========
STRNaming is an algorithm for generating simple, informative names for Short
Tandem Repeat (STR) sequences, such as those used in the field of forensic
genetics, in a standardised and automated manner.


Requirements
------------
STRNaming requires Python version 3.5 or later.


Installation
------------
The recommended way to install STRNaming is by using the `pip` package
installer. If you have `pip` installed, you can easily install STRNaming by
running the following command:

    pip install strnaming

Alternatively, STRNaming can be installed by running:

    python setup.py install


Usage
-----
This initial version of STRNaming allows generating allele names for sequence
data using the ranges and sequence orientation of the "Flanking Region Report"
of the Universal Analysis Software for the ForenSeq DNA Signature Prep Kit
(Verogen).

### Command-line interface
The command-line help can be accessed by running `strnaming --help`. In short,
an STRNaming command looks like this:

    strnaming name-sequences --ranges uas-frr inputfile.txt outputfile.txt

The input file should have a marker name and a sequence on each line, separated
by whitespace (i.e., tabs or spaces).

If no output file is given, the output is written to `stdout`, which normally
shows up in your command line window. If no input file is given either,
STRNaming will read input from `stdin`, allowing you to type the input one line
at a time.

### Programming interface
It is **not recommended** to `import` and use parts of this version of
STRNaming directly from other Python code, because the internal API is not
stable yet. Instead, use the `subprocess` module if you want to use STRNaming
in your Python project at this time. As an added benefit, it will run in a
concurrent process, meaning your code does not (necessarily) have to wait for
STRNaming to finish.

To use STRNaming in other software projects, regardless of the programming
language, it can be run as a separate subprocess. Write a marker name, a
whitespace character, the DNA sequence, and a newline character (`\n`) to its
standard input stream (`stdin`), and STRNaming will write the same marker name,
a tab character, the allele name and a newline character to its standard output
stream (`stdout`). Any errors are reported on the standard error stream
(`stderr`) and will cause the STRNaming process to terminate. By specifying the
`--unbuffered` command-line switch, STRNaming will immediately flush its output
stream after every line of output.

### Offline use
STRNaming will automatically download and cache portions of reference sequence
from the Ensembl REST API (http://rest.ensembl.org). If you are running
STRNaming on a system without internet access, and you need a piece of
reference sequence that was not bundled with the STRNaming package, a message
will be displayed to manually store the reference sequence in a specific
location. To this end, run the following command (on a system with internet
access) to download the sequence:

    strnaming refseq-cache chr2:1489653..1489689

Upon success, the location of the downloaded cache files will be displayed.
These are the files to be copied to the offline system for STRNaming to work.


Release Notes
-------------
### Version 1.2.0 (TBD)
Naming of some loci has been updated as a result of bug fixes and improvements
to the algorithm. Most notably, reference sequence analysis has been redesigned
in such a way that it is no longer affected by the range of reference sequence
analysed at once.
* Updated CE allele numbering of D6S474 (-1 unit).

### Version 1.1.4 (7 February 2023)
* Repeat stretches that fall completely in the prefix or suffix are now ignored.
* Sequences that follow the same repeat pattern as the reference sequence are
  now named much more quickly while consistently using the same structure.
* Added capability to load reference structures from many locations on one
  chromosome in a single pass.

### Version 1.1.3 (18 August 2022)
* Fixed an issue that caused STRNaming to sometimes favour a longer name with
  the same score.

### Version 1.1.2 (10 May 2022)
* Updated CE allele numbering of DYS612 (+6 units).
* Added reference structure for SE33.
* Updated hardcoded reference length adjustment table to suppress second
  structure 5' of DYS522.
* Added double-click-to-toggle-text-alignment feature to HTML output.
* Fixed mtDNA reference sequence download URL.

### Version 1.1.1 (19 July 2021)
* Fixed an issue with CE allele numbering that occurred for reporting ranges
  that started or ended halfway into a structure with a hardcoded reference
  length adjustment.
* Updated table of hardcoded reference length adjustments to include more loci.

### Version 1.1.0 (15 July 2021)
Naming of some loci has been updated as a result of bug fixes and improvements
to the algorithm. Scoring criteria have been updated to minimize unintended
side-effects of these changes.
* Fixed a major issue with HPRTB allele numbering: previously, the CE allele
  number calculated for a given sequence was one higher than it should be.
* Allele names are now permitted to contain repeats of a unit that exceeds the
  dominant unit length of a locus. This change greatly improves naming of some
  complex Y-STRs.
* Short repeat stretches that only partially overlap with a significant repeat
  of a longer unit are no longer discarded. This change may introduce short
  repeats adjacent to longer repeats of a longer unit, which were previously
  'missed' by STRNaming.
* Fixed bug that disallowed making interruptions which could be filled exactly
  with an 'orphan' repeat, thereby forcing the use of a compatible 'anchor'.
* Reference sequence analysis now guarantees that all repeat units in the
  final result are actually repeated.
* Reference repeat units only found outside the reported range are now included
  in the list of preferred units when generating allele names. This change
  improves naming stability when a significant part of the reference STR
  structure lies outside the reported range.
* STRNaming will no longer consider names that include an interruption of which
  the sequence is equal to an adjacent repeat unit (e.g., CCTA[2]CCTA[1]TCTA[2]).

New features:
* The built-in reference sequence cache was introduced, along with the new
  mandatory ACTION command-line argument.
* Colored output in HTML format is now available by using the --html
  command-line argument.
* Reference sequence analysis results of almost the entire human genome have
  been embedded into the package.

### Version 1.0.0 (21 December 2020)
Initial release of STRNaming.
