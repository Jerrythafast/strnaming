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

    strnaming --ranges uas-frr inputfile.txt outputfile.txt

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


Release Notes
-------------
### Version 1.0.0 (21 december 2020)
Initial release of STRNaming.
