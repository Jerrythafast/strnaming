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
"""
1	230769605	230769681	-8
4	154587721	154587804	-12
5	150076325	150076397	-20
7	84160205	84160277	-20
9	77073818	77073874	-8
12	5983958	5984046	-20
13	82148025	82148089	-20
18	63281668	63281757	-17
19	29926234	29926286	+4
21	19181973	19182102	-13
21	43636196	43636290	-29
X	9338319	9338399	+12
X	67757315	67757402	-31
X	134284910	134284930	-20
X	134284959	134285039	-4
X	134481506	134481562	-4
X	150542522	150542590	-8
Y	9684380	9684490	-54
Y	11982075	11982159	-40
Y	12500423	12500495	-24
Y	12500542	12500614	-4
Y	13640728	13640836	-18
Y	15163053	15163169	-20
Y	16631624	16631652	-28
Y	16631673	16631757	-36
Y	18639703	18639875	-132
Y	18888808	18888852	+4
Y	22218923	22219079	-46
Y	23785334	23785514	-43
"""
LENGTH_ADJUSTMENTS = {}

def get_adjustment(chromosome, start, end, multiple_structures):
    """
    Gets the adjustment for obtaining correct CE allele numbers at the
    given range.

    The end position is exclusive. Set multiple_structures to True if
    the input range covers multiple STR structures. If set to False, no
    adjustment will be made if the only included STR structure is one
    that was historically ignored in CE allele numbering.
    """
    adjustment = 0
    if chromosome in LENGTH_ADJUSTMENTS:
        for range_start, range_end, range_adjustment in LENGTH_ADJUSTMENTS[chromosome]:
            if start < range_end and end > range_start:
                if multiple_structures or range_adjustment != range_start - range_end:
                    adjustment += range_adjustment
    return adjustment

def _initialize():
    for line in __doc__.split("\n"):
        if not line:
            continue
        chromosome, start, end, adjustment = line.strip().split()
        if chromosome not in LENGTH_ADJUSTMENTS:
            LENGTH_ADJUSTMENTS[chromosome] = []
        LENGTH_ADJUSTMENTS[chromosome].append((int(start), int(end), int(adjustment)))
_initialize()
del _initialize
