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
"""
1	163590026	163590094	-8	#	D1S1677	15
1	230769605	230769681	-8	#	D1S1656	17
2	168788875	168788937	-18	#	D2S1776	11
4	154587721	154587821	-12	#	FGA	22
5	59403131	59403191	+8	#	D5S2800	17
5	150076325	150076397	-20	#	CSF1PO	13
6	88277144	88277356	-110	#	SE33	25.2
6	112557951	112558019	-4	#	D6S474	16
7	84160205	84160277	-20	#	D7S820	13
9	77073818	77073874	-8	#	D9S1122	12
12	5983958	5984046	-20	#	vWA	17
13	82148025	82148089	-20	#	D13S317	11
18	63281668	63281757	-17	#	D18S51	18
19	29926234	29926286	+4	#	D19S433	14
21	19181973	19182102	-13	#	D21S11	29
21	43636196	43636290	-29	#	PentaD	13
X	67496061	67496165	-20	#	DXS10079	21
X	67757315	67757402	-31	#	DXS10074	14
X	134284910	134284930	-20	#	(before DXS10103)	0
X	134284969	134285041	+4	#	DXS10103	19
X	134481506	134481562	-4	#	HPRTB	14
X	150542522	150542590	-8	#	DXS7423	15
Y	7547447	7547470	-23	#	(before DYS522)	0
Y	7999821	7999903	-18	#	DYS458	16
Y	8558264	8558284	-20	#	(before DYS481)	0
Y	8781946	8782056	-26	#	DYS627	21
Y	9684376	9684490	-54	#	DYS19	15
Y	11982075	11982159	-40	#	DYS391	11
Y	12500423	12500495	-24	#	DYS389I	12
Y	12500543	12500623	-12	#	DYS389II	17
Y	15163053	15163169	-20	#	DYS390	24
Y	16631624	16631652	-28	#	(before Y-GATA-H4)	0
Y	16631673	16631757	-36	#	Y-GATA-H4	12
Y	18639700	18639876	-132	#	DYS385b	11
Y	18680513	18680701	-132	#	DYS385a	14
Y	18888808	18888852	+4	#	DYS461	12
Y	22218923	22219083	-46	#	DYS448	19
Y	23785331	23785514	-43	#	DYF387S1 fragment 1	35
Y	25884568	25884755	-43	#	DYF387S1 fragment 2	36
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
        chromosome, start, end, adjustment = line.split("#", maxsplit=1)[0].strip().split()
        if chromosome not in LENGTH_ADJUSTMENTS:
            LENGTH_ADJUSTMENTS[chromosome] = []
        LENGTH_ADJUSTMENTS[chromosome].append((int(start), int(end), int(adjustment)))
_initialize()
del _initialize
