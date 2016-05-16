# Haplofinder
#
# Copyright (c) 2001, 2016, Roslin Institute and Andy Law
#
#     This program is free software; you can redistribute it and/or modify
#     it under the terms of the GNU General Public License as published by
#     the Free Software Foundation; either version 3 of the License, or
#     (at your option) any later version.
# 
#     This program is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License for more details.
# 
#     You should have received a copy of the GNU General Public License
#     along with this program; if not, write to the Free Software
#     Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
#
#
# If you have problems or questions, then send an email to:
#   Roslin.Bioinformatics@roslin.ed.ac.uk
#
# or write to:
#
# Roslin Bioinformatics Group,
# The Roslin Institute,
# The University of Edinburgh,
# Easter Bush,
# Midlothian,
# EH25 9RG,
# Scotland, UK.
#
#
# Details
#
#
# A python script to identify which pair of known sequence haplotypes a
# given ambiguous sequence matches


The python script associated with this README file implements an incredibly
dumb, string-matching algorithm designed to dis-ambiguate possible pairs of
bovine DRB3 sequences from within sequence reads derived from sequencing PCR
products from diploid individuals. The original code was used in the paper:

Establishment of a sequence-based typing system for BoLA-DRB3 exon 2 (2003)
Miltiadou D, Law AS, Russell GC.
Tissue Antigens. 2003 Jul;62(1):55-65.


Given that there are a defined set of sequence haplotypes for the region of
interest and that each can be aligned one to the other, the predicted ambiguous
sequence reads resulting from each possible pair of haplotypes can be
calculated and then used as a reference against which to search diploid
sequences obtained from samples from individual animals. In order to do so, the
known haplotypes must first be aligned and stored as a single file containing
each haplotype reported on a common coordinate system. Novel sequence files -
also based on the same coordinate system - can then be queried against the
reference database looking for matches.

There are no fancy graph-based memory structures to behold here - simple string
matching does the job. As a result, haplotypes that contain deletions will not
be found using this method. However, at the time of writing the original code
those variants were so few and far between that the very act of identifying
that one such haplotype might be contained within a sample genotype (because no
matches to other known haplotypes were seen) was sufficient to narrow down the
search space to a tractable number of alternative options. In fact, the
deletion was usually obvious from the sequence read which descended into
consistent ambiguity once the deletion was passed so sequences that looked like
they contained anything that the program could not deal with were usually
spotted early in the analysis.

As examples, we have included an early DRB3 database and some example sequence
files generated at random from pairs of haplotypes extracted from that database.

The files, which are all in the examples directory, are as follows...

(MD5 checksum in brackets)

DRB3Ex2DNA.fas   (df224cae8e6d58c7e2dfb2854842fccd)

file-0.seq       (332d17f894f02e69cfc79288570f5e6c)
file-1.seq       (58148e4c56bcfeddff311312ef715c47)
file-2.seq       (475aa53afc862af93200c16ea2233b96)
file-3.seq       (55a36829ea69ad07bf70c5eaf4701425)
file-4.seq       (0ce35ba6e6c4ca702f69b46bf58d45fc)
file-5.seq       (e25f645b5761df28207302931276a90c)

offsets.seq      (03439c2ff6e3689b3d7bea8d47801722)


The DRB3Ex2DNA.fas file contains the aligned reference sequences. Each is
padded with dash characters so that all sequences have the same overall length
and all share a common coordinate system.

Running the haplofinder code for the first time supplying the DRB3Ex2DNA.fas
file as the alignment file will result in the generation of a file called
DRB3Ex2DNA.fas.haplos alongside the original reference file. This contains all
possible pairs of haplotypes as aligned sequences containing ambiguity codes.
These are the sequence strings against which other test sequences are aligned.

From within the project top-level directory (alongside this README.txt file),
running:

> python bin/haplofinder.py -a examples/DRB3Ex2DNA.fas examples/file-*.seq

should provide detail of the sequences matched and the pairs of haplotypes that
each of those files contain. Hopefully the output should be self-explanatory.


To demonstrate the 'offset' functionality, we have also included an edited
example file called offsets.seq. Running this file as the test sequence should
report "nothing found, even after trimming". That is because this sequence does
not directly align against the reference sequences. Rather it has had the
'missing' sequence removed from the start. This means that the sequence will
not align against the references and no matches will be found unless we tell
the code to offset the test sequence against the reference coordinates. We can
do this by using the --offset flag.

Try running:

> python bin/haplofinder.py -a examples/DRB3Ex2DNA.fas examples/offsets.seq

...to test that we are correct that nothing can be found. Then run the same
command but this time with an offset of 19 which is the number of dashes/bases
that have neen trimmed from the start of the sequence. This should allow the
sequence to align properly with the references and the underlying pair of
haplotypes can be identified.


As one further note, the sequences compared are checked to the fullest extent
that the combination of reference and test sequence allow. If the reference
sequence is shorter than the test sequence or vice-versa, each will only be
compared across the region that overlaps. For some haplotypes, this will result
in more than one possible pair of alleles being reported.




