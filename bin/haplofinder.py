#!/usr/local/bin/python

#
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
#
# Author: Andy Law
#
# November 2001
#
# Inputs:
#
# -a|--alignfile     The path to a FASTA-formatted file of unambiguous
#                    haplotype sequences.
#
# -o|--offset        The offset to use when comparing the input ambiguous
#                    sequence against the haplotype combinations
#
# <file>             The path to a FASTA-formatted file of ambiguous
#                    sequences

#
#
# Version 1.1
#
# Andy Law
# 2003-10-01
#
# Changed regex references to re objects as regex is deprecated in Python 2.2
#
# Version 1.2
#
# Andy Law
# 2016-05-13
#
# Updated contact details.
# Reformatted to conform with PEP8
#

import copy
import getopt
import os
import re
import string
import sys

MAX_CLIP_START = 20


# ----------------------------------------------------------------------------

class FastaFile:
    """
    A file of FASTA-Formatted sequences. Can be read or written.
    """

    def __init__(self, filename):
        """
        Initialise the sequence list and associated index and store the
        filename for later use

        :param filename: The file to be read from or written to
        :return:
        """

        self._sequence_array = []
        self._next_sequence_index = 0
        self._filename = filename

    # ----------------
    # In the read_file routine, we open the file (which we assume to be an
    # aligned FASTA-format file and read in each sequence in turn. We check
    # the lengths of the sequences and also keep track of the first and last
    # fully-identified (i.e. not leading or trailing '-') base

    def read_file(self):

        # Try to open the file. If we can open it, read the whole thing in in
        # one go. Split the file on \r then rejoin using \n, then split on
        # \n\n and rejoin using \n so that we get all lines delimited with a
        # single newline character

        input_file = open(self._filename, 'r')
        if input_file:
            file_contents = input_file.read()
            string.joinfields(string.splitfields(file_contents, "\r"), "\n")
            string.joinfields(string.splitfields(file_contents, "\n\n"), "\n")

            # Split the file on '>' and throw away the first record (the blank
            # space before the first '>'). We now have an array of sequence
            # records
            records = (string.splitfields(file_contents, ">"))[1:]

            # Create a regular expression object to look for quoted names in
            # the file
            my_regex = re.compile("^\s*\"([^\"]+)\"")

            # For each sequence in turn, split the record into lines (on '\n'),
            # and split the first line, keeping the first entry as the title.
            # Then join all the other lines back up as the sequence element.
            # Create a new sequence object using these details and stick it in
            # the list of known sequences

            for this_record in records:
                lines = string.splitfields(this_record, "\n")
                title = lines[0]
                my_match_object = my_regex.search(title)
                if my_match_object is None:
                    title = (string.split(lines[0]))[0]
                else:
                    title = my_match_object.group(1)
                sequence = string.joinfields(lines[1:], "")

                new_seq = Sequence(title, sequence)

                self.add_sequence(new_seq)

            # Close the file afterwards
            input_file.close()

        else:
            print "Cannot read the file '" + self._filename + "'"
            sys.exit(1)

            # ----------------
            # Return the next unread sequence from the list read in by the
            # read_file routine

    def get_next_sequence(self):

        if self._next_sequence_index >= len(self._sequence_array):
            return None
        else:
            next_seq = self._sequence_array[self._next_sequence_index]
            self._next_sequence_index += 1
            return next_seq

            # ----------------
            # Add a new sequence to the list of sequences. Sequence is assumed
            # to be a sequence object

    def add_sequence(self, sequence):
        self._sequence_array.append(sequence)

    # ----------------
    # Write the contents of the sequence array out to the file

    def write_file(self):
        output_file = open(self._filename, 'w')
        if output_file:
            for each_sequence in self._sequence_array:
                output_file.write(">\"" + each_sequence.title + "\"\n")
                output_file.write(each_sequence.sequence + "\n")
            output_file.close()


# ----------------------------------------------------------------------------

class HaplotypeDetector:
    """
    A HaplotypeDetector manages the reference haplotypes file and search it
    for matches against sequences supplied
    """

    def __init__(self, filename, offset):
        """
        Construct a new HaplotypeDetector. First check for a pre-computed
        haplotypes file that has a date stamp later than the plain reference
        haplotypes file. If one exists then use that file instead, otherwise
        compute the pairs from the original list and save them for future use.

        :param filename: The name of the file that contains the reference
        haplotypes
        :param offset: The offset to apply when comparing the test sequences
        against the reference
        :return:
        """

        # Initialise the storage structures and store the filename

        self._haplotype_list = []
        self._haplotype_dict = {}

        self._combination_list = []
        self._filename = filename

        # Look for a pre-computed file of haplotypes (filename + ".haplos")

        using_haplo_file = 0

        if os.path.exists(self._filename + ".haplos"):
            real_file_stat = os.stat(self._filename)[8]
            haplotype_file_stat = os.stat(self._filename + ".haplos")[8]
            if haplotype_file_stat > real_file_stat:
                self._filename += ".haplos"
                print "Reading pre-computed haplotype combinations from '" +\
                      self._filename + "'"
                using_haplo_file = 1

                # Read the fasta file and get each sequence in turn. Store them
                # in the Haplotype list and dictionary, having stripped off any
                # leading stuff if required

        fasta_file = FastaFile(self._filename)
        if fasta_file:
            fasta_file.read_file()
            while True:
                new_seq = fasta_file.get_next_sequence()
                if new_seq is None:
                    break
                if new_seq.title in self._haplotype_dict:
                    print "Sequence '" + new_seq.title + "' already exists"
                    sys.exit(1)

                if offset > 0:
                    new_seq = new_seq[offset:]
                self._haplotype_dict[new_seq.title] = new_seq
                self._haplotype_list.append(new_seq)

        else:
            print "Cannot read the file '" + filename + "'"
            sys.exit(1)

        if not using_haplo_file:
            print "There are " + str(
                len(self._haplotype_list)) + " original sequences"
            self.build_combinations()
        else:
            self._combination_list = self._haplotype_list

        print "Length of a sequence is " + str(
            len(self._combination_list[0].sequence))
        print "There are " + str(len(self._combination_list)) + " combinations"

    # ------------------------------------------------------------------------
    #

    def build_combinations(self):
        """
        Build a set of sequence combinations from the original set of sequences
        using just the portion of sequence specified by the start and length
        parameters

        :return:
        """
        print "Building combinations of known haplotypes"
        self._combination_list = self._haplotype_list
        original_count = len(self._combination_list)

        # Once we have all the sequences, add them to each other to give the
        # possible combinations. Add them to the list as well.

        for i in range(0, original_count - 1):
            for j in range(i + 1, original_count):
                new_seq = self._combination_list[i] + self._combination_list[j]
                self._combination_list.append(new_seq)
                sys.stdout.write(".")
                sys.stdout.flush()
        print
        print "We have " + str(
            len(self._combination_list)) + " combinations to compare"
        fasta_file = FastaFile(self._filename + ".haplos")

        for each_seq in self._combination_list:
            fasta_file.add_sequence(each_seq)
        fasta_file.write_file()

    # ------------------------------------------------------------------------
    #

    def find_match(self, sequence):

        print
        print "Analysing '" + sequence.title + "'"
        print "Start of sequence is  '" + sequence.sequence[:20] + "'"
        print "Start of reference is '" + self._haplotype_list[0].sequence[
                                          :20] + "'"

        got_match = 0
        clip_start = 0
        clipped = ''
        while (got_match == 0) and (clip_start <= MAX_CLIP_START):
            if clipped != '':
                print clipped
            for comb_seq in self._combination_list:
                if sequence[clip_start:] == comb_seq[clip_start:]:
                    got_match = 1
                    print "Seems to match '" + comb_seq.title + "' " + clipped
            clip_start += 1
            clipped = "(clipping %d bases from the start)" % clip_start

        if not got_match:
            print "Nothing matched! (even after clipping %d bases from the " \
                  "start" % MAX_CLIP_START


# ----------------------------------------------------------------------------

class IUPAC:
    def __init__(self):

        self.list_by_symbol = {}
        self.list_by_score = {}
        self._next_bit = 1

        self.seed('A')
        self.seed('C')
        self.seed('G')
        self.seed('T')

        self.add('R', self.get_bit_score(['A', 'G']))
        self.add('Y', self.get_bit_score(['T', 'C']))

        self.add('M', self.get_bit_score(['A', 'C']))
        self.add('K', self.get_bit_score(['G', 'T']))

        self.add('S', self.get_bit_score(['G', 'C']))
        self.add('W', self.get_bit_score(['A', 'T']))

        self.add('B', self.get_bit_score(['C', 'G', 'T']))
        self.add('D', self.get_bit_score(['A', 'G', 'T']))
        self.add('H', self.get_bit_score(['A', 'C', 'T']))
        self.add('V', self.get_bit_score(['A', 'C', 'G']))

        self.add('N', self.get_bit_score(['A', 'C', 'G', 'T']))

    def seed(self, label):
        score = self._next_bit
        self._next_bit *= 2
        self.add(label, score)

    def add(self, label, score):
        self.list_by_symbol[label] = score
        self.list_by_score[score] = label

    def get_bit_score(self, base_list):
        score = 0
        for base in base_list:
            score = score | self.list_by_symbol[base]
        return score

    def get_merged_code(self, base_list):
        if '-' in base_list:
            return '-'
        else:
            return self.list_by_score[self.get_bit_score(base_list)]


# ----------------------------------------------------------------------------

class Sequence:
    # Class-based IUPAC object. This is just used as a translator so we don't
    # need separate instances for each Sequence object.
    classIUPAC = IUPAC()

    def __init__(self, title, sequence):
        self.title = title
        self.sequence = string.upper(sequence)
        self.parents = []

    def __add__(self, other):
        if len(self.sequence) != len(other.sequence):
            raise Exception("Can't add 2 sequences of differing length")

        newseq = ''
        for i in range(0, len(self.sequence)):
            newseq = newseq + self.classIUPAC.get_merged_code(
                [self.sequence[i:i + 1], other.sequence[i:i + 1]])

        sequence_object = Sequence(self.title + " + " + other.title, newseq)
        sequence_object.add_parents([self, other])
        return sequence_object

    def __getslice__(self, i, j):
        new = copy.copy(self)
        new.sequence = new.sequence[i:j]
        return new

    def __cmp__(self, other):
        if other is None:
            return 1
        cmplen = min(self.length(), other.length())
        seq1 = self.sequence[:cmplen]
        seq2 = other.sequence[:cmplen]
        while (seq1[:1] == '-') or (seq2[:1] == '-'):
            seq1 = seq1[1:]
            seq2 = seq2[1:]
        while (seq1[-1:] == '-') or (seq2[-1:] == '-'):
            seq1 = seq1[:-1]
            seq2 = seq2[:-1]
        if seq1 == seq2:
            return 0
        else:
            return 1

    def add_parents(self, parents):
        for parent in parents:
            self.parents.append(parent)

    def length(self):
        return len(self.sequence)


class HaploFinder:
    def __init__(self):
        # Get the command line arguments
        #
        try:
            optlist, self.remains = getopt.getopt(sys.argv[1:], 'a:o:cw',
                                                  ['alignfile=', 'offset='])
        except getopt.GetoptError as err:
            print "There was a problem with the command line arguments - " \
                  "exiting"
            print str(err)
            sys.exit(1)

        #
        # Decode them so short arguments are converted to long arguments
        #
        short_to_long = {'a': "alignfile", 'o': "offset"}
        self.options = {}

        for arg in optlist:
            flag = arg[0]
            value = arg[1]

            while flag[0] == '-':
                flag = flag[1:]

            if flag in short_to_long:
                flag = short_to_long[flag]

            if flag in self.options:
                err = "Too many instances of the '%s' argument" % flag
                raise Exception(err)

            self.options[flag] = value

        if 'offset' not in self.options:
            self.options["offset"] = 0

        self.options["offset"] = int(self.options["offset"])

        exit_needed = False
        if 'c' in self.options:
            HaploFinder.copying()
            exit_needed = True
        if 'w' in self.options:
            HaploFinder.warranty()
            exit_needed = True
        if exit_needed is True:
            sys.exit(0)

        # flag the conditions and report any other requested information
        HaploFinder.start_details()
        if ('alignfile' not in self.options) or (len(self.remains) == 0):
            HaploFinder.help()

    @staticmethod
    def help():
        print ("""\
 Usage: %s -a|--alignfile <alignmentfile>
             [-o|-offset <offset>] file1 [file2 ...]

        Where <alignmentfile> is a padded alignment in FASTA file
                    format

            <offset> is the offset of the sequences in the files to be
                     analysed, compared against the alignment. Positive
                     offsets mean the sequence starts within the alignment,
                     negative offsets (enclosed in quotes) mean the sequence
                     starts before the sequence in the alignment file

               file1 [file2 ...] are the names of FASTA formatted
                     sequence files to be analysed
""" % (sys.argv[0]))
        sys.exit(0)

    def run(self):

        # Build the haplotype detector by reading in the file

        hd = HaplotypeDetector(self.options['alignfile'],
                               self.options['offset'])

        for input_file in self.remains:
            fasta_file = FastaFile(input_file)
            fasta_file.read_file()

            while True:
                seq = fasta_file.get_next_sequence()
                if seq is None:
                    break
                if self.options["offset"] < 0:
                    seq = seq[(self.options["offset"] * -1):]
                hd.find_match(seq)

    @staticmethod
    def start_details():
        print ("""\

Haplofinder version 1.2, Copyright (C) 2001,2016 Roslin Institute, Andy Law

Haplofinder comes with ABSOLUTELY NO WARRANTY; for details type:
  '%s -w'.

This is free software, and you are welcome to redistribute it under certain
conditions; for details type:
  '%s -c'

""" % (sys.argv[0], sys.argv[0]))

    @staticmethod
    def copying():
        print ("""\

  4. Conveying Verbatim Copies.

  You may convey verbatim copies of the Program's source code as you
receive it, in any medium, provided that you conspicuously and
appropriately publish on each copy an appropriate copyright notice;
keep intact all notices stating that this License and any
non-permissive terms added in accord with section 7 apply to the code;
keep intact all notices of the absence of any warranty; and give all
recipients a copy of this License along with the Program.

  You may charge any price or no price for each copy that you convey,
and you may offer support or warranty protection for a fee.

  5. Conveying Modified Source Versions.

  You may convey a work based on the Program, or the modifications to
produce it from the Program, in the form of source code under the
terms of section 4, provided that you also meet all of these conditions:

    a) The work must carry prominent notices stating that you modified
    it, and giving a relevant date.

    b) The work must carry prominent notices stating that it is
    released under this License and any conditions added under section
    7.  This requirement modifies the requirement in section 4 to
    "keep intact all notices".

    c) You must license the entire work, as a whole, under this
    License to anyone who comes into possession of a copy.  This
    License will therefore apply, along with any applicable section 7
    additional terms, to the whole of the work, and all its parts,
    regardless of how they are packaged.  This License gives no
    permission to license the work in any other way, but it does not
    invalidate such permission if you have separately received it.

    d) If the work has interactive user interfaces, each must display
    Appropriate Legal Notices; however, if the Program has interactive
    interfaces that do not display Appropriate Legal Notices, your
    work need not make them do so.

  A compilation of a covered work with other separate and independent
works, which are not by their nature extensions of the covered work,
and which are not combined with it such as to form a larger program,
in or on a volume of a storage or distribution medium, is called an
"aggregate" if the compilation and its resulting copyright are not
used to limit the access or legal rights of the compilation's users
beyond what the individual works permit.  Inclusion of a covered work
in an aggregate does not cause this License to apply to the other
parts of the aggregate.

  6. Conveying Non-Source Forms.

  You may convey a covered work in object code form under the terms
of sections 4 and 5, provided that you also convey the
machine-readable Corresponding Source under the terms of this License,
in one of these ways:

    a) Convey the object code in, or embodied in, a physical product
    (including a physical distribution medium), accompanied by the
    Corresponding Source fixed on a durable physical medium
    customarily used for software interchange.

    b) Convey the object code in, or embodied in, a physical product
    (including a physical distribution medium), accompanied by a
    written offer, valid for at least three years and valid for as
    long as you offer spare parts or customer support for that product
    model, to give anyone who possesses the object code either (1) a
    copy of the Corresponding Source for all the software in the
    product that is covered by this License, on a durable physical
    medium customarily used for software interchange, for a price no
    more than your reasonable cost of physically performing this
    conveying of source, or (2) access to copy the
    Corresponding Source from a network server at no charge.

    c) Convey individual copies of the object code with a copy of the
    written offer to provide the Corresponding Source.  This
    alternative is allowed only occasionally and noncommercially, and
    only if you received the object code with such an offer, in accord
    with subsection 6b.

    d) Convey the object code by offering access from a designated
    place (gratis or for a charge), and offer equivalent access to the
    Corresponding Source in the same way through the same place at no
    further charge.  You need not require recipients to copy the
    Corresponding Source along with the object code.  If the place to
    copy the object code is a network server, the Corresponding Source
    may be on a different server (operated by you or a third party)
    that supports equivalent copying facilities, provided you maintain
    clear directions next to the object code saying where to find the
    Corresponding Source.  Regardless of what server hosts the
    Corresponding Source, you remain obligated to ensure that it is
    available for as long as needed to satisfy these requirements.

    e) Convey the object code using peer-to-peer transmission, provided
    you inform other peers where the object code and Corresponding
    Source of the work are being offered to the general public at no
    charge under subsection 6d.

  A separable portion of the object code, whose source code is excluded
from the Corresponding Source as a System Library, need not be
included in conveying the object code work.

  A "User Product" is either (1) a "consumer product", which means any
tangible personal property which is normally used for personal, family,
or household purposes, or (2) anything designed or sold for incorporation
into a dwelling.  In determining whether a product is a consumer product,
doubtful cases shall be resolved in favor of coverage.  For a particular
product received by a particular user, "normally used" refers to a
typical or common use of that class of product, regardless of the status
of the particular user or of the way in which the particular user
actually uses, or expects or is expected to use, the product.  A product
is a consumer product regardless of whether the product has substantial
commercial, industrial or non-consumer uses, unless such uses represent
the only significant mode of use of the product.

  "Installation Information" for a User Product means any methods,
procedures, authorization keys, or other information required to install
and execute modified versions of a covered work in that User Product from
a modified version of its Corresponding Source.  The information must
suffice to ensure that the continued functioning of the modified object
code is in no case prevented or interfered with solely because
modification has been made.

  If you convey an object code work under this section in, or with, or
specifically for use in, a User Product, and the conveying occurs as
part of a transaction in which the right of possession and use of the
User Product is transferred to the recipient in perpetuity or for a
fixed term (regardless of how the transaction is characterized), the
Corresponding Source conveyed under this section must be accompanied
by the Installation Information.  But this requirement does not apply
if neither you nor any third party retains the ability to install
modified object code on the User Product (for example, the work has
been installed in ROM).

  The requirement to provide Installation Information does not include a
requirement to continue to provide support service, warranty, or updates
for a work that has been modified or installed by the recipient, or for
the User Product in which it has been modified or installed.  Access to a
network may be denied when the modification itself materially and
adversely affects the operation of the network or violates the rules and
protocols for communication across the network.

  Corresponding Source conveyed, and Installation Information provided,
in accord with this section must be in a format that is publicly
documented (and with an implementation available to the public in
source code form), and must require no special password or key for
unpacking, reading or copying.

The above is an excerpt of the full GPL license which you should have
received with this program. Please refer to the full document for
more details.""")

    @staticmethod
    def warranty():
        print ("""\

  15. Disclaimer of Warranty.

  THERE IS NO WARRANTY FOR THE PROGRAM, TO THE EXTENT PERMITTED BY
APPLICABLE LAW.  EXCEPT WHEN OTHERWISE STATED IN WRITING THE COPYRIGHT
HOLDERS AND/OR OTHER PARTIES PROVIDE THE PROGRAM "AS IS" WITHOUT WARRANTY
OF ANY KIND, EITHER EXPRESSED OR IMPLIED, INCLUDING, BUT NOT LIMITED TO,
THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
PURPOSE.  THE ENTIRE RISK AS TO THE QUALITY AND PERFORMANCE OF THE PROGRAM
IS WITH YOU.  SHOULD THE PROGRAM PROVE DEFECTIVE, YOU ASSUME THE COST OF
ALL NECESSARY SERVICING, REPAIR OR CORRECTION.

  16. Limitation of Liability.

  IN NO EVENT UNLESS REQUIRED BY APPLICABLE LAW OR AGREED TO IN WRITING
WILL ANY COPYRIGHT HOLDER, OR ANY OTHER PARTY WHO MODIFIES AND/OR CONVEYS
THE PROGRAM AS PERMITTED ABOVE, BE LIABLE TO YOU FOR DAMAGES, INCLUDING ANY
GENERAL, SPECIAL, INCIDENTAL OR CONSEQUENTIAL DAMAGES ARISING OUT OF THE
USE OR INABILITY TO USE THE PROGRAM (INCLUDING BUT NOT LIMITED TO LOSS OF
DATA OR DATA BEING RENDERED INACCURATE OR LOSSES SUSTAINED BY YOU OR THIRD
PARTIES OR A FAILURE OF THE PROGRAM TO OPERATE WITH ANY OTHER PROGRAMS),
EVEN IF SUCH HOLDER OR OTHER PARTY HAS BEEN ADVISED OF THE POSSIBILITY OF
SUCH DAMAGES.

The above is an excerpt of the full GPL license which you should have
received with this program. Please refer to the full document for
more details""")


if __name__ == '__main__':
    HaploFinder().run()
