#!/usr/local/bin/python

#
# Haplofinder
#
# Copyright (c) 2001, 2016, Roslin Institute and Andy Law
#
#     This program is free software; you can redistribute it and/or modify
#     it under the terms of the GNU General Public License as published by
#     the Free Software Foundation; either version 2 of the License, or
#     (at your option) any later version.
# 
#     This program is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License for more details.
# 
#     You should have received a copy of the GNU General Public License
#     along with this program; if not, write to the Free Software
#     Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
#
#
# If you have problems or questions, then send an email to 
# Roslin.Bioinformatics@roslin.ed.ac.uk, or write to Roslin Bioinformatics Group,
# Roslin Institute, Roslin, Midlothian, EH25 9RG, Scotland, UK.
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
#

import copy
import getopt
import os
import re
import string
import sys

MAX_CLIP_START = 20


# --------------------------------------------------------------------------------
#
# FASTA_File
#
# A file of FASTA-Formatted sequences. Can be read or written.
#
# Constructor just creates an object with a set of empty structures

class FastaFile:
    # ----------------
    # In the constructor, we initialise the sequenceArray and indices and store the
    # filename

    def __init__(self, filename):

        self._sequence_array = []
        self._next_sequence_index = 0
        self._filename = filename

    # ----------------
    # In the read_file routine, we open the file (which we assume to be an aligned FASTA-
    # format file and read in each sequence in turn. We check the lengths of the
    # sequences and also keep track of the first and last fully-identified (i.e.
    # not leading or trailing '-') base

    def read_file(self):

        # Try to open the file. If we can open it, read the whole thing in in one
        # go. Split the file on \r then rejoin using \n , then split on \n\n and rejoin
        # using \n so that we get all lines delimited with a single newline character

        input_file = open(self._filename, 'r')
        if input_file:
            file_contents = input_file.read()
            string.joinfields(string.splitfields(file_contents, "\r"), "\n")
            string.joinfields(string.splitfields(file_contents, "\n\n"), "\n")

            # Split the file on '>' and throw away the first record (the blank space before
            # the first '>'). We now have an array of sequence records
            records = (string.splitfields(file_contents, ">"))[1:]

            # Create a regular expression object to look for quoted names in the file
            my_regex = re.compile("^\s*\"([^\"]+)\"")

            # For each sequence in turn, split the record into lines (on '\n'), and split
            # the first line, keeping the first entry as the title. Then join all the other
            # lines back up as the sequence element. Create a new sequence object using these
            # details and stick it in the list of known sequences

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
            print "Cannot read the file '" + self.filename + "'"
            sys.exit(1)

        # ----------------
        #
        # Return the next unread sequence from the list read in by the read_file routine
        #

    def get_next_sequence(self):

        if self._next_sequence_index >= len(self._sequence_array):
            return None
        else:
            next_seq = self._sequence_array[self._next_sequence_index]
            self._next_sequence_index += 1
            return next_seq

        # ----------------
        #
        # Add a new sequence to the list of sequences. Sequence is assumed to be a sequence
        # object
        #

    def add_sequence(self, sequence):
        self._sequence_array.append(sequence)

    # ----------------
    #
    # Write the contents of the sequence array out to the file
    #
    def write_file(self):
        output_file = open(self.filename, 'w')
        if output_file:
            for each_sequence in self._sequence_array:
                output_file.write(">\"" + each_sequence.title + "\"\n")
                output_file.write(each_sequence.sequence + "\n")
            output_file.close()


# --------------------------------------------------------------------------------

class HaplotypeDetector:
    #
    # Constructor
    #

    def __init__(self, filename, offset):

        # Initialise the storage structures and store the filename

        self._haplotype_list = []
        self._haplotype_dict = {}

        self._combination_list = []
        self._filename = filename

        # Look for a pre-computed file of haplotypes (filename + ".haplos")
        # If one exists, and has been modified later then the input filename then
        # use that one

        using_haplo_file = 0

        if os.path.exists(self._filename + ".haplos"):
            real_file_stat = os.stat(self._filename)[8]
            haplotype_file_stat = os.stat(self._filename + ".haplos")[8]
            if haplotype_file_stat > real_file_stat:
                self.filename += ".haplos"
                print "Reading pre-computed haplotype combinations from '" + self._filename + "'"
                using_haplo_file = 1

            # Read the fasta file and get each sequence in turn. Store them in the Haplotype Array
            # and hash, having stripped off any leading stuff if required

        fasta_file = FASTA_File(self._filename)
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
            print "There are " + str(len(self._haplotype_list)) + " original sequences"
            self.build_combinations()
        else:
            self._combination_list = self._haplotype_list

        print "Length of a sequence is " + str(len(self._combination_list[0].sequence))
        print "There are " + str(len(self._combination_list)) + " combinations"

    # --------------------------------------------------------------------------------
    #
    # Build a set of sequence combinations from the original set of sequences
    # using just the portion of sequence specified by the start and length parameters

    def build_combinations(self):
        print "Building combinations of known haplotypes"
        self._combination_list = self._haplotype_list
        original_count = len(self._combination_list)

        # Once we have all the sequences, add them to each other to give the possible
        # combinations. Add them to the Array as well.

        for i in range(0, original_count - 1):
            for j in range(i + 1, original_count):
                new_seq = self._combination_list[i] + self._combination_list[j]
                self._combination_list.append(new_seq)
                sys.stdout.write(".")
                sys.stdout.flush()
        print
        print "We have " + str(len(self._combination_list)) + " combinations to compare"
        fasta_file = FASTA_File(self.filename + ".haplos")

        for each_seq in self._combination_list:
            fasta_file.add_sequence(each_seq)
        fasta_file.write_file()

    # --------------------------------------------------------------------------------
    #

    def find_match(self, sequence):

        print
        print "Analysing '" + sequence.title + "'"
        print "Start of sequence is  '" + sequence.sequence[:20] + "'"
        print "Start of reference is '" + self._haplotype_list[0].sequence[:20] + "'"

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
            print "Nothing matched! (even after clipping %d bases from the start" % MAX_CLIP_START


# --------------------------------------------------------------------------------

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


# --------------------------------------------------------------------------------

class Sequence:
    classIUPAC = IUPAC()

    def __init__(self, title, sequence):
        # self.IUPAC = IUPAC()
        self.title = title
        self.sequence = string.upper(sequence)
        self.parents = []

    def __add__(self, other):
        if len(self.sequence) != len(other.sequence):
            raise Exception("Can't add 2 sequences of differing length")

        newseq = ''
        for i in range(0, len(self.sequence)):
            newseq = newseq + self.classIUPAC.get_merged_code([self.sequence[i:i + 1], other.sequence[i:i + 1]])

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
        seq2 = cmp.sequence[:cmplen]
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
            optlist, self.remains = getopt.getopt(sys.argv[1:], 'a:o:cw', ['alignfile=', 'offset='])
        except getopt.GetoptError as err:
            print "There was a problem with the command line arguments - exiting"
            print str(err)
            sys.exit(1)

        #
        # Decode them so short arguments are converted to long arguments
        #
        short_to_long = {'a': "alignfile",
                         'o': "offset"}
        self.options = {}

        for arg in optlist:
            flag = arg[0]
            value = arg[1]

            while flag[0] == '-':
                flag = flag[1:]

            if flag in short_to_long:
                flag = short_to_long[flag]

            if flag in self.options:
                error_string = "Too many instances of the '" + flag + "' argument"
                raise Exception(error_string)

            self.options[flag] = value

        if 'offset' not in self.options.has_key:
            self.options["offset"] = 0

        self.options["offset"] = int(self.options["offset"])

        exit_needed = 0
        if 'c' in self.options:
            HaploFinder.copyright()
            exit_needed = 1
        if 'w' in self.options:
            HaploFinder.warranty()
            exit_needed = 1
        if exit_needed == 1:
            sys.exit(0)

        # flag the conditions and report any other requested information
        HaploFinder.start_details()
        if ('alignfile' not in self.options) or (len(self.remains) == 0):
            HaploFinder.help()

    @staticmethod
    def help():
        print "Usage: " + sys.argv[0] + " -a|--alignfile <alignmentfile>"
        print "            [-o|-offset <offset>] file1 [file2 ...]"
        print ""
        print "        Where <alignmentfile> is a padded alignment in FASTA file"
        print "                    format"
        print ""
        print "              <offset> is the offset of the sequences in the files to be"
        print "                    analysed, compared against the alignment. Positive"
        print "                    offsets mean the sequence starts within the alignment,"
        print "                    negative offsets (enclosed in quotes) mean the"
        print "                    sequence starts before the sequence in the alignment"
        print "                    file"
        print ""
        print "              file1 [file2 ...] are the names of FASTA formatted"
        print "                    sequence files to be analysed"
        sys.exit(0)

    def run(self):

        # Build the haplotype detector by reading in the file

        hd = HaplotypeDetector(self.options['alignfile'], self.options['offset'])

        for input_file in self.remains:
            fasta_file = FASTA_File(input_file)
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
        print ""
        print "Haplofinder version 1.2, Copyright (C) 2001,2016 Roslin Institute, Andy Law"
        print "Haplofinder comes with ABSOLUTELY NO WARRANTY; for details type"
        print "'" + sys.argv[0] + " -w'.  This is free software, and you are"
        print "welcome to redistribute it under certain conditions; type"
        print "'" + sys.argv[0] + " -c' for details."
        print ""

    @staticmethod
    def copyright():
        print "			    COPYING"
        print ""
        print "You may copy and distribute the Program (or a work based on it, under "
        print "Section 2 of the GPL) in object code or executable form under the terms of"
        print "Sections 1 and 2 of the GPL provided that you also do one of the following:"
        print ""
        print "    a) Accompany it with the complete corresponding machine-readable"
        print "    source code, which must be distributed under the terms of Sections"
        print "    1 and 2 above on a medium customarily used for software interchange; or,"
        print ""
        print "    b) Accompany it with a written offer, valid for at least three"
        print "    years, to give any third party, for a charge no more than your"
        print "    cost of physically performing source distribution, a complete"
        print "    machine-readable copy of the corresponding source code, to be"
        print "    distributed under the terms of Sections 1 and 2  on a medium"
        print "    customarily used for software interchange; or,"
        print ""
        print "    c) Accompany it with the information you received as to the offer"
        print "    to distribute corresponding source code.  (This alternative is"
        print "    allowed only for noncommercial distribution and only if you"
        print "    received the program in object code or executable form with such"
        print "    an offer, in accord with Subsection b above.)	"
        print ""
        print "The above is an excerpt of the full GPL license which you should have"
        print "received with this program. Please refer to the full document for"
        print "more details"

    @staticmethod
    def warranty():
        print "			    NO WARRANTY"
        print ""
        print "  11. BECAUSE THE PROGRAM IS LICENSED FREE OF CHARGE, THERE IS NO WARRANTY"
        print "FOR THE PROGRAM, TO THE EXTENT PERMITTED BY APPLICABLE LAW.  EXCEPT WHEN"
        print "OTHERWISE STATED IN WRITING THE COPYRIGHT HOLDERS AND/OR OTHER PARTIES"
        print "PROVIDE THE PROGRAM \"AS IS\" WITHOUT WARRANTY OF ANY KIND, EITHER EXPRESSED"
        print "OR IMPLIED, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF"
        print "MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.  THE ENTIRE RISK AS"
        print "TO THE QUALITY AND PERFORMANCE OF THE PROGRAM IS WITH YOU.  SHOULD THE"
        print "PROGRAM PROVE DEFECTIVE, YOU ASSUME THE COST OF ALL NECESSARY SERVICING,"
        print "REPAIR OR CORRECTION."
        print ""
        print "  12. IN NO EVENT UNLESS REQUIRED BY APPLICABLE LAW OR AGREED TO IN WRITING"
        print "WILL ANY COPYRIGHT HOLDER, OR ANY OTHER PARTY WHO MAY MODIFY AND/OR"
        print "REDISTRIBUTE THE PROGRAM AS PERMITTED ABOVE, BE LIABLE TO YOU FOR DAMAGES,"
        print "INCLUDING ANY GENERAL, SPECIAL, INCIDENTAL OR CONSEQUENTIAL DAMAGES ARISING"
        print "OUT OF THE USE OR INABILITY TO USE THE PROGRAM (INCLUDING BUT NOT LIMITED"
        print "TO LOSS OF DATA OR DATA BEING RENDERED INACCURATE OR LOSSES SUSTAINED BY"
        print "YOU OR THIRD PARTIES OR A FAILURE OF THE PROGRAM TO OPERATE WITH ANY OTHER"
        print "PROGRAMS), EVEN IF SUCH HOLDER OR OTHER PARTY HAS BEEN ADVISED OF THE"
        print "POSSIBILITY OF SUCH DAMAGES."
        print ""
        print "The above is an excerpt of the full GPL license which you should have"
        print "received with this program. Please refer to the full document for"
        print "more details"


if __name__ == '__main__':
    HaploFinder().run()
