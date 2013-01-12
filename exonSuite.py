# ExonSuite: A program to find the optimized binding regions of the PUF protein.
#
# Copyright (C) <textyear> <name of author>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.
#
#
#
#
# Authors: Dilan Ustek, Abraham Kohrman, Bogdan Krstic, Karissa Fernandez
# Date: December 8th, 2011
# Skipped Exon Suite
# exonSuite.py

import datetime
import time

def runFullGenome(filename):
	"""Driver program for Skipped Exon Suite.
	Given a file containing annotated skipped exon data,
	returns a file indicating the 8mer with the lowest
	off-target score for each skipped exon."""
	print "Starting at " + getLocalTime()
	skips = skipsGen (filename)
	print "Finished skips_gen at " + getLocalTime()
	frequencies = buildFrequencies()
	frequencies = freqCount(skips, frequencies)
	tmp2 = bestExon(skips, frequencies)
	print "Finished bestexon at " + getLocalTime()

	# print "tmp2 is " + str(tmp2)
	tmp3 = fileMaker (tmp2)
	print "Finished fileMaker at " + getLocalTime()
return tmp3

def skipsGen(filename):
"""Returns a dictionary with the keys being the headers of each exon in the given file,
and the value for each key being the sequence contained in the file under each header."""
# initialize dictionary of skipped exons
skips = {}
# open the file, initialize header string
infile = open(filename)
header= ’’
# go through file and perform the following actions:
# for each occurence of ’>’ in the file,
# which signifies the header for each skipped exon,
# all lines until the next header are
# added to the string which is the value for the header.
# Also removes newlines and replaces
# them with empty characters.
for line in infile.readlines():
line = line.replace(’\n’, ’’)
if line[0] == ’>’:
header = line
elif skips.has_key(header):
skips[header] = skips[header] + line
else:
skips[header]=line
return skips
def buildFrequencies():
"""Builds a generic frequency dictionary for PUF binding 8mers (any 8mer
with a ’c’ or a ’t’ in the fifth region).
Keys are 8mers, values are initialized to zero so they can be updated
in freqCount.
Written by Samuel Rebelsky, November 2011."""
frequencies = {}
for n0 in [’a’,’c’,’g’,’t’]:
for n1 in [’a’, ’c’, ’g’, ’t’]:
for n2 in [’a’, ’c’, ’g’, ’t’]:
for n3 in [’a’, ’c’, ’g’, ’t’]:
for n4 in [ ’c’, ’t’]:
for n5 in [’a’, ’c’, ’g’, ’t’]:
for n6 in [’a’, ’c’, ’g’, ’t’]:
for n7 in [’a’, ’c’, ’g’, ’t’]:
mer = "".join([n0,n1,n2,n3,n4,n5,n6,n7])
frequencies[mer] = 0
return frequencies

def freqCount(skips, frequencies):
"""Counts the number of times each possible PUF binding 8mer appears
in the skipped exons. Returns an updated dictionary with exon identifiers
as keys and number of genome-wide occurences as values."""
for skip in skips:
skipkey = skips[skip]
top = len(skipkey) - 7
for i in range(top):
if frequencies.has_key(skipkey[i:i+8]):
frequencies[skipkey[i:i+8]] += 1
return frequencies

def bestExon(skips, frequencies):
"""Given a dictionary of ’skips’ (for each skip, the key is the name
of the skip and the value is the gene data) and a table of exon
frequencies, creates dictionary that gives the best exon for each
gene."""
#skips is a 2-d array I.E. a "list of lists" in python-speak.
#each element in the outer list contains a list
#which contains contains a header followed by the code for the exon
skipExons = []
cycle=1
geneName = {}
# cycle through skips, and find the best8mer for each skip
# update geneName to include the key ’skip’ whose value is
# the [best, count] list returned by best8mer for the skipped exon.
for skip in skips:
exon = skips[skip]
print "cycle " + str(cycle) + ": " + str(exon)
geneName[skip] = best8mer (exon, frequencies, cycle)
cycle+=1
#returns list of the best exon for each gene
return geneName

def best8mer(exon, frequencies, cycle):
"""Given an exon and a dictionary of 8mer-frequencies, find the least
frequently occuring 8-mer in the list of skips. Return that 8-mer
and its frequency."""
# build localfrequencies dictionary as a copy of the global frequencies
# dictionary to avoid updating frequencies dictionary master
# after best8mer is run
localfrequencies = {}
localfrequencies = dict(frequencies)
top = len(exon)-7
# Guess that the first valid exon is the best
j = 0
# Error message if PUF cannot bind to the exon.
# ’best’ and ’count’ values are modified later if PUF can bind.
best = "PUF cannot bind to this exon"
count = 0
# Find the location j of the first valid 8mer in the skipped exon
while ((not frequencies.has_key(exon[j:j+8])) and j < top):
j += 1
# If a valid 8mer is found, set ’best’ equal to it, and ’count’ equal to
# its initial frequency in the localfrequencies dictionary.
if frequencies.has_key(exon[j:j+8]):
best = exon[j:j+8]
count = localfrequencies[best]
# Loop through the rest of the exon, starting at the jth position,
# update localfrequencies, decrementing the value of the associated key
# every time an 8mer appears.
for q in range(j, top):
locmer = exon[q:q+8]
if frequencies.has_key(locmer):
localfrequencies[locmer] = localfrequencies[locmer] - 1
# If the value of locmer in localfrequencies
# is less than the current
# value for ‘count’, set ‘best’ to locmer,
# and ‘count’ to the localfrequencies
# value for locmer.
if localfrequencies[locmer] < count:
best = locmer
count = localfrequencies[locmer]
print "cycle " + str(cycle) + ", " + str(q) + " of " + str(top)
# And we’re done
return [best,count]

def fileMaker(results):
"""Given a dictionary of results, which has skipped exons as keys and
results from best8mer as values, filemaker returns human-readable output
in a filestream."""
today = datetime.date.today()
# Create a file for the data output table
table = open(’Skipped_Exon_Scan_’ + getLocalTimeFile(), ’w’)
# Make a list of the keys from the results
# (the headers of the skipped exons)
# and sort it
res_keys= []
res_keys = results.keys()
sorted(res_keys)
# Write a line for each key and the corresponding results
for key in res_keys:
table.write(key)
table.write(’ : ’)
table.write(str(results[key]))
table.write(’\n’)
table.close()
# Done!
return "Data table created. File name is Skipped_Exon_Scan_" + getLocalTimeFile()

def getLocalTime():
return str(time.localtime()[1]) + "/" + str(time.localtime()[2])
+ "/" + str(time.localtime()[0]) + " " + str(time.localtime()[3])
+ ":" + str(time.localtime()[4]) + ":" + str(time.localtime()[5])

def getLocalTimeFile():
return str(time.localtime()[1]) + "-" + str(time.localtime()[2])
+ "-" + str(time.localtime()[0]) + "_" + str(time.localtime()[3])
+ "-" + str(time.localtime()[4])

