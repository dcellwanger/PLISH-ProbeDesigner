#!/usr/bin/python
###############################################################################
# Database prep for PLISH probe designer
# written by Daniel C. Ellwanger <dcellwanger.dev@gmail.com>
# last modified Jul 18 2018
###############################################################################
import os, sys, re
#from re import compile, match, group
from argparse import ArgumentParser
from subprocess import call
# import modules in src/
sys.path.insert(0, os.path.dirname(os.path.realpath(__file__)) + '/src')
from plishUtils import get_script_path
from plishDbUtils import write_infoFile, write_exonFile, write_sequenceFile
from plishDbUtils import generate_BLASTdb

###############################################################################
# Set environment vars
###############################################################################
os.environ['BLASTDB'] = get_script_path() + '/database/'

###############################################################################
# Fetch arguments
###############################################################################
parser = ArgumentParser()
parser.add_argument('-gff', dest='gff', required=True,
                    help='annotation GFF file', metavar='FILEPATH')
parser.add_argument('-fna', dest='fna', required=True,
                    help='genome sequence FASTA file', metavar='FILEPATH')
parser.add_argument('-db', dest='db', required=True,
                    help='identifier of database (e.g., mmu_refseq); ' + \
                    'please, avoid white-spaces and special characters.', metavar='ID')
parser.add_argument('-name', dest='name', required=True,
                    help='name of database', metavar='NAME')
parser.add_argument('-comment', dest='comment', required=False,
                    help='any comment to add to the info file (e.g., genome assembly)', 
                    metavar='COMMENT')
args = parser.parse_args()

###############################################################################
# Input
###############################################################################
gff_fn = args.gff
fna_fn = args.fna
db_id =  args.db
db_name = args.name
db_comment = args.comment
if db_comment is None:
  db_comment = ''

###############################################################################
# Output
###############################################################################
exon_fn = get_script_path() + '/database/' + db_id + '/' + db_id + '.exons'
txSeq_fn = get_script_path() + '/database/' + db_id + '/' + db_id + '.fna'
info_fn = get_script_path() + '/database/' + db_id + '/' + db_id + '.info'
blastdb = get_script_path() + '/database/' + db_id + '/' + db_id

###############################################################################
# RUN
###############################################################################
write_infoFile(info_fn, db_id, db_name, db_comment)
write_exonFile(gff_fn, exon_fn)
write_sequenceFile(fna_fn, exon_fn, txSeq_fn)
generate_BLASTdb(blastdb, txSeq_fn)
print 'Generation of database "' + db_id + '" is finished.'
