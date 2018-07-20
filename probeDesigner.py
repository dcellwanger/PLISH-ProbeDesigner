#!/usr/bin/python
###############################################################################
# Hybridization probe designer for PLISH
# https://elifesciences.org/articles/30510
# written by Daniel C. Ellwanger <dcellwanger.dev@gmail.com>
# last modified Jun 28 2018
###############################################################################
import os, sys
from argparse import ArgumentParser
from Tkinter import Tk
# import modules in src/
sys.path.insert(0, os.path.dirname(os.path.realpath(__file__)) + '/src')
from plishGUI import GUI
from plishUtils import get_script_path, fetch_dbInfo
from plishMain import main

# Run tool in debug mode? For devel use only.
_debug = True

###############################################################################
# Set environment vars
###############################################################################
os.environ['DATAPATH'] = get_script_path() + '/tools/RNAstructure/data_tables/'
os.environ['BLASTDB'] = get_script_path() + '/database/'

###############################################################################
# Fetch dbs
###############################################################################
dbnames, dbids = fetch_dbInfo()

###############################################################################
# Read arguments and RUN!
###############################################################################
if len(sys.argv) == 1: #run w/ GUI
  root = Tk()
  my_gui = GUI(root, dbnames, dbids, _debug)
  root.mainloop()
else: #run via command line
  parser = ArgumentParser()
  parser.add_argument('-db', '--database', dest='db', required=True,
                      help='database to use: ' + str(dbids), metavar='ID')
  parser.add_argument('-tx', '--transcript', dest='tx', required=True,
                      help='transcript database id', metavar='ID')
  args = parser.parse_args()
  db = args.db
  inputId = args.tx
  main(inputId, db, debug=_debug)
  write_probesCSV(inputId, inputName, hprobes, results_fn)
