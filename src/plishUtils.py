#!/usr/bin/python
###############################################################################
# Utility methods
# written by Daniel C. Ellwanger <dcellwanger.dev@gmail.com>
# last modified Jun 28 2018
###############################################################################
import os
from math import log10, log
from subprocess import call

###############################################################################
# Get tool path
###############################################################################
def get_script_path():
  pth = os.path.dirname(os.path.realpath(__file__))
  pth = os.path.abspath(os.path.join(pth, os.pardir))
  return pth
  #return os.path.dirname(os.path.realpath(sys.argv[0]))

###############################################################################
# Calculate reverse complement
###############################################################################
def revcompl(seq):
  rc = compl(seq[::-1])
  return rc

###############################################################################
# Calculate complement
###############################################################################
def compl(seq):
  rc = seq.upper(). \
       replace('A','t').replace('T','a'). \
       replace('C','g').replace('G','c')
  rc = rc.upper()
  return rc

###############################################################################
# Calculate melt temperature
# Basic equation: Marmur and Doty, J. Mol. Biol  1962
# Salt adjusted equation: Nakano et. al., PNAS 1999
# Formaldehyde equation: Meinkoth and Wahl, Anal Biochem 1984 
# https://academic.oup.com/bioinformatics/article/21/6/711/199347
# http://biotools.nubic.northwestern.edu/OligoCalc.html#helpadjusted
###############################################################################
def melttemp(seq, c_salt, p_formamide): 
  l = len(seq)
  wA = float(seq.count('A'))
  xT = float(seq.count('T'))
  yG = float(seq.count('G'))
  zC = float(seq.count('C'))
  if c_salt == None and p_formamide == None:
    tm = 64.9+41*(yG+zC-16.4)/(wA+xT+yG+zC)
  elif p_formamide == None:
    tm = 100.5+(41*(yG+zC)/(wA+xT+yG+zC))
    tm = tm-(820/(wA+xT+yG+zC))+16.6*log10(c_salt)
  else:
    tm = 81+16.6*log(c_salt)+0.41*(yG/l+zC/l) 
    tm = tm-500/l-0.61*p_formamide
  return tm
  
###############################################################################
# Returns free energy of RNA fold and duplices
###############################################################################
def fold(infile, temperature):
  devnull = open(os.devnull, 'w')
  exe = get_script_path() + "/tools/RNAstructure/exe/oligoscreen "
  outfile = infile + ".out"
  t = "--temperature " + str(temperature)
  cmd = exe + infile + " " + outfile + " " + t
  call([cmd], shell=True, stdout=devnull, stderr=devnull)
  with open(outfile) as f:
    content = f.readlines()
  f.close()
  dg = [x.strip() for x in content]
  os.remove(infile)
  os.remove(outfile)
  devnull.close()
  return(dg)

###############################################################################
# Returns exon lengths for transcript id (considering strand)
############################################################################### 
def fetch_exonLen(transcript_id, db):
  filepath = get_script_path() + '/database/' + db + '/' + db + '.exons'
  exon_lens = list()
  transcript_id = transcript_id.split('.')[0]
  
  with open(filepath) as fh:
    for line in fh:
      edat = line.split("\t")
      edat_id = edat[0].split(".")[0]
      if transcript_id == edat_id:
        exon_lens = edat[1].split(",")
        exon_lens = list(map(int, exon_lens))
        break
    fh.close()
  return exon_lens

###############################################################################
# Returns sequence and gene name for transcript id (considering strand)
############################################################################### 
def fetch_seq(transcript_id, db):
  filepath = get_script_path() + '/database/' + db + '/' + db + '.fa'
  exon_lens = list()
  transcript_id = transcript_id.split('.')[0]
  seq = name = None
  
  with open(filepath) as fh:
    match = False
    for line in fh:
      if line.startswith(">"):
        sdat = line.split("|")
        sdat_id = sdat[0][1:].split(".")[0]
        if transcript_id == sdat_id:
          match = True
          name = sdat[1].strip()
      elif match:
        seq = line.strip().upper()
        break
    fh.close()
  return name, seq
  
###############################################################################
# Returns BLAST results
############################################################################### 
def blast(infile, db):#, taxonid): #nr
  db = get_script_path() + "/database/" + db + "/" + db
  devnull = open(os.devnull, 'w')
  exe = get_script_path() + "/tools/ncbi-blast/bin/blastn "
  outfile = infile + ".blast"
  #org = '-entrez_query "txid' + str(taxonid) + ' [ORGN]" '
  query = '-query ' + infile + ' '  
  task = '-task megablast '
  db = '-db ' + db + ' '
  out = '-out ' + outfile + ' '
  outfmt = '-outfmt 6 '
  evalue = '-evalue 10 '
  cmd = exe + query + task + db + out + outfmt + evalue #+ org + '-remote'
  call([cmd], shell=True, stdout=devnull, stderr=devnull)
  with open(outfile) as f:
    content = f.readlines()
  f.close()
  res = [x.strip() for x in content]
  os.remove(infile)
  os.remove(outfile)
  return(res)

###############################################################################
# Fetch database info
############################################################################### 
def fetch_dbInfo():
  dbdir = get_script_path() + '/database/'
  dbsub = next(os.walk(dbdir))[1]
  dbnames = list()
  dbids = list()
  for db in dbsub:
    with open(dbdir + db + '/' + db + '.info') as f:
      content = f.read().splitlines()
    for line in content:
      if line.startswith("dbname"):
        dbnames.append(line.split("=")[1])
      elif line.startswith("dbid"):
        dbids.append(line.split("=")[1])
  return dbnames, dbids
