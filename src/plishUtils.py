#!/usr/bin/python
###############################################################################
# Utility methods
# written by Daniel C. Ellwanger <dcellwanger.dev@gmail.com>
# last modified Jul 18 2018
###############################################################################
import os
from math import log10, log
from subprocess import call
from Tkinter import Text, DISABLED, NORMAL, END
from itertools import compress

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
        exon_lens = edat[6].split(",")
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

###############################################################################
# Shows progress log
############################################################################### 
def show_log(pr, txt):
    print pr
    if txt is not None:
      txt.config(state=NORMAL)
      txt.insert(END, pr + "\n")
      txt.see(END)
      txt.config(state=DISABLED)
      txt.update()
      
###############################################################################
# Filter probes
###############################################################################    
def filter_probes(hprobes, spec, mingc, multiexon, mintm, maxtm, 
                  mindimer, minfold, maxduplex):     
  if multiexon:
    f = [x.exons[0] != x.exons[1] for x in hprobes]
    hprobes = list(compress(hprobes, f))
  
  specCat = ['none','gene','isoform']  
  d = dict([(y,x+1) for x,y in enumerate(set(specCat))])

  f = [d[x.spec] >= d[spec] and
       x.gc > mingc and
       mintm < x.larm.tm < maxtm and
       mintm < x.rarm.tm < maxtm and
       x.larm.dg_bimol > mindimer and
       x.rarm.dg_bimol > mindimer and
       x.larm.dg_uimol > minfold and
       x.rarm.dg_uimol > minfold and
       x.larm.dg_duplex < maxduplex and
       x.rarm.dg_duplex < maxduplex for x in hprobes]
  return list(compress(hprobes, f))
      
###############################################################################
# Write result CSV file
###############################################################################
def write_probesCSV(inputId, inputName, hprobes, txt):
  inputName = inputName.replace('"', '')
  result_fn = get_script_path() + "/results/" + inputName
  result_fn += "-" + inputId + "_hprobe.csv"
  
  n_probes = str(len(hprobes))
  show_log('Writing result file for ' + n_probes + ' filtered probes ...', txt)
  with open(result_fn, "w") as fh:
    h = 'Hprobe: Id\t'
    h += 'Hprobe: Target sequence\t'
    h += 'Hprobe: %GC\t'
    h += 'Hprobe: Multipe exons?\t'
    h += 'Hprobe: Exons\t'
    h += 'Hprobe: Specificity\t'
    h += 'Hprobe: Blast Hits (Ident%)\t'
    h += 'Left: Seq\t' + 'Left: Tm\t'
    h += 'Left: Bimol.\t' + 'Left: Unimol.\t' + 'Left: Duplex\t'
    h += 'Left: Open5\t' + 'Left: Open3\t'
    h += 'Right: Seq\t' + 'Right: Tm\t'
    h += 'Right: Bimol.\t' + 'Right: Unimol.\t' + 'Right: Duplex\t'
    h += 'Right: Open5\t' + 'Right: Open3'
    fh.write(h + '\n')
    for hp in hprobes:
      l = inputName + '-' + inputId + '-' + str(hp.index) + '\t'
      l += hp.seq + '\t'
      l += str(round(hp.gc, 1)) + '\t'
      l += str(hp.exons[0] != hp.exons[1]) + '\t'
      l += str(hp.exons).replace('[', '').replace(']', '') + '\t'
      l += hp.spec + '\t'
      l += str(hp.bhits).replace('[', '').replace(']', '') + '\t'
      l += hp.larm.seq + '\t'
      l += str(hp.larm.tm) + '\t'
      l += str(hp.larm.dg_bimol) + '\t'
      l += str(hp.larm.dg_uimol) + '\t'
      l += str(hp.larm.dg_duplex) + '\t'
      l += str(hp.larm.dg_2bpat5) + '\t'
      l += str(hp.larm.dg_2bpat3) + '\t'
      l += hp.rarm.seq + '\t'
      l += str(hp.rarm.tm) + '\t'
      l += str(hp.rarm.dg_bimol) + '\t'
      l += str(hp.rarm.dg_uimol) + '\t'
      l += str(hp.rarm.dg_duplex) + '\t'
      l += str(hp.rarm.dg_2bpat5) + '\t'
      l += str(hp.rarm.dg_2bpat3) + '\t'
      fh.write(l + '\n')
  fh.close()
  show_log('Your probe information is available at: \n' + result_fn, txt)
  return result_fn
  
###############################################################################
# Write result FASTA file
###############################################################################
def write_probesFNA(inputId, inputName, hprobes, txt):
  inputName = inputName.replace('"', '')
  result_fn = get_script_path() + "/results/" + inputName
  result_fn += "-" + inputId + "_hprobe.fna"
  
  # connector and bridge sequences
  hl = { '2X' : 'TCGTACGTCTAACTTACGTCGTTATG',
         '3X' : 'TATTCGTTCGAACTTACGTCGTTATG',
         '4X' : 'TTAGTAGGCGAACTTACGTCGTTATG',
         '5X' : 'TAGCGCTAACAACTTACGTCGTTATG',
         '6X' : 'TAGGTCAGGAAACTTACGTCGTTATG'}
  hr = { '2X' : 'TTATACGTCGAGTTGAAGAACAACCTG',
         '3X' : 'TTATACGTCGAGTTGACCGACGTATTG',
         '4X' : 'TTATACGTCGAGTTGAACATAAGTGCG',
         '5X' : 'TTATACGTCGAGTTGAACGTCGTAACA',
         '6X' : 'TTATACGTCGAGTTGAATAGCCAGGTT'}
  
  hprobeId = inputName + '-' + inputId
  with open(result_fn, "w") as fh:
    for hp in hprobes:
      probeName = hprobeId + '-' + str(hp.index)
      for x in sorted(hl.keys()):
        lhead = '>' + 'HL' + x + '-' + probeName
        lseq = hl[x] + hp.larm.seq
        rhead = '>' + 'HR' + x + '-' + probeName
        rseq = hp.rarm.seq + hr[x]
        fh.write(lhead + '\n' + lseq + '\n')
        fh.write(rhead + '\n' + rseq + '\n')
  fh.close()
  show_log('Your probe sequences are available at: \n' + result_fn, txt)  
  return result_fn
