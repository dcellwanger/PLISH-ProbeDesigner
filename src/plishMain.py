#!/usr/bin/python
###############################################################################
# Entry point of probe detection and feature computation
# written by Daniel C. Ellwanger <dcellwanger.dev@gmail.com>
# last modified Jun 28 2018
###############################################################################
import os, sys
from re import findall, finditer
from itertools import compress
from plishHprobe import Hprobe
from plishUtils import fetch_exonLen, fetch_seq, get_script_path, show_log

def main(inputId, db, debug=False, txt=None):
  if not os.path.exists(get_script_path() + '/database/' + db):
    show_log('ERROR: Database ' + db + ' does not exist.', txt)
    return None, None, None
  
  inputId = inputId.replace(' ', '')
  inputName, inputSeq = fetch_seq(inputId, db)
  
  if inputName is None:
    show_log('ERROR: Transcript ID ' + inputId + ' does not exist.', txt)
    return None, None, None
  
  blastdb = db
  
  show_log('Target: ' + inputId + " (" + inputName + ")", txt)

  ###############################################################################
  # 1. Find anchors 
  ###############################################################################
  l = len(inputSeq)
  anc_index = [m.start() for m in finditer('(?:AG|TA)', inputSeq)]
  f = [x > 18 and x+19 < l for x in anc_index]
  anc_index = list(compress(anc_index, f))
  anc_seq = [inputSeq[x-19:x+21] for x in anc_index]
  
  show_log('#Candidates: ' + str(len(anc_index)), txt)

  ###############################################################################
  # 2. Init hybridization probes
  ###############################################################################
  hprobes = [ Hprobe(anc_seq[i], inputId, inputName, anc_index[i]) for i in range(len(anc_index)) ]
  
  if debug:
    hprobes = hprobes[0:9]
  
  ###############################################################################
  # 3. Splice Junction Sites
  ###############################################################################
  show_log('Step 1/4: Analyzing splice junction sites...', txt)
  
  elen = fetch_exonLen(inputId, db)
  
  for x in hprobes:
    hpStart_exon = -1
    hpEnd_exon = -1
    csum = 0
    for i in range(len(elen)):
      if(x.start > csum):
        hpStart_exon = i + 1
      if(x.end > csum):
        hpEnd_exon = i + 1
      csum = csum + elen[i]
    x.exons = [hpStart_exon, hpEnd_exon]
  
  #f = [x.exons[0] != x.exons[1] for x in hprobes]
  #hprobes = list(compress(hprobes, f))
  #print('#Filtered by splice junction: ' + str(len(hprobes)))
  
  ###############################################################################
  # 4. Calculate melting temperature
  ###############################################################################
  show_log('Step 2/4: Calculating melting temperature...', txt)

  for x in hprobes:
    x.calc_tm(c_salt=0.05, p_formamide=None) #1, 0.5
  
  # Filter
  #fL = [x.larm.tm > 55.0 and x.larm.tm < 65.0 for x in hprobes] #45 - 65
  #fR = [x.rarm.tm >= 55.0 and x.rarm.tm <= 65.0 for x in hprobes]
  #f = [all(tup) for tup in zip(fL, fR)]
  #hprobes = list(compress(hprobes, f))
  #print('#Filtered by melting temp.: ' + str(len(hprobes)))
  
  ###############################################################################
  # 5. Calculate thermodynamic features
  ###############################################################################
  show_log('Step 3/4: Calculating thermodynamics...', txt)
  
  for x in hprobes:
    x.do_fold(temperature=310.15) #310.15 K = 37 C
  
  ###############################################################################
  # 6. BLAST
  ###############################################################################
  show_log('Step 4/4: Assessing specificity...', txt)
  
  for x in hprobes:
    x.do_blast(blastdb)
    
  ###############################################################################
  # Return
  ###############################################################################
  show_log('------------------[ DONE ]------------------', txt)
  return inputId, inputName, hprobes
