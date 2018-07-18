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
from plishUtils import fetch_exonLen, fetch_seq, get_script_path
from Tkinter import Text, DISABLED, NORMAL, END

def _showLog(pr, txt):
    print pr
    if txt is not None:
      txt.config(state=NORMAL)
      txt.insert(END, pr + "\n")
      txt.see(END)
      txt.config(state=DISABLED)
      txt.update()

def main(inputId, db, debug=False, txt=None):
  if not os.path.exists(get_script_path() + '/database/' + db):
    _showLog('ERROR: Database ' + db + ' does not exist.', txt)
    return
  
  inputId = inputId.replace(' ', '')
  inputName, inputSeq = fetch_seq(inputId, db)
  
  if inputName is None:
    _showLog('ERROR: Transcript ID ' + inputId + ' does not exist.', txt)
    return
  
  result_fn = get_script_path() + "/results/" + inputId
  result_fn += "-" + inputName.replace('"', '') + "_hprobe.csv"
  blastdb = db
  
  _showLog('Target: ' + inputId + " (" + inputName + ")", txt)

  ###############################################################################
  # 1. Find anchors 
  ###############################################################################
  l = len(inputSeq)
  anc_index = [m.start() for m in finditer('(?:AG|TA)', inputSeq)]
  f = [x > 18 and x+19 < l for x in anc_index]
  anc_index = list(compress(anc_index, f))
  anc_seq = [inputSeq[x-19:x+21] for x in anc_index]
  
  _showLog('#Candidates: ' + str(len(anc_index)), txt)

  ###############################################################################
  # 2. Init hybridization probes
  ###############################################################################
  hprobes = [ Hprobe(anc_seq[i], inputId, inputName, anc_index[i]) for i in range(len(anc_index)) ]
  
  if debug:
    hprobes = hprobes[0:9]
  
  ###############################################################################
  # 3. Splice Junction Sites
  ###############################################################################
  _showLog('Step 1/5: Analyzing splice junction sites...', txt)
  
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
  # 4. Calculate melting temperature and filter
  ###############################################################################
  _showLog('Step 2/5: Calculating melting temperature...', txt)

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
  _showLog('Step 3/5: Calculating thermodynamics...', txt)
  
  for x in hprobes:
    x.do_fold(temperature=310.15) #310.15 K = 37 C
  
  ###############################################################################
  # 6. BLAST
  ###############################################################################
  _showLog('Step 4/5: Assessing specificity...', txt)
  
  for x in hprobes:
    x.do_blast(blastdb) #310.15 K = 37 C
  
  ###############################################################################
  # 7. Write CSV file
  ###############################################################################
  _showLog('Step 5/5: Writing result file...', txt)

  with open(result_fn, "w") as fh:
    h = 'Hprobe: Target seq (' + inputId + ' ' + inputName + ')\t'
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
      l = hp.seq + '\t'
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
  _showLog('Your results are available at: \n' + result_fn, txt)
