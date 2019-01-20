#!/usr/bin/python
###############################################################################
# Container for hybridization probe
# written by Daniel C. Ellwanger <dcellwanger.dev@gmail.com>
# last modified Jun 28 2018
###############################################################################
from plishUtils import revcompl, melttemp, get_script_path, fold, blast

###############################################################################
# Public container for hybridization probe data
# Attributes
# start:    start index of probe in transcript sequence (int)
# end:      end index of probe in transcript sequence (int)
# seq:      probe sequence (str)
# seq_id:   target sequence id
# seq_name: target sequence/gene name
# larm:     left arm of probe (_HprobeArm)
# rarm:     right arm of probe (_HprobeArm)
# exons:    exons that are spanned by probe (list)
# bhits:    blast hits <'acc|"name":identity'>
###############################################################################
class Hprobe:
  # Constructor
  def __init__(self, seq, seq_id, seq_name, index):
    seq_rc = revcompl(seq)
    self.start = index - 19
    self.end = index + 21
    self.rarm = _HprobeArm(seq_rc[:20])
    self.larm = _HprobeArm(seq_rc[20:])
    self.seq = seq
    self.seq_id = seq_id.split('.')[0]
    self.seq_name = seq_name
    self.index = index
    self.gc = (seq.count("G") + seq.count("C")) * 100.0 / len(seq) 

  # Attributes
  start = end = rarm = larm = index = gc = exons = bhits = spec = None
  seq = seq_id = seq_name = None
  
  # Methods
  # Calculate melting temperature
  def calc_tm(self, c_salt, p_formamide):
    self.larm.tm = round(melttemp(self.larm.seq, c_salt, p_formamide), 1)
    self.rarm.tm = round(melttemp(self.rarm.seq, c_salt, p_formamide), 1)

  # Calculate folding and duplices
  def do_fold(self, temperature=310.15): #310.15 K = 37 C
    infile = get_script_path() + "/tmp/" + str(id(self)) + ".seq"
    with open(infile, "w") as f:
      f.write(self.larm.seq + "\n" + self.rarm.seq) 
    f.close()
    dg = fold(infile, temperature)
    larm_dg = dg[1].split("\t")
    self.larm.dg_bimol = float(larm_dg[1])
    self.larm.dg_uimol = float(larm_dg[2])
    self.larm.dg_duplex = float(larm_dg[3])
    self.larm.dg_2bpat5 = float(larm_dg[4])
    self.larm.dg_2bpat3 = float(larm_dg[5])
    rarm_dg = dg[2].split("\t")
    self.rarm.dg_bimol = float(rarm_dg[1])
    self.rarm.dg_uimol = float(rarm_dg[2])
    self.rarm.dg_duplex = float(rarm_dg[3])
    self.rarm.dg_2bpat5 = float(rarm_dg[4])
    self.rarm.dg_2bpat3 = float(rarm_dg[5])
  
  # Calculate blast hits
  def do_blast(self, blastdb):
    infile = get_script_path() + "/tmp/" + str(id(self)) + ".seq"
    with open(infile, "w") as f:
      f.write(self.seq) 
    f.close()
    hits = blast(infile, blastdb)
    self.spec = 'isoform'
    bhits = list()
    if(len(hits) == 0):
      self.spec = 'none'
    else:
      for h in hits:
        hdat = h.split("\t")
        acc = hdat[1]
        ide = hdat[2]
        
        adat = acc.split('|')
        txid = adat[0].split('.')[0]
        symbol = adat[1]
        bhits.append(txid + "|" + symbol + ":" + ide)
        
        if txid != self.seq_id:
          self.spec = 'gene'
        if symbol != self.seq_name:
          self.spec = 'none'
    self.bhits = bhits  

###############################################################################
# Private container for hybridization probe arm
# Attributes
# seq:       probe arm sequence (str)
# tm:        melting temp (float)
# dg_bimol:  bimolecular folding free energy change for two identical
#            oligonucleotides interacting (float)
# dg_uimol:  folding free energy change for unimolecular 
#            self-structure (float)
# dg_duplex: folding free energy change for duplex formation with the 
#            complementary sequence (float)
# dg_2bpat5: cost of opening the two base pairs at the 5' end of the 
#            oligonucleotide in a duplex with the complementary 
#            sequence (float)
# dg_2bpat3: cost of opening the two base pairs at the 3' end of the 
#            oligonucleotide in a duplex with the complementary 
#            sequence (float)
###############################################################################
class _HprobeArm:
  # Constructor
  def __init__(self, seq):
    self.seq = seq
    
  #Attributes
  seq = tm = dg_bimol = dg_uimol = dg_duplex = dg_2bpat5 = dg_2bpat3 = None
