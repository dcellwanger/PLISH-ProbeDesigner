#!/usr/bin/python
###############################################################################
# Database prep for PLISH probe designer
# written by Daniel C. Ellwanger <dcellwanger.dev@gmail.com>
# last modified Jun 28 2018
# !! Has TODOs !!
###############################################################################
from re import findall
from plishUtils import compl

###############################################################################
# Global variables (need to be adjusted)
# TODO: Replace with command line arguments
###############################################################################
# Input: gtf file and matching genome sequence fasta file
gtf_fn = 'annotation.gtf' 
genomeSeq_fn = 'genomeseq.fa'
refseq_fa = 'gga_refseq_NCBI.fa'
refseq_gff3 = 'gga_refseq_NCBI.gff3'

# Output
exon_fn = 'refseq_gga.exons'
txSeq_fn = 'refseq_gga.fa'

###############################################################################
# Generate exon coodinate file for each transcript
############################################################################### 
def exon_coordinates(gtf_fn, out_fn): 
  print "Extracting exon info ..."
  txs2elow = {}
  txs2eup = {}
  txs2strand = {}
  txs2chrom = {}
  txs2name = {}
  lcount = 0
  with open(gtf_fn) as fh:
    for line in fh:
      lcount = lcount + 1
      if 'exon' in line:
        tx_id = findall('transcript_id "(.+?)";|$', line)[0]
        tx_name = findall('gene_name "(.+?)";|$', line)[0]
        if tx_name == '':
          tx_name = findall('gene_id "(.+?)";|$', line)[0]
        edat = line.split("\t")
        chrom = edat[0]
        elow = int(edat[3])
        eup = int(edat[4])
        strand = edat[6]
        if tx_id in txs2name:
          txs2elow[tx_id].append(elow)
          txs2eup[tx_id].append(eup)
        else:
          txs2elow[tx_id] = [elow]
          txs2eup[tx_id] = [eup]
          txs2strand[tx_id] = strand
          txs2chrom[tx_id] = chrom
          txs2name[tx_id] = tx_name
      if lcount % 10000 == 0:
        print "Processed " + str(lcount) + " lines ..."
  fh.close()

  print "Calculating exon lengths ..."
  txs2elen = {}
  for tx_id in txs2name.keys():
    strand = txs2strand[tx_id]
    txs2elow[tx_id].sort(reverse=(strand == "-"))
    txs2eup[tx_id].sort(reverse=(strand == "-"))
    if strand == "-":
      starts = txs2eup[tx_id]
      ends = [-x for x in txs2elow[tx_id]]
    else:
      starts = [-x for x in txs2elow[tx_id]]
      ends = txs2eup[tx_id]
    txs2elen[tx_id] = [sum(x)+1 for x in zip(starts, ends)]
    
  print "Writing file ..."
  with open(out_fn, "w") as fh:
    for tx_id in txs2name.keys():
      elow = txs2elow[tx_id]
      elow = str(elow).replace('[', '').replace(']', '').replace(' ', '')
      eup = txs2eup[tx_id]
      eup = str(eup).replace('[', '').replace(']', '').replace(' ', '')
      elen = txs2elen[tx_id]
      elen = str(elen).replace('[', '').replace(']', '').replace(' ', '')
      strand = txs2strand[tx_id]
      chrom = txs2chrom[tx_id]
      name = txs2name[tx_id]
      fh.write(tx_id + '\t' + '"' + name + '"' + '\t' + chrom + '\t' + 
              strand + '\t' + elow + '\t' + eup + '\t' + elen + '\n')
  fh.close()
      
###############################################################################
# Generate fasta sequence file for each transcript
############################################################################### 
def tx_sequences(genomeSeq_fn, exon_fn, out_fn):
  def einfo(eline, cseq):
    eseq = ''
    edat = eline.split("\t")
    txid = edat[0]
    txname = edat[1]
    strand = edat[3]
    elow = edat[4].split(",")
    eup = edat[5].split(",")
    for i in range(len(elow)):              
      s_from = int(elow[i])
      s_to = int(eup[i])
      if strand == '-': #minus strand
        eseq = eseq + cseq[s_to-1:s_from-2:-1]
      else: #plus strand
        eseq = eseq + cseq[s_from-1:s_to:1]
    if strand == '-':
      eseq = compl(eseq)
    return txid, txname, eseq
  
  chr2info = {}
  with open(exon_fn) as fh:
    for line in fh:
      chrom = line.split("\t")[2]
      if chrom in chr2info:
        chr2info[chrom].append(line)
      else:
        chr2info[chrom] = [line]
  fh.close()

  chrom = None
  cseq = ''
  with open(out_fn, 'w') as fout, open(genomeSeq_fn) as fin:
    for line in fin:
      if line.startswith(">") and chrom != None:
        if chrom in chr2info:
          for eline in chr2info[chrom]:
            txid, txname, eseq = einfo(eline, cseq)
            fout.write('>' + txid + '|' + txname + '\n' + eseq + '\n')
        #reset  
        cseq = ''
        chrom = line[1:].strip()
        
      elif line.startswith(">") and chrom == None:
        chrom = line[1:].strip()
      else:
        cseq = cseq + line.strip()
    
    if chrom in chr2info: #last chr seq
      for eline in chr2info[chrom]:
        txid, txname, eseq = einfo(eline, cseq)
        fout.write('>' + txid + '|' + txname + '\n' + eseq + '\n')
  fin.close()    
  fout.close()

###############################################################################
# Generate fasta sequence file for each transcript from NCBI export
############################################################################### 
def transform_refseq_NCBI(refseq_fa, out_fn):
  head = None
  seq = ''
  with open(out_fn, 'w') as fout, open(refseq_fa) as fin:
    for line in fin:
      if line == '':
        continue
      elif line.startswith(">"):
        if head != None:
          if flag == True:
            fout.write(head + "\n" + seq + "\n")
          seq = ''
          
        ldat = line.split(" ")
        head = ldat[0]
        flag = False
        for l in ldat:
          if ")," in l:
            head += l.replace('(', '|"').replace('),', '"')
            flag = True

      else:
        seq += line.strip()
    if flag == True:
      fout.write(">" + head + "\n" + seq + "\n")
  fout.close()
  fin.close()

###############################################################################
# Generate exon file for each transcript from NCBI export
############################################################################### 
def transform_exon_NCBI(refseq_gff3, out_fn):
  txs2elow = {}
  txs2eup = {}
  txs2strand = {}
  with open(refseq_gff3) as fh:
    for line in fh:
      if line.startswith('#') or line == '\n':
        continue
      else:
        ldat = line.split('\t')
        entry = ldat[2]
        annos = ldat[8].split(";")
        acc = None
        for anno in annos:
          if 'transcript_id' in anno:
            acc = anno.split("=")[1].strip()
        if entry == 'exon' and acc is not None:
          lo = int(ldat[3])
          up = int(ldat[4])
          if acc in txs2elow:
            txs2elow[acc].append(lo)
            txs2eup[acc].append(up)
          else:
            txs2elow[acc] = [lo]
            txs2eup[acc] = [up]
            txs2strand[acc] = ldat[6]
  fh.close()
  print("Found information for " + str(len(txs2strand.keys())) + " transcripts.")

  print "Calculating exon lengths ..."
  txs2elen = {}
  for tx_id in txs2strand.keys():
    strand = txs2strand[tx_id]
    txs2elow[tx_id].sort(reverse=(strand == "-"))
    txs2eup[tx_id].sort(reverse=(strand == "-"))
    if strand == "-":
      starts = txs2eup[tx_id]
      ends = [-x for x in txs2elow[tx_id]]
    else:
      starts = [-x for x in txs2elow[tx_id]]
      ends = txs2eup[tx_id]
    txs2elen[tx_id] = [sum(x)+1 for x in zip(starts, ends)]
    
  print "Writing file ..."
  with open(out_fn, "w") as fh:
    for tx_id in txs2strand.keys():
      elow = txs2elow[tx_id]
      elow = str(elow).replace('[', '').replace(']', '').replace(' ', '')
      eup = txs2eup[tx_id]
      eup = str(eup).replace('[', '').replace(']', '').replace(' ', '')
      elen = txs2elen[tx_id]
      elen = str(elen).replace('[', '').replace(']', '').replace(' ', '')
      strand = txs2strand[tx_id]
      fh.write(tx_id + '\t' + elen + '\n')
  fh.close()

###############################################################################
# RUN
# TODO: Automate file generation; write *.info file.
############################################################################### 
#exon_coordinates(gtf_fn, exon_fn)
#tx_sequences(genomeSeq_fn, exon_fn, txSeq_fn)
#transform_refseq_NCBI(refseq_fa, txSeq_fn)
#transform_exon_NCBI(refseq_gff3, exon_fn)

# Generate BLAST db
#from plishUtils import get_script_path
#dbid = 'refseq_gga'
#cmd = get_script_path() + '/tools/ncbi-blast/bin/makeblastdb ' + \
#  '-in ' + get_script_path() + '/database/' + dbid + '/' + dbid + '.fa ' + \
#  '-dbtype nucl ' + \
#  '-out ' + get_script_path() + '/database/' + dbid + '/' + dbid
#print cmd  
