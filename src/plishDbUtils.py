#!/usr/bin/python
###############################################################################
# Database prep utils for PLISH probe designer
# written by Daniel C. Ellwanger <dcellwanger.dev@gmail.com>
# last modified Jul 18 2018
###############################################################################
import os, sys, re
from argparse import ArgumentParser
from subprocess import call
# import modules in src/
sys.path.insert(0, os.path.dirname(os.path.realpath(__file__)) + '/src')
from plishUtils import compl, get_script_path

###############################################################################
# Generate info file
###############################################################################
def write_infoFile(info_fn, db_id, db_name, db_comment):
  print 'Writing info file...'
  if not os.path.exists(os.path.dirname(info_fn)):
    try:
      os.makedirs(os.path.dirname(info_fn))
    except OSError as exc: # guard against race condition
      if exc.errno != errno.EEXIST:
        raise
  with open(info_fn, 'w') as fh:
    fh.write('#' + db_comment + '\n')
    fh.write('dbid=' + db_id + '\n')
    fh.write('dbname=' + db_name + '\n')
  fh.close()

###############################################################################
# Generate exon coodinate file for each transcript
###############################################################################
def write_exonFile(gff_fn, out_fn):
  print "Extracting exon info..."
  txs2elow = {}
  txs2eup = {}
  txs2strand = {}
  txs2chrom = {}
  txs2name = {}
  lcount = 0
  pid = re.compile('.*transcript_id=(.+?);') #transcript_id=(.+);?
  pname = re.compile('.*gene(?:_id)?=(.+?);')
  #pname2 = re.compile('.*gene_id=(.+);?)
  with open(gff_fn) as fh:
    for line in fh:
      if line.startswith('#'):
        continue
      lcount = lcount + 1
      ldata = line.split('\t')
      line = line.strip() + ";"
      if ldata[2] == 'exon':
        tx_id = pid.match(line)
        if tx_id is None:
          continue
        tx_id = tx_id.group(1)
        tx_name = pname.match(line).group(1)
        #tx_id = findall('transcript_id=(.+);?|$', line)[0]
        #tx_name = findall('gene=(.+?);?|$', line)[0]
        #if tx_name == '':
        #tx_name = findall('gene_id=(.+?);?|$', line)[0]
        chrom = ldata[0]
        elow = int(ldata[3])
        eup = int(ldata[4])
        strand = ldata[6]
        if tx_id in txs2name:
          txs2elow[tx_id].append(elow)
          txs2eup[tx_id].append(eup)
        else:
          txs2elow[tx_id] = [elow]
          txs2eup[tx_id] = [eup]
          txs2strand[tx_id] = strand
          txs2chrom[tx_id] = chrom
          txs2name[tx_id] = tx_name
      if lcount % 500000 == 0:
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
  
  print "Writing exon file ..."
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
def write_sequenceFile(genomeSeq_fn, exon_fn, out_fn):
  print 'Writing sequence file ...'
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
      if line.startswith(">") and chrom is not None:
        if chrom in chr2info:
          for eline in chr2info[chrom]:
            txid, txname, eseq = einfo(eline, cseq)
            fout.write('>' + txid + '|' + txname + '\n' + eseq + '\n')
          #reset
          cseq = ''
          chrom = line[1:].split(" ")[0]
      elif line.startswith(">") and chrom is None:
        chrom = line[1:].split(" ")[0]
      else:
        cseq += line.strip()
    if chrom in chr2info: #last chr seq
      for eline in chr2info[chrom]:
        txid, txname, eseq = einfo(eline, cseq)
        fout.write('>' + txid + '|' + txname + '\n' + eseq + '\n')
  fin.close()
  fout.close()

###############################################################################
# Generate BLAST+ database
###############################################################################
def generate_BLASTdb(blastdb, txSeq_fn):
  print 'Generating BLAST+ database ...'
  exe = get_script_path() + "/tools/ncbi-blast/bin/makeblastdb "
  cmd = exe + ' -in ' + txSeq_fn + ' -dbtype nucl ' + ' -out ' + blastdb
  devnull = open(os.devnull, 'w')
  call([cmd], shell=True)#, stdout=devnull, stderr=devnull)
