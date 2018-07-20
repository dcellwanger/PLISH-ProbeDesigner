#!/usr/bin/python
###############################################################################
# Graphical user interface
# written by Daniel C. Ellwanger <dcellwanger.dev@gmail.com>
# last modified Jul 3 2018
###############################################################################
import os, tkMessageBox
from plishMain import main
from plishUtils import get_script_path, filter_probes, write_probesCSV, write_probesFNA
from Tkinter import Label, Listbox, Button, Entry, Text, StringVar, OptionMenu
from Tkinter import mainloop, PhotoImage, IntVar, BooleanVar, DoubleVar
from Tkinter import PanedWindow, LabelFrame, Checkbutton
from Tkinter import DISABLED, NORMAL, END, W, N, E, LEFT, VERTICAL

# Version and default values
_version = '0.3.0 (2018)'
_defaultGC = 45
_defaultMultiExon = True
_defaultMinTm = 55
_defaultMaxTm = 65
_defaultMinDimer = -10
_defaultMinFold = -10
_defaultMaxDuplex = -30

###############################################################################
# Container
###############################################################################
class GUI:
  def _run(self):
    self.saveBtn.config(state=DISABLED)
    self.progressTxt.config(state=NORMAL)
    self.progressTxt.delete('1.0', END)
    self.progressTxt.update()
    self.progressTxt.config(state=DISABLED)
    inputId = self.txEntry.get()
    item = map(int, self.dbLbox.curselection())
    db = self.dbids[item[0]]
    self.runBtn.config(state=DISABLED)
    self.inputId, self.inputName, self.hprobes = main(inputId, db, debug=self.debug, txt=self.progressTxt)
    self.runBtn.config(state=NORMAL)
    if self.hprobes is not None:
      self.saveBtn.config(state=NORMAL)
    
  def _quitGUI(self):
    #rpath = self.progressTxt.get('8.0','end-1c')
    #if rpath.startswith('Your results'):
    #  tkMessageBox.showinfo("Quit", self.progressTxt.get('8.0','end-1c'))
    self.master.destroy()
    
  def _save(self):
    hps = filter_probes(self.hprobes, self.spec.get(), self.mingc.get(), 
                        self.multiexon.get(), self.mintm.get(), self.maxtm.get(),
                        self.mindimer.get(), self.minfold.get(), self.maxduplex.get())
    result_csv = write_probesCSV(self.inputId, self.inputName, hps, self.progressTxt)
    result_fna = write_probesFNA(self.inputId, self.inputName, hps, self.progressTxt)
    tkMessageBox.showinfo('Result file', 
                          'Details on ' + str(len(hps)) + \
                          ' hybridization probe(s) were exported to ' + \
                          result_csv + "\n\n" + \
                          'Sequences of '+ str(len(hps)) + \
                          ' hybridization probe(s) were exported to ' + \
                          result_fna)
                          
  def __init__(self, master, dbnames, dbids, debug):
    self.dbids = dbids
    self.debug = debug
    
    self.master = master
    master.title('Plish Probe Designer')
    
    self.logoImg = PhotoImage(file=get_script_path() + '/img/plishLogo.gif')
    self.logoLbl = Label(master, image=self.logoImg)
    self.logoLbl.grid(row=0, columnspan=3)
    self.logoLbl.img = self.logoImg
    
    self.dbLbl = Label(master, text='Database')
    self.dbLbl.grid(row=1, sticky=W+N)
    self.dbLbox = Listbox(master, width=63, height=4)
    self.dbLbox.configure(exportselection=False)
    
    for i in range(len(dbnames)):
      self.dbLbox.insert(i, dbnames[i])
    self.dbLbox.select_set(0)
    self.dbLbox.grid(row=1, column=1, columnspan=2, sticky=N)
    
    self.txLbl = Label(master, text='Transcript ID')
    self.txLbl.grid(row=2, sticky=W+N)
    self.txEntry = Entry(master, width=39)
    self.txEntry.grid(row=2, column=1, sticky=W+N)
  
    self.runBtn = Button(master, text='Run', command=self._run, width=15)
    self.runBtn.grid(row=2, column=2)
  
    self.progressLbl = Label(master, text='Progress')
    self.progressLbl.grid(row=4, sticky=W+N)
    self.progressTxt = Text(bg="#263238", fg="#ffffff", state=DISABLED, width=51, height=16)
    self.progressTxt.grid(row=4, column=1)
  
    self.saveBtn = Button(master, text='Save', command=self._save, state=DISABLED, width=15)
    self.saveBtn.grid(row=5, column=2, sticky=N)
  
    self.quitBtn = Button(master, text='Quit', command=self._quitGUI, width=15)
    self.quitBtn.grid(row=6, column=2, sticky=N)
    
    self.aboutLF = LabelFrame(master, text='About', width=300)
    self.aboutLF.grid(row=5, column=0, rowspan=2, columnspan=2, sticky=N+W)
    self.versionLbl = Label(self.aboutLF, text='PLISH Probe Designer, Version ' + _version + '\n' + \
                      '(c) Heller lab, Stanford University School of Medicine\n' + \
                      '     Daniel C. Ellwanger <dcellwanger.dev@gmail.com>                       ', 
                      justify=LEFT)
    self.versionLbl.grid(row=0, column=0, sticky=N)
    
    # Filter
    self.filterLF = LabelFrame(master, text='Filter')
    self.filterLF.grid(row=4, column=2, rowspan=2, sticky=N+W)
    
    self.mingc = DoubleVar()
    self.mingc.set(_defaultGC)
    self.mingcLbl = Label(self.filterLF, text='Min. GC')
    self.mingcLbl.grid(row=0, column=0, sticky=N+W)
    self.mingcEntry = Entry(self.filterLF, width=5, text=self.mingc)
    self.mingcEntry.grid(row=0, column=1, sticky=N+W)
    self.mingcLbl2 = Label(self.filterLF, text='%')
    self.mingcLbl2.grid(row=0, column=2, sticky=N+W)

    self.spec = StringVar(master)
    self.spec.set("isoform")
    self.specLbl = Label(self.filterLF, text='Specificity')
    self.specLbl.grid(row=1, column=0, sticky=N+W)
    self.specOm = OptionMenu(self.filterLF, self.spec, "isoform", "gene", "none")
    self.specOm.grid(row=1, column=1, sticky=N+W, columnspan=2)
    
    self.mintm = DoubleVar()
    self.mintm.set(_defaultMinTm)
    self.mintmLbl = Label(self.filterLF, text='Min. Tm')
    self.mintmLbl.grid(row=2, column=0, sticky=N+W)
    self.mintmEntry = Entry(self.filterLF, width=5, text=self.mintm)
    self.mintmEntry.grid(row=2, column=1, sticky=N+W)
    self.mintmLbl2 = Label(self.filterLF, text=u'\N{DEGREE SIGN}'+ 'C')
    self.mintmLbl2.grid(row=2, column=2, sticky=N+W)
    
    self.maxtm = DoubleVar()
    self.maxtm.set(_defaultMaxTm)
    self.maxtmLbl = Label(self.filterLF, text='Max. Tm')
    self.maxtmLbl.grid(row=3, column=0, sticky=N+W)
    self.maxtmEntry = Entry(self.filterLF, width=5, text=self.maxtm)
    self.maxtmEntry.grid(row=3, column=1, sticky=N+W)
    self.maxtmLbl2 = Label(self.filterLF, text=u'\N{DEGREE SIGN}'+ 'C')
    self.maxtmLbl2.grid(row=3, column=2, sticky=N+W)
    
    self.minfold = DoubleVar()
    self.minfold.set(_defaultMinFold)
    self.minfoldLbl = Label(self.filterLF, text='Min. Fold')
    self.minfoldLbl.grid(row=4, column=0, sticky=N+W)
    self.minfoldEntry = Entry(self.filterLF, width=5, text=self.minfold)
    self.minfoldEntry.grid(row=4, column=1, sticky=N+W)
    self.minfoldLbl2 = Label(self.filterLF, text='kcal/mol')
    self.minfoldLbl2.grid(row=4, column=2, sticky=N+W)
    
    self.mindimer = DoubleVar()
    self.mindimer.set(_defaultMinDimer)
    self.mindimerLbl = Label(self.filterLF, text='Min. Dimer')
    self.mindimerLbl.grid(row=5, column=0, sticky=N+W)
    self.mindimerEntry = Entry(self.filterLF, width=5, text=self.mindimer)
    self.mindimerEntry.grid(row=5, column=1, sticky=N+W)
    self.mindimerLbl2 = Label(self.filterLF, text='kcal/mol')
    self.mindimerLbl2.grid(row=5, column=2, sticky=N+W)
  
    self.maxduplex = DoubleVar()
    self.maxduplex.set(_defaultMaxDuplex)
    self.maxduplexLbl = Label(self.filterLF, text='Max. Duplex')
    self.maxduplexLbl.grid(row=6, column=0, sticky=N+W)
    self.maxduplexEntry = Entry(self.filterLF, width=5, text=self.maxduplex)
    self.maxduplexEntry.grid(row=6, column=1, sticky=N+W)
    self.maxduplexLbl2 = Label(self.filterLF, text='kcal/mol')
    self.maxduplexLbl2.grid(row=6, column=2, sticky=N+W)    

    self.multiexon = BooleanVar()
    self.multiexon.set(_defaultMultiExon)
    self.multiexonCb = Checkbutton(self.filterLF, text='Multi-exon', 
                                   variable=self.multiexon, 
                                   onvalue=True, offvalue=False)
    self.multiexonCb.grid(row=7, column=0, sticky=N+W)
