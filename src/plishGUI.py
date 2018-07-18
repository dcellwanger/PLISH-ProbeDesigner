#!/usr/bin/python
###############################################################################
# Graphical user interface
# written by Daniel C. Ellwanger <dcellwanger.dev@gmail.com>
# last modified Jul 3 2018
###############################################################################
import os, tkMessageBox
from plishMain import main
from plishUtils import get_script_path
from Tkinter import Label, Listbox, Button, Entry, Text, mainloop, DISABLED, NORMAL, END, W, N, PhotoImage

###############################################################################
# Container
###############################################################################
class GUI:
  def _run(self):
    self.progressTxt.config(state=NORMAL)
    self.progressTxt.delete('1.0', END)
    self.progressTxt.update()
    self.progressTxt.config(state=DISABLED)
    inputId = self.txEntry.get()
    item = map(int, self.dbLbox.curselection())
    db = self.dbids[item[0]]
    print db
    self.runBtn.config(state=DISABLED)
    main(inputId, db, debug=self.debug, txt=self.progressTxt)
    self.runBtn.config(state=NORMAL)
    
  def _quitGUI(self):
    rpath = self.progressTxt.get('8.0','end-1c')
    if rpath.startswith('Your results'):
      tkMessageBox.showinfo("Quit", self.progressTxt.get('8.0','end-1c'))
    self.master.destroy()
    
  def __init__(self, master, dbnames, dbids, debug):
    self.dbids = dbids
    self.debug = debug
    
    self.master = master
    master.title('Plish Probe Designer')
    
    self.logoImg = PhotoImage(file=get_script_path() + '/img/plishLogo.gif')
    self.logoLbl = Label(master, image=self.logoImg)
    self.logoLbl.grid(row=0, columnspan=2)
    self.logoLbl.img = self.logoImg
    
    self.dbLbl = Label(master, text='Database')
    self.dbLbl.grid(row=1, sticky=W+N)
    self.dbLbox = Listbox(master, width=40)
    self.dbLbox.configure(exportselection=False)
    
    for i in range(len(dbnames)):
      self.dbLbox.insert(i, dbnames[i])
    self.dbLbox.select_set(0)
    self.dbLbox.grid(row=1, column=1)
    
    self.txLbl = Label(master, text='Transcript ID')
    self.txLbl.grid(row=2, sticky=W+N)
    self.txEntry = Entry(master, width=40)
    self.txEntry.grid(row=2, column=1)
  
    self.runBtn = Button(master, text='Run', command=self._run, width=30)
    self.runBtn.grid(row=3, column=1)
  
    self.progressLbl = Label(master, text='Progress')
    self.progressLbl.grid(row=4, sticky=W+N)
    self.progressTxt = Text(bg="#263238", fg="#ffffff", state=DISABLED, width=51, height=12)
    self.progressTxt.grid(row=4, column=1)
  
    self.quitBtn = Button(master, text='Quit', command=self._quitGUI, width=30)
    self.quitBtn.grid(row=5, column=1)
    
    self.versionLbl = Label(master, text='Version 0.2.0, 2018.\n' + \
                      '(c) Heller lab, Stanford School of Medicine\n' + \
                      'Daniel C. Ellwanger <dcellwanger.dev@gmail.com>')
    self.versionLbl.grid(row=6, column=1)
    
    #Attributes
    dbids = debug = None
    master = None
    logoImg = None
    logoLbl = dbLbl = txLbl = progressLbl = versionLbl = None
    dbLbox = None
    TxEntry = None 
    runBtn = quitBtn = None
    progressTxt = None
