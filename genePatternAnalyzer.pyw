#!/usr/bin/env python

"""
Tools Name : Base Analyzer Tool
Original Developer: Afsana Rahman Snigdha
Original Development Date : 01-01-2012
Version : 1.0
Purpose : Draw the main window with manubar and import all the UI files
			that are used in this tool
"""
import pdb
import sys
import fileinput
import math
import re
from PyQt4 import QtGui
from PyQt4 import QtCore
from UiFiles.CodonUsaseWithComperisionUI import codonUsaseWithComperisionUI
from UiFiles.CodonUsaseWithoutComperisionUI import codonUsaseWithoutComperisionUI
from UiFiles.GCContentUI import gcContentUI
from UiFiles.ATRichRegionUI import ATRichRegionUI
from UiFiles.GCRichRegionUI import GCRichRegionUI
from UiFiles.BothRichRegionUI import BothRichRegionUI

class GenePatternAnalyzer(QtGui.QMainWindow):
	
		def __init__(self, parent=None):
			super(GenePatternAnalyzer, self).__init__()
			self.initUI()
						        
		def withComperisionCodonUsase(self):
			self.CUW.show()
			
		def withoutComperisionCodonUsase(self):
			self.CUWO.show()
			
		def gcContent(self):
			self.GC.show()		
				
		def atRichRegion(self):
			self.ATR.show()
		
		def bothRichRegion(self):
			self.ATGCR.show()
		
		def gcRichRegion(self):
			self.GCR.show()
						
		def initUI(self): 
			self.statusBar()
			withCompersion = QtGui.QAction('With Comparision', self)
			withCompersion.setShortcut('Ctrl+C')
			withCompersion.setStatusTip('With Comparision')
			self.CUW = codonUsaseWithComperisionUI()
			withCompersion.triggered.connect(self.withComperisionCodonUsase)
			
			withoutCompersion = QtGui.QAction('Without Comparison', self)
			withoutCompersion.setShortcut('Ctrl+W')
			withoutCompersion.setStatusTip('Without Comparision')
			self.CUWO = codonUsaseWithoutComperisionUI()
			withoutCompersion.triggered.connect(self.withoutComperisionCodonUsase)
			
			gcCounting = QtGui.QAction('GC Count', self)
			gcCounting.setShortcut('Ctrl+O')
			gcCounting.setStatusTip('GC Count')
			self.GC = gcContentUI()
			gcCounting.triggered.connect(self.gcContent)

			atRich = QtGui.QAction('AT Rich', self)
			atRich.setShortcut('Ctrl+A')
			atRich.setStatusTip('AT Rich')
			self.ATR = ATRichRegionUI()
			atRich.triggered.connect(self.atRichRegion)
			
			gcRich = QtGui.QAction('GC Rich', self)
			gcRich.setShortcut('Ctrl+G')
			gcRich.setStatusTip('GC Rich')
			self.GCR = GCRichRegionUI()
			gcRich.triggered.connect(self.gcRichRegion)        

			bothRich = QtGui.QAction('AT and GC Rich', self)
			bothRich.setShortcut('Ctrl+R')
			bothRich.setStatusTip('AT and GC Rich')
			self.ATGCR = BothRichRegionUI()
			bothRich.triggered.connect(self.bothRichRegion) 
		
			menubar = self.menuBar()
			CodonUsase = menubar.addMenu('Codon Usage')
			CodonUsase.addAction(withCompersion)
			CodonUsase.addAction(withoutCompersion)
			
			GCContent = menubar.addMenu('GC Count')
			GCContent.addAction(gcCounting)
			
			RichRegion = menubar.addMenu('AT/GC Rich Region')
			RichRegion.addAction(atRich)
			RichRegion.addAction(gcRich)
			RichRegion.addAction(bothRich)
			
			self.setWindowTitle('Gene Pattern Analyzer')
			self.setFixedSize(600,550)
			self.setGeometry(350,150,30,30)
			#self.setStyleSheet("background-image: url(/home/afsana/Pictures/Bio-Informatics/dna_rgb.gif);")
			#self.showMaximized()
			self.show()
         
def main():
    app = QtGui.QApplication(sys.argv)
    ex = GenePatternAnalyzer()
    sys.exit(app.exec_())


if __name__ == '__main__':
    main()
