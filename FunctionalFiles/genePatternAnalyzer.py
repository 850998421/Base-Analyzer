#!/usr/bin/env python

import pdb
import sys
import fileinput
import math
import re
from PyQt4 import QtGui
from PyQt4 import QtCore
from CodonUsaseWithComperisionUI import codonUsaseWithComperisionUI
from CodonUsaseWithoutComperisionUI import codonUsaseWithoutComperisionUI
from GCContentUI import gcContentUI
from ATRichRegionUI import ATRichRegionUI
from GCRichRegionUI import GCRichRegionUI
from BothRichRegionUI import BothRichRegionUI

class GenePatternAnalyzer(QtGui.QMainWindow):
	
		def __init__(self, parent=None):
			super(GenePatternAnalyzer, self).__init__()
			self.initUI()      
		def bothRichRegion(self):
			self.ATGCR.show()
								
		def withComperisionCodonUsase(self):
			self.CUW.show()
			
		def withoutComperisionCodonUsase(self):
			self.CUWO.show()
							
		def gcContent(self):
			self.GC.show()					
						
		def atRichRegion(self):
			self.ATR.show()

		def gcRichRegion(self):
			self.GCR.show()		
			
		def initUI(self): 
			self.statusBar()
			withCompersion = QtGui.QAction('With Comparision', self)
			withCompersion.setShortcut('Ctrl+C')
			withCompersion.setStatusTip('With Comparision')
			self.CUW = codonUsaseWithComperisionUI()
			withCompersion.triggered.connect(self.withComperisionCodonUsase)
			
			withoutCompersion = QtGui.QAction('Without Comparision', self)
			withoutCompersion.setShortcut('Ctrl+w')
			withoutCompersion.setStatusTip('Without Comparision')
			self.CUWO = codonUsaseWithoutComperisionUI()
			withoutCompersion.triggered.connect(self.withoutComperisionCodonUsase)
			
			gcCounting = QtGui.QAction('GC Content', self)
			gcCounting.setShortcut('Ctrl+g')
			gcCounting.setStatusTip('GC Content')
			self.GC = gcContentUI()
			gcCounting.triggered.connect(self.gcContent)

			atRich = QtGui.QAction('AT Rich', self)
			atRich.setShortcut('Ctrl+a')
			atRich.setStatusTip('AT Rich')
			self.ATR = ATRichRegionUI()
			atRich.triggered.connect(self.atRichRegion)
			
			gcRich = QtGui.QAction('GC Rich', self)
			gcRich.setShortcut('Ctrl+c')
			gcRich.setStatusTip('GC Rich')
			self.GCR = GCRichRegionUI()
			gcRich.triggered.connect(self.gcRichRegion)        

			bothRich = QtGui.QAction('AT & GC Rich', self)
			bothRich.setShortcut('Ctrl+t')
			bothRich.setStatusTip('AT & GC Rich')
			self.ATGCR = BothRichRegionUI()
			bothRich.triggered.connect(self.bothRichRegion) 
		
			menubar = self.menuBar()
			CodonUsase = menubar.addMenu('Codon Usase')
			CodonUsase.addAction(withCompersion)
			CodonUsase.addAction(withoutCompersion)
			
			GCContent = menubar.addMenu('GC Counting')
			GCContent.addAction(gcCounting)
			
			RichRegion = menubar.addMenu('Rich Region')
			RichRegion.addAction(atRich)
			RichRegion.addAction(gcRich)
			RichRegion.addAction(bothRich)
			
			self.setWindowTitle('Gene Pattern Analyzer')
			self.setFixedSize(600,500)
			self.setGeometry(400,200,90,80)
			#self.setStyleSheet("background-image: url(/home/afsana/Pictures/Bio-Informatics/dna_rgb.gif);")
			#self.showMaximized()
			self.show()
         
def main():
    app = QtGui.QApplication(sys.argv)
    ex = GenePatternAnalyzer()
    sys.exit(app.exec_())


if __name__ == '__main__':
    main()
