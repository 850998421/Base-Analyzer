"""
Tools Name : Base Analyzer Tool
Original Developer: Afsana Rahman Snigdha
Original Development Date : 01-01-2012
Version : 1.0
Purpose : Draw the Both Rich Region window and import BothRichRegion file to perform its functional part
"""

#!/usr/bin/python -d
import pdb
import os 
import sys
import fileinput
import math
import re
import csv
import webbrowser
from PyQt4 import QtCore, QtGui
from FunctionalFiles.BothRichRegion import BothRichRegion


class BothRichRegionUI(QtGui.QDialog):

	def __init__(self, parent=None):
		QtGui.QWidget.__init__(self, parent)
		
#create a label		
		self.label_input_file=QtGui.QLabel(self)
		self.label_input_file.setText('Select a file of fasta sequence')
		self.label_input_file.setGeometry(QtCore.QRect(35,40,200,24))
		self.label_input_file.move(35, 40)
		self.label_input_file.show()

#create a line edit			
		self.input_file_browse = QtGui.QLineEdit(self)
		self.input_file_browse.move(35, 70)
		self.input_file_browse.setGeometry(QtCore.QRect(35, 70,335,28))
		self.input_file_browse.setEnabled(False)
		self.input_file_browse.show()

#create a button		
		self.input_file_browse_button = QtGui.QPushButton('Browse', self)
		self.input_file_browse_button.move(270, 100)
		self.input_file_browse_button.setGeometry(QtCore.QRect(270, 100,100,28))
		self.input_file_browse_button.show()
		self.input_file_browse_button.clicked.connect(self.inputFile)

#create a label
		self.label_percentage_at=QtGui.QLabel(self)
		self.label_percentage_at.setText('Considerable AT %')
		self.label_percentage_at.setGeometry(QtCore.QRect(35, 130,200,24))
		self.label_percentage_at.move(35, 130)
		self.label_percentage_at.show()

#create a combo box
		self.combo_percentage_at=QtGui.QComboBox(self)
		self.combo_percentage_at.insertItem(0,'--select--')
		self.combo_percentage_at.insertItem(1,'60')
		self.combo_percentage_at.insertItem(2,'70')
		self.combo_percentage_at.insertItem(3,'80')
		self.combo_percentage_at.insertItem(4,'90')
		self.combo_percentage_at.setGeometry(QtCore.QRect(170, 130,80,24))
		self.combo_percentage_at.move(170, 130)
		self.combo_percentage_at.show()
		self.combo_percentage_at.activated[str].connect(self.atPercetageCombo) 

#create a label				
		self.label_output_file_at=QtGui.QLabel(self)
		self.label_output_file_at.setText('AT Rich Output File')
		self.label_output_file_at.setGeometry(QtCore.QRect(35, 160,200,24))
		self.label_output_file_at.move(35, 160)
		self.label_output_file_at.show()

#create a line edit			
		self.file_browse_at = QtGui.QLineEdit(self)
		self.file_browse_at.move(35, 190)
		self.file_browse_at.setGeometry(QtCore.QRect(35,190,335,28))
		self.file_browse_at.setEnabled(False)
		self.file_browse_at.show()

#create a button		
		self.browse_button_at = QtGui.QPushButton('Browse', self)
		self.browse_button_at.move(270, 220)
		self.browse_button_at.setGeometry(QtCore.QRect(270, 220,100,28))
		self.browse_button_at.show()
		self.browse_button_at.clicked.connect(self.atOutputFile)

#create a label
		self.label_percentage_gc=QtGui.QLabel(self)
		self.label_percentage_gc.setText('Considerable GC %')
		self.label_percentage_gc.setGeometry(QtCore.QRect(35, 250,200,24))
		self.label_percentage_gc.move(35, 250)
		self.label_percentage_gc.show()
		
#create a combo box
		self.combo_percentage_gc=QtGui.QComboBox(self)
		self.combo_percentage_gc.insertItem(0,'--select--')
		self.combo_percentage_gc.insertItem(1,'60')
		self.combo_percentage_gc.insertItem(2,'70')
		self.combo_percentage_gc.insertItem(3,'80')
		self.combo_percentage_gc.insertItem(4,'90')
		self.combo_percentage_gc.setGeometry(QtCore.QRect(170, 250,80,24))
		self.combo_percentage_gc.move(170, 250)
		self.combo_percentage_gc.show()
		self.combo_percentage_gc.activated[str].connect(self.gcPercetageCombo) 

#create a label			
		self.label_output_file_gc=QtGui.QLabel(self)
		self.label_output_file_gc.setText('GC Rich Output File')
		self.label_output_file_gc.setGeometry(QtCore.QRect(35, 280,200,24))
		self.label_output_file_gc.move(35, 280)
		self.label_output_file_gc.show()
		
#create a line edit			
		self.file_browse_gc = QtGui.QLineEdit(self)
		self.file_browse_gc.move(35, 310)
		self.file_browse_gc.setGeometry(QtCore.QRect(35,310,335,28))
		self.file_browse_gc.setEnabled(False)
		self.file_browse_gc.show()

#create a button		
		self.browse_button_gc = QtGui.QPushButton('Browse', self)
		self.browse_button_gc.move(270, 340)
		self.browse_button_gc.setGeometry(QtCore.QRect(270, 340,100,28))
		self.browse_button_gc.show()
		self.browse_button_gc.clicked.connect(self.gcOutputFile)
		
#create a button				
		self.save_button = QtGui.QPushButton('Ok', self)
		self.save_button.move(150, 400)
		self.save_button.setGeometry(QtCore.QRect(150, 400,100,28))
		self.save_button.show()
		self.save_button.clicked.connect(self.saveFile)

#create a button			
		self.close_button = QtGui.QPushButton('Close', self)
		self.close_button.move(270, 400)
		self.close_button.setGeometry(QtCore.QRect(270,400,100,28))
		self.close_button.show()
		self.close_button.clicked.connect(self.close)

#create a window			
		self.setFixedSize(400,450)
		self.setGeometry(450, 220,160, 150)
		self.setWindowTitle('Both Rich Region')
		self.setModal(True)
		
		self.sourcefile = ""
		self.updated_file_name = ""
		self.percentage_of_at="--select--"
		self.percentage_of_gc="--select--"

#function to get the value at percentage combobox	
	def atPercetageCombo(self,text):
		self.percentage_of_at = text
		#print self.percentage_of_at	

#function to get the value gc percentage combobox
	def gcPercetageCombo(self, text):
		self.percentage_of_gc = text
		#print self.percentage_of_gc	

#function to perform writing the output file			
	def saveFile(self):
		if str(self.sourcefile)== "":
			self.msgBox=QtGui.QMessageBox()
			self.msgBox.setWindowTitle("!!!!!!!!!Warnings!!!!!!!!!")
			self.msgBox.setText("Please select a input file")
			self.msgBox.move(500, 330)
			self.msgBox.exec_()
			
		if str(self.updated_file_name_at)== "":
			self.msgBox=QtGui.QMessageBox()
			self.msgBox.setWindowTitle("!!!!!!!!!Warnings!!!!!!!!!")
			self.msgBox.setText("Please give a output file AT Rich")
			self.msgBox.move(500, 330)
			self.msgBox.exec_()	
			
		if str(self.percentage_of_at)== "--select--":
			self.msgBox=QtGui.QMessageBox()
			self.msgBox.setWindowTitle("!!!!!!!!!Warnings!!!!!!!!!")
			self.msgBox.setText("Please select GC percentage")
			self.msgBox.move(500, 330)
			self.msgBox.exec_()
			
		if str(self.percentage_of_gc)== "--select--":
			self.msgBox=QtGui.QMessageBox()
			self.msgBox.setWindowTitle("!!!!!!!!!Warnings!!!!!!!!!")
			self.msgBox.setText("Please select GC percentage")
			self.msgBox.move(500, 330)
			self.msgBox.exec_()
												
		if str(self.updated_file_name_gc)== "--select--":
			self.msgBox=QtGui.QMessageBox()
			self.msgBox.setWindowTitle("!!!!!!!!!Warnings!!!!!!!!!")
			self.msgBox.setText("Please give a output file GC Rich")
			self.msgBox.move(500, 330)
			self.msgBox.exec_()
		
		if str(self.sourcefile)!= ""  and str(self.updated_file_name_at)!= "" and str(self.updated_file_name_gc)!= "" and str(self.percentage_of_at)!="--select--" and str(self.percentage_of_gc)!="--select--":
			self.W = BothRichRegion()
			atrich = self.W.ATRichRegion(self.sourcefile,self.updated_file_name_at,self.percentage_of_at);
			gcrich = self.W.GCRichRegion(self.sourcefile,self.updated_file_name_gc,self.percentage_of_gc);
			
			if str(atrich) == "yes" and str(gcrich) == "yes":
				self.msgBox=QtGui.QMessageBox()
				self.msgBox.setWindowTitle("             Success Message             ")
				self.msgBox.setText("Result file save successfully.")
				self.msgBox.move(500, 330)
				ret = self.msgBox.exec_()
				if ret :
					self.close()
					self.input_file_browse.setText("")
					self.file_browse_gc.setText("")
					self.file_browse_at.setText("")
					webbrowser.open(self.updated_file_name_at)
					webbrowser.open(self.updated_file_name_gc)
					
			elif str(atrich) != "yes" or str(gcrich) != "yes":
				if str(atrich) != "yes":					
					self.msgBox=QtGui.QMessageBox()
					self.msgBox.setWindowTitle("!!!!!!!!!Warnings!!!!!!!!!")
					self.msgBox.setText(str(atrich))
					self.msgBox.move(500, 330)
					self.msgBox.exec_()
					self.input_file_browse.setText("")
					self.file_browse_gc.setText("")
					self.file_browse_at.setText("")
					
				if str(gcrich) != "yes":					
					self.msgBox=QtGui.QMessageBox()
					self.msgBox.setWindowTitle("!!!!!!!!!Warnings!!!!!!!!!")
					self.msgBox.setText(str(gcrich))
					self.msgBox.move(500, 330)
					self.msgBox.exec_()
					self.input_file_browse.setText("")
					self.file_browse_gc.setText("")
					self.file_browse_at.setText("")

#function to take the output file name for gc rich from user									
	def gcOutputFile(self):								
		output_file_name = QtGui.QFileDialog.getSaveFileName(self, 'Save File', '/home',"HTML Files (*.html)")
		if(int(str(output_file_name).find(".")) != -1):
			self.updated_file_name_gc=output_file_name[0:int(str(output_file_name).find("."))] +".html"
			self.file_browse_gc.setText(str(self.updated_file_name_gc));
		elif str(output_file_name)!= "":
			self.updated_file_name_gc=output_file_name +".html"
			self.file_browse_gc.setText(str(self.updated_file_name_gc)) 			

#function to take the output file name for at rich from user
	def atOutputFile(self):								
		output_file_name = QtGui.QFileDialog.getSaveFileName(self, 'Save File', '/home',"HTML Files (*.html)")
		if(int(str(output_file_name).find(".")) != -1):
			self.updated_file_name_at=output_file_name[0:int(str(output_file_name).find("."))] +".html"
			self.file_browse_at.setText(str(self.updated_file_name_at));
		else:
			self.updated_file_name_at=output_file_name +".html"
			self.file_browse_at.setText(str(self.updated_file_name_at))
			
#function to take the input file name from user					
	def inputFile(self):
		self.sourcefile = QtGui.QFileDialog.getOpenFileName(self, 'Open file', '/home',"Text Files (*.txt)")
		file_size = os.path.getsize(self.sourcefile)
		if file_size > 1073741824:
			self.msgBox=QtGui.QMessageBox()
			self.msgBox.setWindowTitle("!!!!!!!!!Warnings!!!!!!!!!")
			self.msgBox.setText("Input file is too Large.")
			self.msgBox.move(500, 330)
			ret = self.msgBox.exec_()
		elif file_size < 1073741824:	
			self.input_file_browse.setText(str(self.sourcefile))
		
