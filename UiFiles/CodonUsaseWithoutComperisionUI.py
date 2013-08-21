"""
Tools Name : Base Analyzer Tool
Original Developer: Afsana Rahman Snigdha
Original Development Date : 01-01-2012
Version : 1.0
Purpose : Draw the Codon Usage Table window and import CodonUsaseWithoutComperision file to perform its functional part
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
from FunctionalFiles.CodonUsaseWithoutComperision import codonUsaseTableWithoutComperision


class codonUsaseWithoutComperisionUI(QtGui.QDialog):

	def __init__(self, parent=None):
		QtGui.QWidget.__init__(self, parent)

#create a lebel		
		self.label_input_file=QtGui.QLabel(self)
		self.label_input_file.setText('Select a file of fasta sequence')
		self.label_input_file.setGeometry(QtCore.QRect(35,40,200,24))
		self.label_input_file.move(35, 40)
		self.label_input_file.show()

#create a edit line			
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
		self.label_output_file_format=QtGui.QLabel(self)
		self.label_output_file_format.setText('Select a output file format')
		self.label_output_file_format.setGeometry(QtCore.QRect(35, 130,300,24))
		self.label_output_file_format.move(35, 130)
		self.label_output_file_format.show()

#create a radio button			
		self.csv_radio_button = QtGui.QRadioButton('csv', self)
		self.csv_radio_button.setCheckable(True)
		self.csv_radio_button.move(30, 150)
		self.csv_radio_button.show()
		self.csv_radio_button.clicked.connect(self.csvFile)

#create a radio button			
		self.txt_radio_button = QtGui.QRadioButton('txt', self)
		self.txt_radio_button.setCheckable(True)
		self.txt_radio_button.move(90, 150)
		self.txt_radio_button.show()
		self.txt_radio_button.clicked.connect(self.txtFile)

#create a label		
		self.label_output_file=QtGui.QLabel(self)
		self.label_output_file.setText('Save File To ')
		self.label_output_file.setGeometry(QtCore.QRect(35, 180,200,24))
		self.label_output_file.move(35, 180)
		self.label_output_file.show()

#create a line edit			
		self.file_browse = QtGui.QLineEdit(self)
		self.file_browse.move(35, 210)
		self.file_browse.setGeometry(QtCore.QRect(35,210,335,28))
		self.file_browse.setEnabled(False)
		self.file_browse.show()
		
#create a button		
		self.browse_button = QtGui.QPushButton('Browse', self)
		self.browse_button.move(270, 240)
		self.browse_button.setGeometry(QtCore.QRect(270, 240,100,28))
		self.browse_button.show()
		self.browse_button.clicked.connect(self.outputFile)

#create a button		
		self.save_button = QtGui.QPushButton('Ok', self)
		self.save_button.move(150, 330)
		self.save_button.setGeometry(QtCore.QRect(150, 330,100,28))
		self.save_button.show()
		self.save_button.clicked.connect(self.saveFile)

#create a button		
		self.close_button = QtGui.QPushButton('Close', self)
		self.close_button.move(270, 330)
		self.close_button.setGeometry(QtCore.QRect(270,330,100,28))
		self.close_button.show()
		self.close_button.clicked.connect(self.close)

#create a window		
		self.setFixedSize(400,400)
		self.setGeometry(450, 250,160, 150)
		self.setWindowTitle('Codon Usage Table (Without Comparison)')
		self.setModal(True)
		
		self.sourcefile = ""
		self.updated_file_name = ""
		self.txt_radiobutton = "0"
		self.csv_radiobutton="0"
		
#function to get the value of txt radio button			
	def txtFile(self):
		self.txt_radiobutton="1"
		self.csv_radiobutton="0"

#function to get the value of csv radio button			
	def csvFile(self):
		self.csv_radiobutton="1"
		self.txt_radiobutton="0"		  	

#function to perform writing the output file			
	def saveFile(self):
		if str(self.sourcefile)== "":
			self.msgBox=QtGui.QMessageBox()
			self.msgBox.setWindowTitle("!!!!!!!!!Warnings!!!!!!!!!")
			self.msgBox.setText("Please select a input file")
			self.msgBox.move(500, 330)
			self.msgBox.exec_()
			
		if str(self.updated_file_name)== "":
			self.msgBox=QtGui.QMessageBox()
			self.msgBox.setWindowTitle("!!!!!!!!!Warnings!!!!!!!!!")
			self.msgBox.setText("Please give a output file")
			self.msgBox.move(500, 330)
			self.msgBox.exec_()				
		
		if str(self.sourcefile)!= ""  and str(self.updated_file_name)!= "":
			self.W = codonUsaseTableWithoutComperision()
			codonusase = self.W.codonUsaseTableWithoutComperision(self.sourcefile,self.updated_file_name);
			if str(codonusase) == "yes":
				self.msgBox=QtGui.QMessageBox()
				self.msgBox.setWindowTitle("             Success Message             ")
				self.msgBox.setText("Result file save successfully.")
				self.msgBox.move(500, 330)
				ret = self.msgBox.exec_()
				if ret :
					self.close()
					self.input_file_browse.setText("")
					self.file_browse.setText("")
					self.csv_radio_button.setChecked(False)
					self.txt_radio_button.setChecked(False)
					file_size = os.path.getsize(self.updated_file_name)
					if self.csv_radiobutton == "1":
						if file_size > 40960:
							self.msgBox=QtGui.QMessageBox()
							self.msgBox.setWindowTitle("!!!!!!!!!Warnings!!!!!!!!!")
							self.msgBox.setText("Result file is too large to open. File is saved in " + self.updated_file_name + " location")
							self.msgBox.move(500, 330)
							ret = self.msgBox.exec_()
						elif file_size <= 40960:								
							webbrowser.open(self.updated_file_name)
					elif self.txt_radiobutton =="1":
						if file_size > 40960:
							self.msgBox=QtGui.QMessageBox()
							self.msgBox.setWindowTitle("!!!!!!!!!Warnings!!!!!!!!!")
							self.msgBox.setText("Result file is too large to open. File is saved in " + self.updated_file_name + " location")
							self.msgBox.move(500, 330)
							ret = self.msgBox.exec_()
						elif file_size <= 40960:								
							webbrowser.open(self.updated_file_name)
			elif str(codonusase) != "yes":
				self.msgBox=QtGui.QMessageBox()
				self.msgBox.setWindowTitle("!!!!!!!!!Warnings!!!!!!!!!")
				self.msgBox.setText(str(codonusase))
				self.msgBox.move(500, 330)
				self.msgBox.exec_()

#function to take the output file name from user				
	def outputFile(self):
		if str(self.csv_radiobutton)== "0" and str(self.txt_radiobutton)== "0":
			self.msgBox=QtGui.QMessageBox()
			self.msgBox.setWindowTitle("!!!!!!!!!Warnings!!!!!!!!!")
			self.msgBox.setText("Please select a output file format")
			self.msgBox.move(500, 330)
			self.msgBox.exec_()
			
		elif str(self.csv_radiobutton)== "1":
			output_file_name = QtGui.QFileDialog.getSaveFileName(self, 'Save File', '/home',"CSV Files (*.csv)")
			if(int(str(output_file_name).find(".")) != -1):
				self.updated_file_name=output_file_name[0:int(str(output_file_name).find("."))] +".csv"
				self.file_browse.setText(str(self.updated_file_name));
			elif(str(output_file_name) != ""):
				self.updated_file_name=output_file_name  +".csv"
				self.file_browse.setText(str(self.updated_file_name))	
			
		elif str(self.txt_radiobutton)== "1":
			output_file_name = QtGui.QFileDialog.getSaveFileName(self, 'Save File', '/home',"Text Files (*.txt)")
			if(int(str(output_file_name).find(".")) != -1):
				self.updated_file_name=output_file_name[0:int(str(output_file_name).find("."))] +".txt"
				self.file_browse.setText(str(self.updated_file_name));
			elif(str(output_file_name) != ""):
				self.updated_file_name=output_file_name +".txt"
				self.file_browse.setText(str(self.updated_file_name)) 			

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
		
