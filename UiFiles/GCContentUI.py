"""
Tools Name : Base Analyzer Tool
Original Developer: Afsana Rahman Snigdha
Original Development Date : 01-01-2012
Version : 1.0
Purpose : Draw the GC Count window and import GCContent file to perform its functional part
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
from FunctionalFiles.GCContent import gcContent


class gcContentUI(QtGui.QDialog):

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
		self.label_output_file_format=QtGui.QLabel(self)
		self.label_output_file_format.setText('Select an option of bellow')
		self.label_output_file_format.setGeometry(QtCore.QRect(35, 130,300,24))
		self.label_output_file_format.move(35, 130)
		self.label_output_file_format.show()
		
#create a check box			
		self.first_positon_radio_button = QtGui.QCheckBox('First Position', self)
		self.first_positon_radio_button.setCheckable(True)
		self.first_positon_radio_button.move(30, 160)
		self.first_positon_radio_button.show()
		self.first_positon_radio_button.clicked.connect(self.firstPosition)

#create a check box			
		self.second_positon_radio_button = QtGui.QCheckBox('Second Position', self)
		self.second_positon_radio_button.setCheckable(True)
		self.second_positon_radio_button.move(160, 160)
		self.second_positon_radio_button.show()
		self.second_positon_radio_button.clicked.connect(self.secondPosition)

#create a check box			
		self.third_position_radio_button = QtGui.QCheckBox('Third Position', self)
		self.third_position_radio_button.setCheckable(True)
		self.third_position_radio_button.move(30, 190)
		self.third_position_radio_button.show()
		self.third_position_radio_button.clicked.connect(self.thirdPosition)

#create a check box			
		self.full_genome_radio_button = QtGui.QCheckBox('Full Genome', self)
		self.full_genome_radio_button.setCheckable(True)
		self.full_genome_radio_button.move(160, 190)
		self.full_genome_radio_button.show()
		self.full_genome_radio_button.clicked.connect(self.fullGenomePosition)

#create a label			
		self.label_output_file=QtGui.QLabel(self)
		self.label_output_file.setText('Save File To ')
		self.label_output_file.setGeometry(QtCore.QRect(35, 220,200,24))
		self.label_output_file.move(35, 220)
		self.label_output_file.show()

#create a line edit		
		self.file_browse = QtGui.QLineEdit(self)
		self.file_browse.move(35, 250)
		self.file_browse.setGeometry(QtCore.QRect(35,250,335,28))
		self.file_browse.setEnabled(False)
		self.file_browse.show()
		
#create a button		
		self.browse_button = QtGui.QPushButton('Browse', self)
		self.browse_button.move(270, 280)
		self.browse_button.setGeometry(QtCore.QRect(270, 280,100,28))
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
		self.setWindowTitle('GC Count')
		self.setModal(True)
		
		self.sourcefile = ""
		self.updated_file_name = ""
		self.first_Position_checkbox="0"
		self.second_Position_checkbox="0"
		self.third_Position_checkbox="0"
		self.full_genome_Position_checkbox="0"

#function to get the value of first position checkbox		
	def firstPosition(self):
		self.first_Position_checkbox="1"
		
#function to get the value of second position checkbox		
	def secondPosition(self):
		self.second_Position_checkbox="1"

#function to get the value of third position checkbox			
	def thirdPosition(self):
		self.third_Position_checkbox="1"
		
#function to get the value of full genome checkbox					  	
	def fullGenomePosition(self):
		self.full_genome_Position_checkbox="1"

#function to perform writing the output file			
	def saveFile(self):		
		if str(self.sourcefile)== "":
			self.msgBox=QtGui.QMessageBox()
			self.msgBox.setWindowTitle("!!!!!!!!!Warnings!!!!!!!!!")
			self.msgBox.setText("Please select a input file")
			self.msgBox.move(500, 330)
			self.msgBox.exec_()

		elif str(self.first_Position_checkbox)== "0" and str(self.second_Position_checkbox)== "0" and str(self.third_Position_checkbox)== "0" and str(self.full_genome_Position_checkbox)== "0":
			self.msgBox=QtGui.QMessageBox()
			self.msgBox.setWindowTitle("!!!!!!!!!Warnings!!!!!!!!!")
			self.msgBox.setText("Please select an option")
			self.msgBox.move(500, 330)
			self.msgBox.exec_()			
		
		elif str(self.updated_file_name)== "":
			self.msgBox=QtGui.QMessageBox()
			self.msgBox.setWindowTitle("!!!!!!!!!Warnings!!!!!!!!!")
			self.msgBox.setText("Please give a output file")
			self.msgBox.move(500, 330)
			self.msgBox.exec_()	
			
		if str(self.sourcefile)!= ""  and str(self.updated_file_name)!= "":
			self.W = gcContent()
			codonusase = self.W.gcContent(self.sourcefile,self.updated_file_name,self.first_Position_checkbox,self.second_Position_checkbox,self.third_Position_checkbox,self.full_genome_Position_checkbox);
			if str(codonusase) == "yes":
				self.msgBox=QtGui.QMessageBox()
				self.msgBox.setWindowTitle("             Success Message             ")
				self.msgBox.setText("Result file save successfully.")
				self.msgBox.move(500, 330)
				ret = self.msgBox.exec_()
				if ret :
					self.close()
					file_size = os.path.getsize(self.updated_file_name)
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
		if str(self.first_Position_checkbox)== "0" and str(self.second_Position_checkbox)== "0" and str(self.third_Position_checkbox)== "0" and str(self.full_genome_Position_checkbox)== "0":
			self.msgBox=QtGui.QMessageBox()
			self.msgBox.setWindowTitle("!!!!!!!!!Warnings!!!!!!!!!")
			self.msgBox.setText("Please select an option")
			self.msgBox.move(500, 330)
			self.msgBox.exec_()	
			
		else:
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
		
