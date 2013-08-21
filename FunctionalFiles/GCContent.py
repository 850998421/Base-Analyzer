"""
Tools Name : Base analyzer Tool
Original Developer: afsana Rahman Snigdha
Original Development Date : 01-01-2012
Version : 1.0
Purpose : 
		Input File: a text file containing one or more fasta sequence.
		Output File : a text file with the GC Count result with the percentage.
		Methods : Count GC in First, second, Third position in a codon and count GC in whole sequence.

"""

import pdb
import os 
import sys
import fileinput
import math
import re
import csv
from PyQt4 import QtCore, QtGui

class gcContent(QtGui.QMainWindow):

  def __init__(self, parent=None):
    super(gcContent, self).__init__()
  
  def gcContent(self,sourcepath,output_file,fp,sp,tp,fg):
	  #print "fp===" + str(fp) +"sp===" + str(sp) + "tp===" + str(tp) +"fg===" + str(fg)
	  sequence_start = 0;
	  sequence = "";
	  sequence_name = "";
	  G_WholeGenome = 0;
	  C_WholeGenome = 0;
	  G_FirstCodon  = 0;
	  C_FirstCodon  = 0;
	  G_SecondCodon = 0;
	  C_SecondCodon = 0;
	  G_ThirdCodon = 0;
	  C_ThirdCodon = 0;
      
	  Total_GC_WholeGenome = 0;
	  Total_GC_FirstCodon  = 0;
	  Total_GC_SecondCodon = 0;
	  Total_GC_ThirdCodon  = 0;
	  Percentage_WholeGenome = 0.00;
	  Percentage_FirstCodon = 0.00;
	  Percentage_SecondCodon = 0.00;
	  Percentage_ThirdCodon = 0.00;
	  Sequence_Divisible_By_Three = 0
	  ofile_path =output_file ;
	  linenum=0
	  writer = open(ofile_path, 'w')
	  for line in open(sourcepath,'r'): #reading the input file
		  linenum=linenum+1
		  line=line.rstrip('\n')
		  #if line != "" :
		  

		  if(line[0:1]=='>' and sequence_start==0):#checking the fasta format
			  sequence_start+=1;
			  sequence_name=line[1:]
			  
			  #print "if(line[0:1]=='>' and sequence_start==0):linenum==="+ str(linenum)
			  #print line
			  
		  elif(line[0:1]=='>'):
			  #linenum=linenum+1
			  if(sequence_start>=1 and sequence!=""):
				  sequence=re.sub(' ', '', sequence)
				  error_found=2
				  for z in range(0,len(sequence)):
					  test_string=sequence[z:]
					  ch=test_string[0:1]
					  if ch!='A' and ch!='T' and ch!='C' and ch!='G':
						error_found = 1 			  
				  if error_found ==1:
					  #print "if(sequence_start>=1 and sequence!=""):linenum==="+ str(linenum)
					  
					  success= "Error in Input File line no " + str(linenum)
					  return success
				  
				  Length = len(sequence)
				  Sequence_Divisible_By_Three = math.floor((len(sequence)/3)*3)
				  sequence_length=len(sequence)
				  sequence_start+=1;
				  
				  G_WholeGenome=sequence.count('G')
				  C_WholeGenome=sequence.count('C')
				  i = 0
				  while int(i)<Sequence_Divisible_By_Three :
					  substring=sequence[i:]
					  
					  substring_first=substring
					  substring_second=substring[1:]
					  substring_third=substring[2:]
					  
					  first_ch = substring_first[0:1]
					  second_ch = substring_second[0:1]
					  third_ch = substring_third[0:1]
					  
#counting ATCG===start						
					  if first_ch=='G' : 
						G_FirstCodon += 1
					  if first_ch=='C' :
						C_FirstCodon += 1
					  if second_ch=='G' :
						G_SecondCodon += 1
					  if second_ch=='C' :
						C_SecondCodon += 1
					  if third_ch=='G' :
						G_ThirdCodon += 1
					  if third_ch=='C' :	
						C_ThirdCodon += 1
#counting ATCG===finish						  
					  i=i+3 
					    
				  Total_GC_WholeGenome=G_WholeGenome+C_WholeGenome;
				  Percentage_WholeGenome =float((float(Total_GC_WholeGenome)/float(Length))*100)
					
				  Total_GC_FirstCodon = G_FirstCodon + C_FirstCodon
				  Percentage_FirstCodon =float((float(Total_GC_FirstCodon)/float(Sequence_Divisible_By_Three))*100);
					
				  Total_GC_SecondCodon =G_SecondCodon + C_SecondCodon;
				  Percentage_SecondCodon=float((float(Total_GC_SecondCodon)/float(Sequence_Divisible_By_Three))*100);

				  Total_GC_ThirdCodon = G_ThirdCodon + C_ThirdCodon;
				  Percentage_ThirdCodon =float((float(Total_GC_ThirdCodon)/float(Sequence_Divisible_By_Three))*100);
				  
				  
				  output_line = "Sequence name :  " + sequence_name + "\n\n"
				  writer.write(output_line)
				  
				  if fp=="1" :
					  output_line = "Total GC in First Position  " + str(Total_GC_FirstCodon) + "\n"
					  writer.write(output_line)
					  output_line = "Parcentage of GC  " + str(Percentage_FirstCodon) + "\n\n"
					  writer.write(output_line)
				  if sp=="1" :
					  output_line = "Total GC in Second Position  " + str(Total_GC_SecondCodon) + "\n"
					  writer.write(output_line)
					  output_line = "Parcentage of GC  " + str(Percentage_SecondCodon) + "\n\n"
					  writer.write(output_line)	
				  if tp=="1" :
					  output_line = "Total GC in Third Position  " + str(Total_GC_ThirdCodon) + "\n"
					  writer.write(output_line)
					  output_line = "Parcentage of GC  " + str(Percentage_ThirdCodon) + "\n\n"
					  writer.write(output_line)	
				  if fg=="1" :
					  output_line = "Total GC in Whole Genome  " + str(Total_GC_WholeGenome) + "\n"
					  writer.write(output_line)
					  output_line = "Parcentage of GC  " + str(Percentage_WholeGenome) + "\n\n"
					  writer.write(output_line)		
				  
				  sequence_name=line[1:]				  					  			  
				  sequence = ""
				  success = "yes"
				  
			  elif(sequence_start==1 and sequence==""):
				  #print "elif(sequence_start==1 and sequence==""):linenum===elif"+ str(linenum)
				  success = "Sequence is not found line no " + str(linenum)
				  break;
			  
		  elif(line[0:1]!='>' and sequence_start>=1 and sequence==""):
			  #linenum=linenum+1
			  error_found=2
			  for z in range(0,len(sequence)):
				  test_string=sequence[z:]
				  ch=test_string[0:1]
				  if ch!='A' and ch!='T' and ch!='C' and ch!='G':
					error_found = 1 			  
			  if error_found ==1:
				  #print "elif(line[0:1]!='>' and sequence_start>=1):linenum==="+ str(linenum)
				 
				  success= "Error in Input File line no " + str(linenum)
				  return success
			  else :	  			  
				  sequence=sequence + line 
			  
	  if sequence != "":
			  sequence=re.sub(' ', '', sequence)
			  Length = len(sequence)
			  Sequence_Divisible_By_Three = math.floor((len(sequence)/3)*3)
			  sequence_length=len(sequence)
			  sequence_start+=1;
			  
			  G_WholeGenome=sequence.count('G',0,Length)
			  C_WholeGenome=sequence.count('C',0,Length)
			  i = 0
			  while int(i)<int(Sequence_Divisible_By_Three) :
				  substring=sequence[i:]
				  
				  substring_first=substring
				  substring_second=substring[1:]
				  substring_third=substring[2:]
				  
				  first_ch = substring_first[0:1]
				  second_ch = substring_second[0:1]
				  third_ch = substring_third[0:1]
				  					
				  if first_ch=='G' :
					G_FirstCodon += 1
				  if first_ch=='C' :
					C_FirstCodon += 1
				  if second_ch=='G' :
					G_SecondCodon += 1
				  if second_ch=='C' :
					  C_SecondCodon += 1
				  if third_ch=='G' :
					  G_ThirdCodon += 1
				  if third_ch=='C' :
					  C_ThirdCodon += 1
				  i=i+3  
			  
			  Total_GC_WholeGenome=G_WholeGenome+C_WholeGenome;
			  Percentage_WholeGenome =float((float(Total_GC_WholeGenome)/float(Length))*100)
				
			  Total_GC_FirstCodon = G_FirstCodon + C_FirstCodon
			  Percentage_FirstCodon =float((float(Total_GC_FirstCodon)/float(Sequence_Divisible_By_Three))*100);
				
			  Total_GC_SecondCodon =G_SecondCodon + C_SecondCodon;
			  Percentage_SecondCodon=float((float(Total_GC_SecondCodon)/float(Sequence_Divisible_By_Three))*100);

			  Total_GC_ThirdCodon = G_ThirdCodon + C_ThirdCodon;
			  Percentage_ThirdCodon =float((float(Total_GC_ThirdCodon)/float(Sequence_Divisible_By_Three))*100);
			  
			  
			  output_line = "Sequence name :  " + sequence_name + "\n\n"
			  writer.write(output_line)
			  
			  if fp=="1" :
				  output_line = "Total GC in First Position  " + str(Total_GC_FirstCodon) + "\n"
				  writer.write(output_line)
				  output_line = "Parcentage of GC  " + str(Percentage_FirstCodon) + "\n\n"
				  writer.write(output_line)
			  if sp=="1" :
				  output_line = "Total GC in Second Position  " + str(Total_GC_SecondCodon) + "\n"
				  writer.write(output_line)
				  output_line = "Parcentage of GC  " + str(Percentage_SecondCodon) + "\n\n"
				  writer.write(output_line)	
			  if tp=="1" :
				  output_line = "Total GC in Third Position  " + str(Total_GC_ThirdCodon) + "\n"
				  writer.write(output_line)
				  output_line = "Parcentage of GC  " + str(Percentage_ThirdCodon) + "\n\n"
				  writer.write(output_line)	
			  if fg=="1" :
				  output_line = "Total GC in Whole Genome  " + str(Total_GC_WholeGenome) + "\n"
				  writer.write(output_line)
				  output_line = "Parcentage of GC  " + str(Percentage_WholeGenome) + "\n\n"
				  writer.write(output_line)		
			  sequence = ""
			  success = "yes" 
	  
	  elif sequence == "" and sequence_name != "" :
		  success = "Sequence is not found line no " + str(linenum)
	  
	  writer.close()
	  return success
		  
