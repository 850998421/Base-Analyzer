"""
Tools Name : Base Analyzer Tool
Original Developer: Afsana Rahman Snigdha
Original Development Date : 01-01-2012
Version : 1.0
Purpose : 
		Input File: A text file containing one or more fasta sequence.
		Output File : A csv,txt file according to codon usage table for each sequence including a RSCU Value compare with base sequence.
		Methods :
			Calculation of RSCU = (total number of indivisual codon found in the sequence*total number of codon that produce same amino acid)/
							total number of different codon in whole sequence.
					Example : UUA is a codon that produce Leu amino acid and there are five other codons which produce the same amino acid.
								if in a sequence UUA codon found 10 times and total number of different codon are 50
								then RSCU = (10*6)/50 
								N.B: 6 is the number of codons which produce same amino acid.
															
								
"""

import pdb
import os 
import sys
import fileinput
import math
import re
import csv
from PyQt4 import QtCore, QtGui

class codonUsaseTableWithoutComperision(QtGui.QMainWindow):

  def __init__(self, parent=None):
    super(codonUsaseTableWithoutComperision, self).__init__()
  
  def codonUsaseTableWithoutComperision(self,sourcepath,output_file):
	  base_sequence_start = 0;
	  compare_sequence_start = 0;
	  base_sequence = "";
	  compare_sequence = "";
	  base_sequence_name = "";
	  compare_sequence_name = "";
	  ofile_path =output_file ;
	  success = "" 
	  linenum = 0; 
	  for line in open(sourcepath,'r'):
		  line=line.rstrip('\n')
		  linenum = linenum + 1; 
		  if(line[0:1]=='>' and base_sequence_start==0):
			  base_sequence_start+=1;
			  base_sequence_name=line[1:]
		  elif(line[0:1]=='>'):
			  if(base_sequence_start==1 and base_sequence!=""):
				  base_sequence=re.sub('T', 'U', base_sequence)
				  base_sequence=re.sub(' ', '', base_sequence)
				  base_sequence=base_sequence.upper()
				  base_sequence_length=math.floor((len(base_sequence)/3)*3)
				  base_sequence=base_sequence[0:int(base_sequence_length)]
				  base_sequence_start+=1;
				  
				  # declaring codon table array
				  base_array={
							'U_array':{
										'UUU':{'number' : 0,'RNA' : 'Phe','RSCU' : 0.00},
										'UUC':{'number' : 0,'RNA' : 'Phe','RSCU' : 0.00},
										'UUA':{'number' : 0,'RNA' : 'Leu','RSCU' : 0.00},
										'UUG':{'number' : 0,'RNA' : 'Leu','RSCU' : 0.00},
										'UCU':{'number' : 0,'RNA' : 'Ser','RSCU' : 0.00},
										'UCC':{'number' : 0,'RNA' : 'Ser','RSCU' : 0.00},
										'UCA':{'number' : 0,'RNA' : 'Ser','RSCU' : 0.00},
										'UCG':{'number' : 0,'RNA' : 'Ser','RSCU' : 0.00},
										'UAU':{'number' : 0,'RNA' : 'Tyr','RSCU' : 0.00},
										'UAC':{'number' : 0,'RNA' : 'Tyr','RSCU' : 0.00},
										'UAA':{'number' : 0,'RNA' : 'Stop','RSCU' : 0.00},
										'UAG':{'number' : 0,'RNA' : 'Stop','RSCU' : 0.00},
										'UGU':{'number' : 0,'RNA' : 'Cys','RSCU' : 0.00},
										'UGC':{'number' : 0,'RNA' : 'Cys','RSCU' : 0.00},
										'UGA':{'number' : 0,'RNA' : 'Stop','RSCU' : 0.00},
										'UGG':{'number' : 0,'RNA' : 'Trp','RSCU' : 0.00}
										},
							'G_array':{
										'GUU':{'number' : 0,'RNA' : 'Val','RSCU' :  0.00}, 
										'GUC':{'number' : 0,'RNA' : 'Val','RSCU' :  0.00}, 
										'GUA':{'number' : 0,'RNA' : 'Val','RSCU' :  0.00}, 
										'GUG':{'number' : 0,'RNA' : 'Val','RSCU' :  0.00}, 
										'GCU':{'number' : 0,'RNA' : 'Ala','RSCU' :  0.00}, 
										'GCC':{'number' : 0,'RNA' : 'Ala','RSCU' :  0.00}, 
										'GCA':{'number' : 0,'RNA' : 'Ala','RSCU' :  0.00}, 
										'GCG':{'number' : 0,'RNA' : 'Ala','RSCU' :  0.00}, 
										'GAU':{'number' : 0,'RNA' : 'Asp','RSCU' :  0.00}, 
										'GAC':{'number' : 0,'RNA' : 'Asp','RSCU' :  0.00}, 
										'GAA':{'number' : 0,'RNA' : 'Glu','RSCU' :  0.00}, 
										'GAG':{'number' : 0,'RNA' : 'Glu','RSCU' :  0.00}, 
										'GGU':{'number' : 0,'RNA' : 'Gly','RSCU' :  0.00}, 
										'GGC':{'number' : 0,'RNA' : 'Gly','RSCU' :  0.00}, 
										'GGA':{'number' : 0,'RNA' : 'Gly','RSCU' :  0.00}, 
										'GGG':{'number' : 0,'RNA' : 'Gly','RSCU' :  0.00}
										},
							'C_array':{
										'CUU':{'number' : 0,'RNA' : 'Leu','RSCU' :  0.00},
										'CUC':{'number' : 0,'RNA' : 'Leu','RSCU' :  0.00},
										'CUA':{'number' : 0,'RNA' : 'Leu','RSCU' :  0.00},
										'CUG':{'number' : 0,'RNA' : 'Leu','RSCU' :  0.00},
										'CCU':{'number' : 0,'RNA' : 'Pro','RSCU' :  0.00},
										'CCC':{'number' : 0,'RNA' : 'Pro','RSCU' :  0.00},
										'CCA':{'number' : 0,'RNA' : 'Pro','RSCU' :  0.00},
										'CCG':{'number' : 0,'RNA' : 'Pro','RSCU' :  0.00},
										'CAU':{'number' : 0,'RNA' : 'His','RSCU' :  0.00},
										'CAC':{'number' : 0,'RNA' : 'His','RSCU' :  0.00},
										'CAA':{'number' : 0,'RNA' : 'Gln','RSCU' :  0.00},
										'CAG':{'number' : 0,'RNA' : 'Gln','RSCU' :  0.00},
										'CGU':{'number' : 0,'RNA' : 'Arg','RSCU' :  0.00},
										'CGC':{'number' : 0,'RNA' : 'Arg','RSCU' :  0.00},
										'CGA':{'number' : 0,'RNA' : 'Arg','RSCU' :  0.00},
										'CGG':{'number' : 0,'RNA' : 'Arg','RSCU' :  0.00}
										},
							'A_array':{
										'AUU':{'number' : 0,'RNA' : 'Ile','RSCU' :  0.00} ,
										'AUC':{'number' : 0,'RNA' : 'Ile','RSCU' :  0.00} ,
										'AUA':{'number' : 0,'RNA' : 'Ile','RSCU' :  0.00} ,
										'AUG':{'number' : 0,'RNA' : 'Met','RSCU' :  0.00} ,
										'ACU':{'number' : 0,'RNA' : 'Thr','RSCU' :  0.00} ,
										'ACC':{'number' : 0,'RNA' : 'Thr','RSCU' :  0.00} ,
										'ACA':{'number' : 0,'RNA' : 'Thr','RSCU' :  0.00} ,
										'ACG':{'number' : 0,'RNA' : 'Thr','RSCU' :  0.00} ,
										'AAU':{'number' : 0,'RNA' : 'Asn','RSCU' :  0.00} ,
										'AAC':{'number' : 0,'RNA' : 'Asn','RSCU' :  0.00} ,
										'AAA':{'number' : 0,'RNA' : 'Lys','RSCU' :  0.00} ,
										'AAG':{'number' : 0,'RNA' : 'Lys','RSCU' :  0.00} ,
										'AGU':{'number' : 0,'RNA' : 'Ser','RSCU' :  0.00} ,
										'AGC':{'number' : 0,'RNA' : 'Ser','RSCU' :  0.00} ,
										'AGA':{'number' : 0,'RNA' : 'Arg','RSCU' :  0.00} ,
										'AGG':{'number' : 0,'RNA' : 'Arg','RSCU' :  0.00} 
										},
							'Amino_acid_array':{
												'Phe':{'alt_codon_num' : 2,'total_codon' : 0} ,
												'Leu':{'alt_codon_num' : 6,'total_codon' : 0} ,
												'Ile':{'alt_codon_num' : 3,'total_codon' : 0} ,
												'Met':{'alt_codon_num' : 1,'total_codon' : 0} ,
												'Val':{'alt_codon_num' : 4,'total_codon' : 0} ,
												'Ser':{'alt_codon_num' : 6,'total_codon' : 0} ,
												'Pro':{'alt_codon_num' : 4,'total_codon' : 0} ,
												'Thr':{'alt_codon_num' : 4,'total_codon' : 0} ,
												'Ala':{'alt_codon_num' : 4,'total_codon' : 0} ,
												'Tyr':{'alt_codon_num' : 2,'total_codon' : 0} ,
												'Stop':{'alt_codon_num' : 3,'total_codon' : 0} ,
												'His':{'alt_codon_num' : 2,'total_codon' : 0} ,
												'Gln':{'alt_codon_num' : 2,'total_codon' : 0} ,
												'Asn':{'alt_codon_num' : 2,'total_codon' : 0} ,
												'Lys':{'alt_codon_num' : 2,'total_codon' : 0} ,
												'Asp':{'alt_codon_num' : 2,'total_codon' : 0} ,
												'Glu':{'alt_codon_num' : 2,'total_codon' : 0} ,
												'Cys':{'alt_codon_num' : 2,'total_codon' : 0} ,
												'Trp':{'alt_codon_num' : 1,'total_codon' : 0} ,
												'Arg':{'alt_codon_num' : 6,'total_codon' : 0} ,
												'Gly':{'alt_codon_num' : 4,'total_codon' : 0} 
												}
						
							}
				  i=0;			
				  while int(i)<base_sequence_length :
							substring=base_sequence[i:]
							codon = substring[0:3]
							first_ch=codon[0:1]
							
							if(first_ch=='U'):
								if (base_array['U_array'][codon]):
									base_array['U_array'][codon]['number']=base_array['U_array'][codon]['number']+1
									
							if(first_ch=='G'):
									if (base_array['G_array'][codon]):
										base_array['G_array'][codon]['number']=base_array['G_array'][codon]['number']+1
										
							if(first_ch=='C'):
									if (base_array['C_array'][codon]):
										base_array['C_array'][codon]['number']=base_array['C_array'][codon]['number']+1
										
							if(first_ch=='A'):
									if (base_array['A_array'][codon]):
										base_array['A_array'][codon]['number']=base_array['A_array'][codon]['number']+1
							i=i+3
							
				  base_array['Amino_acid_array']['Phe']['total_codon']=base_array['U_array']['UUU']['number']+base_array['U_array']['UUC']['number'];
				  base_array['Amino_acid_array']['Leu']['total_codon']=base_array['U_array']['UUA']['number']+base_array['U_array']['UUG']['number']+base_array['C_array']['CUU']['number']+base_array['C_array']['CUC']['number']+base_array['C_array']['CUA']['number']+base_array['C_array']['CUG']['number'];
				  base_array['Amino_acid_array']['Ile']['total_codon']=base_array['A_array']['AUU']['number']+base_array['A_array']['AUC']['number']+base_array['A_array']['AUA']['number'];
				  base_array['Amino_acid_array']['Met']['total_codon']=base_array['A_array']['AUG']['number'];
				  base_array['Amino_acid_array']['Val']['total_codon']=base_array['G_array']['GUU']['number']+base_array['G_array']['GUC']['number']+base_array['G_array']['GUA']['number']+base_array['G_array']['GUG']['number'];
				  base_array['Amino_acid_array']['Ser']['total_codon']=base_array['U_array']['UCU']['number']+base_array['U_array']['UCC']['number']+base_array['U_array']['UCA']['number']+base_array['U_array']['UCG']['number']+base_array['A_array']['AGU']['number']+base_array['A_array']['AGC']['number'];
				  base_array['Amino_acid_array']['Pro']['total_codon']=base_array['C_array']['CCU']['number']+base_array['C_array']['CCC']['number']+base_array['C_array']['CCA']['number']+base_array['C_array']['CCG']['number'];
				  base_array['Amino_acid_array']['Thr']['total_codon']=base_array['A_array']['ACU']['number']+base_array['A_array']['ACC']['number']+base_array['A_array']['ACA']['number']+base_array['A_array']['ACG']['number'];
				  base_array['Amino_acid_array']['Ala']['total_codon']=base_array['G_array']['GCU']['number']+base_array['G_array']['GCC']['number']+base_array['G_array']['GCA']['number']+base_array['G_array']['GCG']['number'];
				  base_array['Amino_acid_array']['Tyr']['total_codon']=base_array['U_array']['UAU']['number']+base_array['U_array']['UAC']['number'];
				  base_array['Amino_acid_array']['Stop']['total_codon']=base_array['U_array']['UAA']['number']+base_array['U_array']['UGA']['number']+base_array['U_array']['UAG']['number'];
				  base_array['Amino_acid_array']['His']['total_codon']=base_array['C_array']['CAU']['number']+base_array['C_array']['CAC']['number'];
				  base_array['Amino_acid_array']['Gln']['total_codon']=base_array['C_array']['CAA']['number']+base_array['C_array']['CAG']['number'];
				  base_array['Amino_acid_array']['Asn']['total_codon']=base_array['A_array']['AAU']['number']+base_array['A_array']['AAC']['number'];
				  base_array['Amino_acid_array']['Lys']['total_codon']=base_array['A_array']['AAA']['number']+base_array['A_array']['AAG']['number'];
				  base_array['Amino_acid_array']['Asp']['total_codon']=base_array['G_array']['GAU']['number']+base_array['G_array']['GAC']['number'];
				  base_array['Amino_acid_array']['Glu']['total_codon']=base_array['G_array']['GAA']['number']+base_array['G_array']['GAG']['number'];
				  base_array['Amino_acid_array']['Cys']['total_codon']=base_array['U_array']['UGU']['number']+base_array['U_array']['UGC']['number'];
				  base_array['Amino_acid_array']['Trp']['total_codon']=base_array['U_array']['UGG']['number'];
				  base_array['Amino_acid_array']['Arg']['total_codon']=base_array['C_array']['CGU']['number']+base_array['C_array']['CGC']['number']+base_array['C_array']['CGA']['number']+base_array['C_array']['CGG']['number']+base_array['A_array']['AGA']['number']+base_array['A_array']['AGG']['number'];
				  base_array['Amino_acid_array']['Gly']['total_codon']=base_array['G_array']['GGU']['number']+base_array['G_array']['GGC']['number']+base_array['G_array']['GGA']['number']+base_array['G_array']['GGG']['number'];
				  
				  i=0;			
				  while int(i)<base_sequence_length :
						substring=base_sequence[i:]
						codon = substring[0:3]
						first_ch=codon[0:1]
						
						if(first_ch=='U'):
							if (base_array['U_array'][codon]):
								amio_acid_name=base_array['U_array'][codon]['RNA'];
								if base_array['Amino_acid_array'][amio_acid_name]['total_codon']>0:								
									base_array['U_array'][codon]['RSCU']=(base_array['U_array'][codon]['number']*base_array['Amino_acid_array'][amio_acid_name]['alt_codon_num'])/base_array['Amino_acid_array'][amio_acid_name]['total_codon']
									base_array['U_array'][codon]['RSCU']=round(base_array['U_array'][codon]['RSCU'],3)
								
						if(first_ch=='G'):
							if (base_array['G_array'][codon]):
								amio_acid_name=base_array['G_array'][codon]['RNA'];
								if base_array['Amino_acid_array'][amio_acid_name]['total_codon']>0:								
									base_array['G_array'][codon]['RSCU']=(base_array['G_array'][codon]['number']*base_array['Amino_acid_array'][amio_acid_name]['alt_codon_num'])/base_array['Amino_acid_array'][amio_acid_name]['total_codon'];
									base_array['G_array'][codon]['RSCU']=round(base_array['G_array'][codon]['RSCU'],3)
								
						if(first_ch=='C'):
							if (base_array['C_array'][codon]):
								amio_acid_name=base_array['C_array'][codon]['RNA'];
								if base_array['Amino_acid_array'][amio_acid_name]['total_codon']>0:								
									base_array['C_array'][codon]['RSCU']=(base_array['C_array'][codon]['number']*base_array['Amino_acid_array'][amio_acid_name]['alt_codon_num'])/base_array['Amino_acid_array'][amio_acid_name]['total_codon'];
									base_array['C_array'][codon]['RSCU']=round(base_array['C_array'][codon]['RSCU'],3)
								
						if(first_ch=='A'):
							if (base_array['A_array'][codon]):
								amio_acid_name=base_array['A_array'][codon]['RNA'];
								if base_array['Amino_acid_array'][amio_acid_name]['total_codon']>0:								
									base_array['A_array'][codon]['RSCU']=(base_array['A_array'][codon]['number']*base_array['Amino_acid_array'][amio_acid_name]['alt_codon_num'])/base_array['Amino_acid_array'][amio_acid_name]['total_codon'];
									base_array['A_array'][codon]['RSCU']=round(base_array['A_array'][codon]['RSCU'],3)
						i=i+3

				  writer = open(ofile_path, 'wb')
				  output_line = "Sequence name :  " + base_sequence_name + "\n\n"
				  writer.write(output_line)
				  
				  output_line = "Phe"+"\t"+"UUU"+"\t"+ str(base_array['U_array']['UUU']['number'])+"\t"+ str(base_array['U_array']['UUU']['RSCU'])
				  output_line =output_line+"\t"+"Ser"+"\t"+"UCU"+"\t"+str(base_array['U_array']['UCU']['number'])+"\t"+str(base_array['U_array']['UCU']['RSCU'])
				  output_line =output_line+"\t"+"Tyr"+"\t"+"UAU"+"\t"+str(base_array['U_array']['UAU']['number'])+"\t"+str(base_array['U_array']['UAU']['RSCU'])
				  output_line =output_line+"\t"+"Cys"+"\t"+"UGU"+"\t"+str(base_array['U_array']['UGU']['number'])+"\t"+str(base_array['U_array']['UGU']['RSCU'])+ "\n"
				  writer.write(output_line)
				  
				  output_line = ""+"\t"+"UUC"+"\t"+str(base_array['U_array']['UUC']['number'])+"\t"+str(base_array['U_array']['UUC']['RSCU'])
				  output_line =output_line+"\t"+"\t"+"UCC"+"\t"+str(base_array['U_array']['UCC']['number'])+"\t"+str(base_array['U_array']['UCC']['RSCU'])
				  output_line =output_line+"\t"+"\t"+"UAC"+"\t"+str(base_array['U_array']['UAC']['number'])+"\t"+str(base_array['U_array']['UAC']['RSCU'])
				  output_line =output_line+"\t"+"\t"+"UGC"+"\t"+str(base_array['U_array']['UGC']['number'])+"\t"+str(base_array['U_array']['UGC']['RSCU'])+ "\n"
				  writer.write(output_line)
				  
				  output_line = "Leu"+"\t"+"UUA"+"\t"+str(base_array['U_array']['UUA']['number'])+"\t"+str(base_array['U_array']['UUA']['RSCU'])
				  output_line =output_line+"\t"+"\t"+"UCA"+"\t"+str(base_array['U_array']['UCA']['number'])+"\t"+str(base_array['U_array']['UCA']['RSCU'])
				  output_line =output_line+"\t"+"Stop"+"\t"+"UAA"+"\t"+str(base_array['U_array']['UAA']['number'])+"\t"+str(base_array['U_array']['UAA']['RSCU'])
				  output_line =output_line+"\t"+"Stop"+"\t"+"UGA"+"\t"+str(base_array['U_array']['UGA']['number'])+"\t"+str(base_array['U_array']['UGA']['RSCU'])+ "\n"
				  writer.write(output_line)
				  
				  output_line = ""+"\t"+"UUG"+"\t"+str(base_array['U_array']['UUG']['number'])+"\t"+str(base_array['U_array']['UUG']['RSCU'])
				  output_line =output_line+"\t"+"\t"+"UCG"+"\t"+str(base_array['U_array']['UCG']['number'])+"\t"+str(base_array['U_array']['UCG']['RSCU'])
				  output_line =output_line+"\t"+"\t"+"UAG"+"\t"+str(base_array['U_array']['UAG']['number'])+"\t"+str(base_array['U_array']['UAG']['RSCU'])
				  output_line =output_line+"\t"+"Trp"+"\t"+"UGG"+"\t"+str(base_array['U_array']['UGG']['number'])+"\t"+str(base_array['U_array']['UGG']['RSCU'])+ "\n"
				  writer.write(output_line)
				  
				  output_line = ""+"\t"+"CUU"+"\t"+str(base_array['C_array']['CUU']['number'])+"\t"+str(base_array['C_array']['CUU']['RSCU'])
				  output_line =output_line+"\t"+"Pro"+"\t"+"CCU"+"\t"+str(base_array['C_array']['CCU']['number'])+"\t"+str(base_array['C_array']['CCU']['RSCU'])
				  output_line =output_line+"\t"+"His"+"\t"+"CAU"+"\t"+str(base_array['C_array']['CAU']['number'])+"\t"+str(base_array['C_array']['CAU']['RSCU'])
				  output_line =output_line+"\t"+"Arg"+"\t"+"CGU"+"\t"+str(base_array['C_array']['CGU']['number'])+"\t"+str(base_array['C_array']['CGU']['RSCU'])+ "\n"
				  writer.write(output_line)
				  
				  output_line = ""+"\t"+"CUC"+"\t"+str(base_array['C_array']['CUC']['number'])+"\t"+str(base_array['C_array']['CUC']['RSCU'])
				  output_line =output_line+"\t"+"\t"+"CCC"+"\t"+str(base_array['C_array']['CCC']['number'])+"\t"+str(base_array['C_array']['CCC']['RSCU'])
				  output_line =output_line+"\t"+"\t"+"CAC"+"\t"+str(base_array['C_array']['CAC']['number'])+"\t"+str(base_array['C_array']['CAC']['RSCU'])
				  output_line =output_line+"\t"+"\t"+"CGC"+"\t"+str(base_array['C_array']['CGC']['number'])+"\t"+str(base_array['C_array']['CGC']['RSCU'])+ "\n"
				  writer.write(output_line)
				  
				  output_line = ""+"\t"+"CUA"+"\t"+str(base_array['C_array']['CUA']['number'])+"\t"+str(base_array['C_array']['CUA']['RSCU'])
				  output_line =output_line+"\t"+"\t"+"CCA"+"\t"+str(base_array['C_array']['CCA']['number'])+"\t"+str(base_array['C_array']['CCA']['RSCU'])
				  output_line =output_line+"\t"+"Gln"+"\t"+"CAA"+"\t"+str(base_array['C_array']['CAA']['number'])+"\t"+str(base_array['C_array']['CAA']['RSCU'])
				  output_line =output_line+"\t"+"\t"+"CAG"+"\t"+str(base_array['C_array']['CAG']['number'])+"\t"+str(base_array['C_array']['CAG']['RSCU'])+ "\n"
				  writer.write(output_line)
				  
				  output_line = ""+"\t"+"CUG"+"\t"+str(base_array['C_array']['CUG']['number'])+"\t"+str(base_array['C_array']['CUG']['RSCU'])
				  output_line = output_line + "\t"+"\t"+"CCG"+"\t"+str(base_array['C_array']['CCG']['number'])+"\t"+str(base_array['C_array']['CCG']['RSCU'])
				  output_line = output_line + "\t"+"\t"+"CAG"+"\t"+str(base_array['C_array']['CAG']['number'])+"\t"+str(base_array['C_array']['CAG']['RSCU'])
				  output_line = output_line + "\t"+"\t"+"CGG"+"\t"+str(base_array['C_array']['CGG']['number'])+"\t"+str(base_array['C_array']['CGG']['RSCU'])+ "\n"
				  writer.write(output_line)
				  
				  output_line = "Iso"+"\t"+"AUU"+"\t"+str(base_array['A_array']['AUU']['number'])+"\t"+str(base_array['A_array']['AUU']['RSCU'])
				  output_line =output_line+"\t"+"Thr"+"\t"+"ACU"+"\t"+str(base_array['A_array']['ACU']['number'])+"\t"+str(base_array['A_array']['ACU']['RSCU'])
				  output_line =output_line+"\t"+"Asn"+"\t"+"AAU"+"\t"+str(base_array['A_array']['AAU']['number'])+"\t"+str(base_array['A_array']['AAU']['RSCU'])
				  output_line =output_line+"\t"+"Ser"+"\t"+"AGU"+"\t"+str(base_array['A_array']['AGU']['number'])+"\t"+str(base_array['A_array']['AGU']['RSCU'])+ "\n"
				  writer.write(output_line)
				  
				  output_line = ""+"\t"+"AUC"+"\t"+str(base_array['A_array']['AUC']['number'])+"\t"+str(base_array['A_array']['AUC']['RSCU'])
				  output_line =output_line+"\t"+"\t"+"ACC"+"\t"+str(base_array['A_array']['ACC']['number'])+"\t"+str(base_array['A_array']['ACC']['RSCU'])
				  output_line =output_line+"\t"+"\t"+"AAC"+"\t"+str(base_array['A_array']['AAC']['number'])+"\t"+str(base_array['A_array']['AAC']['RSCU'])
				  output_line =output_line+"\t"+"\t"+"AGC"+"\t"+str(base_array['A_array']['AGC']['number'])+"\t"+str(base_array['A_array']['AGC']['RSCU'])+ "\n"
				  writer.write(output_line)
				  
				  output_line = ""+"\t"+"AUA"+"\t"+str(base_array['A_array']['AUA']['number'])+"\t"+str(base_array['A_array']['AUA']['RSCU'])
				  output_line =output_line+"\t"+"\t"+"ACA"+"\t"+str(base_array['A_array']['ACA']['number'])+"\t"+str(base_array['A_array']['ACA']['RSCU'])
				  output_line =output_line+"\t"+"Lys"+"\t"+"AAA"+"\t"+str(base_array['A_array']['AAA']['number'])+"\t"+str(base_array['A_array']['AAA']['RSCU'])
				  output_line =output_line+"\t"+"Arg"+"\t"+"AGA"+"\t"+str(base_array['A_array']['AGA']['number'])+"\t"+str(base_array['A_array']['AGA']['RSCU'])+ "\n"
				  writer.write(output_line)
				  
				  output_line = "Met"+"\t"+"AUG"+"\t"+str(base_array['A_array']['AUG']['number'])+"\t"+str(base_array['A_array']['AUG']['RSCU'])
				  output_line =output_line+"\t"+"\t"+"ACG"+"\t"+str(base_array['A_array']['ACG']['number'])+"\t"+str(base_array['A_array']['ACG']['RSCU'])
				  output_line =output_line+"\t"+"\t"+"AAG"+"\t"+str(base_array['A_array']['AAG']['number'])+"\t"+str(base_array['A_array']['AAG']['RSCU'])
				  output_line =output_line+"\t"+"\t"+"AGG"+"\t"+str(base_array['A_array']['AGG']['number'])+"\t"+str(base_array['A_array']['AGG']['RSCU'])+ "\n"
				  writer.write(output_line)
				  
				  output_line = "Val"+"\t"+"GUU"+"\t"+str(base_array['G_array']['GUU']['number'])+"\t"+str(base_array['G_array']['GUU']['RSCU'])
				  output_line =output_line+"\t"+"Ala"+"\t"+"GCU"+"\t"+str(base_array['G_array']['GCU']['number'])+"\t"+str(base_array['G_array']['GCU']['RSCU'])
				  output_line =output_line+"\t"+"Asp"+"\t"+"GAU"+"\t"+str(base_array['G_array']['GAU']['number'])+"\t"+str(base_array['G_array']['GAU']['RSCU'])
				  output_line =output_line+"\t"+"Gly"+"\t"+"GGU"+"\t"+str(base_array['G_array']['GGU']['number'])+"\t"+str(base_array['G_array']['GGU']['RSCU'])+ "\n"
				  writer.write(output_line)
				  
				  output_line = ""+"\t"+"GUC"+"\t"+str(base_array['G_array']['GUC']['number'])+"\t"+str(base_array['G_array']['GUC']['RSCU'])
				  output_line =output_line+"\t"+"\t"+"GCC"+"\t"+str(base_array['G_array']['GCC']['number'])+"\t"+str(base_array['G_array']['GCC']['RSCU'])
				  output_line =output_line+"\t"+"\t"+"GAC"+"\t"+str(base_array['G_array']['GAC']['number'])+"\t"+str(base_array['G_array']['GAC']['RSCU'])
				  output_line =output_line+"\t"+"\t"+"GGC"+"\t"+str(base_array['G_array']['GGC']['number'])+"\t"+str(base_array['G_array']['GGC']['RSCU'])+ "\n"
				  writer.write(output_line)
				  
				  output_line = ""+"\t"+"GUA"+"\t"+str(base_array['G_array']['GUA']['number'])+"\t"+str(base_array['G_array']['GUA']['RSCU'])
				  output_line =output_line+"\t"+"\t"+"GCA"+"\t"+str(base_array['G_array']['GCA']['number'])+"\t"+str(base_array['G_array']['GCA']['RSCU'])
				  output_line =output_line+"\t"+"Glu"+"\t"+"GAA"+"\t"+str(base_array['G_array']['GAA']['number'])+"\t"+str(base_array['G_array']['GAA']['RSCU'])
				  output_line =output_line+"\t"+"\t"+"GGA"+"\t"+str(base_array['G_array']['GGA']['number'])+"\t"+str(base_array['G_array']['GGA']['RSCU'])+ "\n"
				  writer.write(output_line)
				  
				  output_line = ""+"\t"+"GUG"+"\t"+str(base_array['G_array']['GUG']['number'])+"\t"+str(base_array['G_array']['GUG']['RSCU'])
				  output_line =output_line+"\t"+"\t"+"GCG"+"\t"+str(base_array['G_array']['GCG']['number'])+"\t"+str(base_array['G_array']['GCG']['RSCU'])
				  output_line =output_line+"\t"+"\t"+"GAG"+"\t"+str(base_array['G_array']['GAG']['number'])+"\t"+str(base_array['G_array']['GAG']['RSCU'])
				  output_line =output_line+"\t"+"\t"+"GGG"+"\t"+str(base_array['G_array']['GGG']['number'])+"\t"+str(base_array['G_array']['GGG']['RSCU'])+ "\n\n"
				  writer.write(output_line)
				  
				  base_sequence = ""
				  base_sequence_start = 2
			  elif(base_sequence_start==1 and base_sequence==""):
				  success = "Base sequence is not found in line number " + str(linenum)
				  break;
			  elif(compare_sequence_start>=1 and compare_sequence==""): 
				  success = "Compare sequence is not found in line number  " + str(linenum)
				  break;
			  elif(compare_sequence_start>=1 and compare_sequence!=""):
				  compare_sequence=re.sub('T', 'U', compare_sequence)
				  compare_sequence=re.sub(' ', '', compare_sequence)
				  compare_sequence=compare_sequence.upper()
				  compare_sequence_length=math.floor((len(compare_sequence)/3)*3)
				  compare_sequence=compare_sequence[0:int(compare_sequence_length)]
				  compare_sequence_array={
							'U_array':{
										'UUU':{'number' : 0,'RNA' : 'Phe','RSCU' : 0.00},
										'UUC':{'number' : 0,'RNA' : 'Phe','RSCU' : 0.00},
										'UUA':{'number' : 0,'RNA' : 'Leu','RSCU' : 0.00},
										'UUG':{'number' : 0,'RNA' : 'Leu','RSCU' : 0.00},
										'UCU':{'number' : 0,'RNA' : 'Ser','RSCU' : 0.00},
										'UCC':{'number' : 0,'RNA' : 'Ser','RSCU' : 0.00},
										'UCA':{'number' : 0,'RNA' : 'Ser','RSCU' : 0.00},
										'UCG':{'number' : 0,'RNA' : 'Ser','RSCU' : 0.00},
										'UAU':{'number' : 0,'RNA' : 'Tyr','RSCU' : 0.00},
										'UAC':{'number' : 0,'RNA' : 'Tyr','RSCU' : 0.00},
										'UAA':{'number' : 0,'RNA' : 'Stop','RSCU' : 0.00},
										'UAG':{'number' : 0,'RNA' : 'Stop','RSCU' : 0.00},
										'UGU':{'number' : 0,'RNA' : 'Cys','RSCU' : 0.00},
										'UGC':{'number' : 0,'RNA' : 'Cys','RSCU' : 0.00},
										'UGA':{'number' : 0,'RNA' : 'Stop','RSCU' : 0.00},
										'UGG':{'number' : 0,'RNA' : 'Trp','RSCU' : 0.00}
										},
							'G_array':{
										'GUU':{'number' : 0,'RNA' : 'Val','RSCU' :  0.00}, 
										'GUC':{'number' : 0,'RNA' : 'Val','RSCU' :  0.00}, 
										'GUA':{'number' : 0,'RNA' : 'Val','RSCU' :  0.00}, 
										'GUG':{'number' : 0,'RNA' : 'Val','RSCU' :  0.00}, 
										'GCU':{'number' : 0,'RNA' : 'Ala','RSCU' :  0.00}, 
										'GCC':{'number' : 0,'RNA' : 'Ala','RSCU' :  0.00}, 
										'GCA':{'number' : 0,'RNA' : 'Ala','RSCU' :  0.00}, 
										'GCG':{'number' : 0,'RNA' : 'Ala','RSCU' :  0.00}, 
										'GAU':{'number' : 0,'RNA' : 'Asp','RSCU' :  0.00}, 
										'GAC':{'number' : 0,'RNA' : 'Asp','RSCU' :  0.00}, 
										'GAA':{'number' : 0,'RNA' : 'Glu','RSCU' :  0.00}, 
										'GAG':{'number' : 0,'RNA' : 'Glu','RSCU' :  0.00}, 
										'GGU':{'number' : 0,'RNA' : 'Gly','RSCU' :  0.00}, 
										'GGC':{'number' : 0,'RNA' : 'Gly','RSCU' :  0.00}, 
										'GGA':{'number' : 0,'RNA' : 'Gly','RSCU' :  0.00}, 
										'GGG':{'number' : 0,'RNA' : 'Gly','RSCU' :  0.00}
										},
							'C_array':{
										'CUU':{'number' : 0,'RNA' : 'Leu','RSCU' :  0.00},
										'CUC':{'number' : 0,'RNA' : 'Leu','RSCU' :  0.00},
										'CUA':{'number' : 0,'RNA' : 'Leu','RSCU' :  0.00},
										'CUG':{'number' : 0,'RNA' : 'Leu','RSCU' :  0.00},
										'CCU':{'number' : 0,'RNA' : 'Pro','RSCU' :  0.00},
										'CCC':{'number' : 0,'RNA' : 'Pro','RSCU' :  0.00},
										'CCA':{'number' : 0,'RNA' : 'Pro','RSCU' :  0.00},
										'CCG':{'number' : 0,'RNA' : 'Pro','RSCU' :  0.00},
										'CAU':{'number' : 0,'RNA' : 'His','RSCU' :  0.00},
										'CAC':{'number' : 0,'RNA' : 'His','RSCU' :  0.00},
										'CAA':{'number' : 0,'RNA' : 'Gln','RSCU' :  0.00},
										'CAG':{'number' : 0,'RNA' : 'Gln','RSCU' :  0.00},
										'CGU':{'number' : 0,'RNA' : 'Arg','RSCU' :  0.00},
										'CGC':{'number' : 0,'RNA' : 'Arg','RSCU' :  0.00},
										'CGA':{'number' : 0,'RNA' : 'Arg','RSCU' :  0.00},
										'CGG':{'number' : 0,'RNA' : 'Arg','RSCU' :  0.00}
										},
							'A_array':{
										'AUU':{'number' : 0,'RNA' : 'Ile','RSCU' :  0.00} ,
										'AUC':{'number' : 0,'RNA' : 'Ile','RSCU' :  0.00} ,
										'AUA':{'number' : 0,'RNA' : 'Ile','RSCU' :  0.00} ,
										'AUG':{'number' : 0,'RNA' : 'Met','RSCU' :  0.00} ,
										'ACU':{'number' : 0,'RNA' : 'Thr','RSCU' :  0.00} ,
										'ACC':{'number' : 0,'RNA' : 'Thr','RSCU' :  0.00} ,
										'ACA':{'number' : 0,'RNA' : 'Thr','RSCU' :  0.00} ,
										'ACG':{'number' : 0,'RNA' : 'Thr','RSCU' :  0.00} ,
										'AAU':{'number' : 0,'RNA' : 'Asn','RSCU' :  0.00} ,
										'AAC':{'number' : 0,'RNA' : 'Asn','RSCU' :  0.00} ,
										'AAA':{'number' : 0,'RNA' : 'Lys','RSCU' :  0.00} ,
										'AAG':{'number' : 0,'RNA' : 'Lys','RSCU' :  0.00} ,
										'AGU':{'number' : 0,'RNA' : 'Ser','RSCU' :  0.00} ,
										'AGC':{'number' : 0,'RNA' : 'Ser','RSCU' :  0.00} ,
										'AGA':{'number' : 0,'RNA' : 'Arg','RSCU' :  0.00} ,
										'AGG':{'number' : 0,'RNA' : 'Arg','RSCU' :  0.00} 
										},
							'Amino_acid_array':{
												'Phe':{'alt_codon_num' : 2,'total_codon' : 0} ,
												'Leu':{'alt_codon_num' : 6,'total_codon' : 0} ,
												'Ile':{'alt_codon_num' : 3,'total_codon' : 0} ,
												'Met':{'alt_codon_num' : 1,'total_codon' : 0} ,
												'Val':{'alt_codon_num' : 4,'total_codon' : 0} ,
												'Ser':{'alt_codon_num' : 6,'total_codon' : 0} ,
												'Pro':{'alt_codon_num' : 4,'total_codon' : 0} ,
												'Thr':{'alt_codon_num' : 4,'total_codon' : 0} ,
												'Ala':{'alt_codon_num' : 4,'total_codon' : 0} ,
												'Tyr':{'alt_codon_num' : 2,'total_codon' : 0} ,
												'Stop':{'alt_codon_num' : 3,'total_codon' : 0} ,
												'His':{'alt_codon_num' : 2,'total_codon' : 0} ,
												'Gln':{'alt_codon_num' : 2,'total_codon' : 0} ,
												'Asn':{'alt_codon_num' : 2,'total_codon' : 0} ,
												'Lys':{'alt_codon_num' : 2,'total_codon' : 0} ,
												'Asp':{'alt_codon_num' : 2,'total_codon' : 0} ,
												'Glu':{'alt_codon_num' : 2,'total_codon' : 0} ,
												'Cys':{'alt_codon_num' : 2,'total_codon' : 0} ,
												'Trp':{'alt_codon_num' : 1,'total_codon' : 0} ,
												'Arg':{'alt_codon_num' : 6,'total_codon' : 0} ,
												'Gly':{'alt_codon_num' : 4,'total_codon' : 0} 
												}
						
							}
				  i=0;			
				  while int(i)<compare_sequence_length :
							substring=compare_sequence[i:]
							codon = substring[0:3]
							first_ch=codon[0:1]
							"""
							second_ch=codon[1:1]
							third_ch=codon[2:1]
							
							if(first_ch!='U' or first_ch!='G' or first_ch!='C' or first_ch!='A'):
								  success= "Error in Input File in " + compare_sequence_name + " sequence" 
								  return success
							if(second_ch!='U' or second_ch!='G' or second_ch!='C' or second_ch!='A'):
								  success= "Error in Input File in " + compare_sequence_name + " sequence" 
								  return success
							if(third_ch!='U' or third_ch!='G' or third_ch!='C' or third_ch!='A'):
								  success= "Error in Input File in " + compare_sequence_name + " sequence" 
								  return success							
							"""
							
							if(first_ch=='U'):
								if (compare_sequence_array['U_array'][codon]):
									compare_sequence_array['U_array'][codon]['number']=compare_sequence_array['U_array'][codon]['number']+1
									
							if(first_ch=='G'):
									if (compare_sequence_array['G_array'][codon]):
										compare_sequence_array['G_array'][codon]['number']=compare_sequence_array['G_array'][codon]['number']+1
										
							if(first_ch=='C'):
									if (compare_sequence_array['C_array'][codon]):
										compare_sequence_array['C_array'][codon]['number']=compare_sequence_array['C_array'][codon]['number']+1
										
							if(first_ch=='A'):
									if (compare_sequence_array['A_array'][codon]):
										compare_sequence_array['A_array'][codon]['number']=compare_sequence_array['A_array'][codon]['number']+1
							i=i+3
							
				  compare_sequence_array['Amino_acid_array']['Phe']['total_codon']=compare_sequence_array['U_array']['UUU']['number']+compare_sequence_array['U_array']['UUC']['number'];
				  compare_sequence_array['Amino_acid_array']['Leu']['total_codon']=compare_sequence_array['U_array']['UUA']['number']+compare_sequence_array['U_array']['UUG']['number']+compare_sequence_array['C_array']['CUU']['number']+compare_sequence_array['C_array']['CUC']['number']+compare_sequence_array['C_array']['CUA']['number']+compare_sequence_array['C_array']['CUG']['number'];
				  compare_sequence_array['Amino_acid_array']['Ile']['total_codon']=compare_sequence_array['A_array']['AUU']['number']+compare_sequence_array['A_array']['AUC']['number']+compare_sequence_array['A_array']['AUA']['number'];
				  compare_sequence_array['Amino_acid_array']['Met']['total_codon']=compare_sequence_array['A_array']['AUG']['number'];
				  compare_sequence_array['Amino_acid_array']['Val']['total_codon']=compare_sequence_array['G_array']['GUU']['number']+compare_sequence_array['G_array']['GUC']['number']+compare_sequence_array['G_array']['GUA']['number']+compare_sequence_array['G_array']['GUG']['number'];
				  compare_sequence_array['Amino_acid_array']['Ser']['total_codon']=compare_sequence_array['U_array']['UCU']['number']+compare_sequence_array['U_array']['UCC']['number']+compare_sequence_array['U_array']['UCA']['number']+compare_sequence_array['U_array']['UCG']['number']+compare_sequence_array['A_array']['AGU']['number']+compare_sequence_array['A_array']['AGC']['number'];
				  compare_sequence_array['Amino_acid_array']['Pro']['total_codon']=compare_sequence_array['C_array']['CCU']['number']+compare_sequence_array['C_array']['CCC']['number']+compare_sequence_array['C_array']['CCA']['number']+compare_sequence_array['C_array']['CCG']['number'];
				  compare_sequence_array['Amino_acid_array']['Thr']['total_codon']=compare_sequence_array['A_array']['ACU']['number']+compare_sequence_array['A_array']['ACC']['number']+compare_sequence_array['A_array']['ACA']['number']+compare_sequence_array['A_array']['ACG']['number'];
				  compare_sequence_array['Amino_acid_array']['Ala']['total_codon']=compare_sequence_array['G_array']['GCU']['number']+compare_sequence_array['G_array']['GCC']['number']+compare_sequence_array['G_array']['GCA']['number']+compare_sequence_array['G_array']['GCG']['number'];
				  compare_sequence_array['Amino_acid_array']['Tyr']['total_codon']=compare_sequence_array['U_array']['UAU']['number']+compare_sequence_array['U_array']['UAC']['number'];
				  compare_sequence_array['Amino_acid_array']['Stop']['total_codon']=compare_sequence_array['U_array']['UAA']['number']+compare_sequence_array['U_array']['UGA']['number']+compare_sequence_array['U_array']['UAG']['number'];
				  compare_sequence_array['Amino_acid_array']['His']['total_codon']=compare_sequence_array['C_array']['CAU']['number']+compare_sequence_array['C_array']['CAC']['number'];
				  compare_sequence_array['Amino_acid_array']['Gln']['total_codon']=compare_sequence_array['C_array']['CAA']['number']+compare_sequence_array['C_array']['CAG']['number'];
				  compare_sequence_array['Amino_acid_array']['Asn']['total_codon']=compare_sequence_array['A_array']['AAU']['number']+compare_sequence_array['A_array']['AAC']['number'];
				  compare_sequence_array['Amino_acid_array']['Lys']['total_codon']=compare_sequence_array['A_array']['AAA']['number']+compare_sequence_array['A_array']['AAG']['number'];
				  compare_sequence_array['Amino_acid_array']['Asp']['total_codon']=compare_sequence_array['G_array']['GAU']['number']+compare_sequence_array['G_array']['GAC']['number'];
				  compare_sequence_array['Amino_acid_array']['Glu']['total_codon']=compare_sequence_array['G_array']['GAA']['number']+compare_sequence_array['G_array']['GAG']['number'];
				  compare_sequence_array['Amino_acid_array']['Cys']['total_codon']=compare_sequence_array['U_array']['UGU']['number']+compare_sequence_array['U_array']['UGC']['number'];
				  compare_sequence_array['Amino_acid_array']['Trp']['total_codon']=compare_sequence_array['U_array']['UGG']['number'];
				  compare_sequence_array['Amino_acid_array']['Arg']['total_codon']=compare_sequence_array['C_array']['CGU']['number']+compare_sequence_array['C_array']['CGC']['number']+compare_sequence_array['C_array']['CGA']['number']+compare_sequence_array['C_array']['CGG']['number']+compare_sequence_array['A_array']['AGA']['number']+compare_sequence_array['A_array']['AGG']['number'];
				  compare_sequence_array['Amino_acid_array']['Gly']['total_codon']=compare_sequence_array['G_array']['GGU']['number']+compare_sequence_array['G_array']['GGC']['number']+compare_sequence_array['G_array']['GGA']['number']+compare_sequence_array['G_array']['GGG']['number'];
				  
				  i=0			
				  while int(i)<compare_sequence_length :
						substring=compare_sequence[i:]
						codon = substring[0:3]
						first_ch=codon[0:1]
						
						if(first_ch=='U'):
							if (compare_sequence_array['U_array'][codon]) and compare_sequence_array['Amino_acid_array'][amio_acid_name]['total_codon']>0:
								amio_acid_name=compare_sequence_array['U_array'][codon]['RNA'];
								if compare_sequence_array['Amino_acid_array'][amio_acid_name]['total_codon']>0:								
									compare_sequence_array['U_array'][codon]['RSCU']=(compare_sequence_array['U_array'][codon]['number']*compare_sequence_array['Amino_acid_array'][amio_acid_name]['alt_codon_num'])/compare_sequence_array['Amino_acid_array'][amio_acid_name]['total_codon']
									compare_sequence_array['U_array'][codon]['RSCU']=round(compare_sequence_array['U_array'][codon]['RSCU'],3)
						if(first_ch=='G'):
							if (compare_sequence_array['G_array'][codon]) and compare_sequence_array['Amino_acid_array'][amio_acid_name]['total_codon']>0:
								amio_acid_name=compare_sequence_array['G_array'][codon]['RNA'];
								if compare_sequence_array['Amino_acid_array'][amio_acid_name]['total_codon']>0:								
									compare_sequence_array['G_array'][codon]['RSCU']=(compare_sequence_array['G_array'][codon]['number']*compare_sequence_array['Amino_acid_array'][amio_acid_name]['alt_codon_num'])/compare_sequence_array['Amino_acid_array'][amio_acid_name]['total_codon'];
									compare_sequence_array['G_array'][codon]['RSCU']=round(compare_sequence_array['G_array'][codon]['RSCU'],3)
						if(first_ch=='C'):
							if (compare_sequence_array['C_array'][codon]) and compare_sequence_array['Amino_acid_array'][amio_acid_name]['total_codon']>0:
								amio_acid_name=compare_sequence_array['C_array'][codon]['RNA'];
								if compare_sequence_array['Amino_acid_array'][amio_acid_name]['total_codon']>0:								
									compare_sequence_array['C_array'][codon]['RSCU']=(compare_sequence_array['C_array'][codon]['number']*compare_sequence_array['Amino_acid_array'][amio_acid_name]['alt_codon_num'])/compare_sequence_array['Amino_acid_array'][amio_acid_name]['total_codon'];
									compare_sequence_array['C_array'][codon]['RSCU']=round(compare_sequence_array['C_array'][codon]['RSCU'],3)
						if(first_ch=='A'):
							if (compare_sequence_array['A_array'][codon]) :
								amio_acid_name=compare_sequence_array['A_array'][codon]['RNA'];
								if compare_sequence_array['Amino_acid_array'][amio_acid_name]['total_codon']>0:
									compare_sequence_array['A_array'][codon]['RSCU']=(compare_sequence_array['A_array'][codon]['number']*compare_sequence_array['Amino_acid_array'][amio_acid_name]['alt_codon_num'])/compare_sequence_array['Amino_acid_array'][amio_acid_name]['total_codon'];
									compare_sequence_array['A_array'][codon]['RSCU']=round(compare_sequence_array['A_array'][codon]['RSCU'],3)
						i=i+3
						
				  compare_sequence=""
				  output_line = "Sequence name :  " + compare_sequence_name + "\n\n"
				  writer.write(output_line)
				  
				  output_line = "Phe"+"\t"+"UUU"+"\t"+str(compare_sequence_array['U_array']['UUU']['number'])+"\t"+str(compare_sequence_array['U_array']['UUU']['RSCU'])
				  output_line =output_line+"\t"+"Ser"+"\t"+"UCU"+"\t"+str(compare_sequence_array['U_array']['UCU']['number'])+"\t"+str(compare_sequence_array['U_array']['UCU']['RSCU'])
				  output_line =output_line+"\t"+"Tyr"+"\t"+"UAU"+"\t"+str(compare_sequence_array['U_array']['UAU']['number'])+"\t"+str(compare_sequence_array['U_array']['UAU']['RSCU'])
				  output_line =output_line+"\t"+"Cys"+"\t"+"UGU"+"\t"+str(compare_sequence_array['U_array']['UGU']['number'])+"\t"+str(compare_sequence_array['U_array']['UGU']['RSCU'])+ "\n"
				  writer.write(output_line)
				  
				  output_line = ""+"\t"+"UUC"+"\t"+str(compare_sequence_array['U_array']['UUC']['number'])+"\t"+str(compare_sequence_array['U_array']['UUC']['RSCU'])
				  output_line =output_line+"\t"+"\t"+"UCC"+"\t"+str(compare_sequence_array['U_array']['UCC']['number'])+"\t"+str(compare_sequence_array['U_array']['UCC']['RSCU'])
				  output_line =output_line+"\t"+"\t"+"UAC"+"\t"+str(compare_sequence_array['U_array']['UAC']['number'])+"\t"+str(compare_sequence_array['U_array']['UAC']['RSCU'])
				  output_line =output_line+"\t"+"\t"+"UGC"+"\t"+str(compare_sequence_array['U_array']['UGC']['number'])+"\t"+str(compare_sequence_array['U_array']['UGC']['RSCU'])+ "\n"
				  writer.write(output_line)
				  
				  output_line = "Leu"+"\t"+"UUA"+"\t"+str(compare_sequence_array['U_array']['UUA']['number'])+"\t"+str(compare_sequence_array['U_array']['UUA']['RSCU'])
				  output_line =output_line+"\t"+"\t"+"UCA"+"\t"+str(compare_sequence_array['U_array']['UCA']['number'])+"\t"+str(compare_sequence_array['U_array']['UCA']['RSCU'])
				  output_line =output_line+"\t"+"Stop"+"\t"+"UAA"+"\t"+str(compare_sequence_array['U_array']['UAA']['number'])+"\t"+str(compare_sequence_array['U_array']['UAA']['RSCU'])
				  output_line =output_line+"\t"+"Stop"+"\t"+"UGA"+"\t"+str(compare_sequence_array['U_array']['UGA']['number'])+"\t"+str(compare_sequence_array['U_array']['UGA']['RSCU'])+ "\n"
				  writer.write(output_line)
				  
				  output_line = ""+"\t"+"UUG"+"\t"+str(compare_sequence_array['U_array']['UUG']['number'])+"\t"+str(compare_sequence_array['U_array']['UUG']['RSCU'])
				  output_line =output_line+"\t"+"\t"+"UCG"+"\t"+str(compare_sequence_array['U_array']['UCG']['number'])+"\t"+str(compare_sequence_array['U_array']['UCG']['RSCU'])
				  output_line =output_line+"\t"+"\t"+"UAG"+"\t"+str(compare_sequence_array['U_array']['UAG']['number'])+"\t"+str(compare_sequence_array['U_array']['UAG']['RSCU'])
				  output_line =output_line+"\t"+"Trp"+"\t"+"UGG"+"\t"+str(compare_sequence_array['U_array']['UGG']['number'])+"\t"+str(compare_sequence_array['U_array']['UGG']['RSCU'])+ "\n"
				  writer.write(output_line)
				  
				  output_line = ""+"\t"+"CUU"+"\t"+str(compare_sequence_array['C_array']['CUU']['number'])+"\t"+str(compare_sequence_array['C_array']['CUU']['RSCU'])
				  output_line =output_line+"\t"+"Pro"+"\t"+"CCU"+"\t"+str(compare_sequence_array['C_array']['CCU']['number'])+"\t"+str(compare_sequence_array['C_array']['CCU']['RSCU'])
				  output_line =output_line+"\t"+"His"+"\t"+"CAU"+"\t"+str(compare_sequence_array['C_array']['CAU']['number'])+"\t"+str(compare_sequence_array['C_array']['CAU']['RSCU'])
				  output_line =output_line+"\t"+"Arg"+"\t"+"CGU"+"\t"+str(compare_sequence_array['C_array']['CGU']['number'])+"\t"+str(compare_sequence_array['C_array']['CGU']['RSCU'])+ "\n"
				  writer.write(output_line)
				  
				  output_line = ""+"\t"+"CUC"+"\t"+str(compare_sequence_array['C_array']['CUC']['number'])+"\t"+str(compare_sequence_array['C_array']['CUC']['RSCU'])
				  output_line =output_line+"\t"+"\t"+"CCC"+"\t"+str(compare_sequence_array['C_array']['CCC']['number'])+"\t"+str(compare_sequence_array['C_array']['CCC']['RSCU'])
				  output_line =output_line+"\t"+"\t"+"CAC"+"\t"+str(compare_sequence_array['C_array']['CAC']['number'])+"\t"+str(compare_sequence_array['C_array']['CAC']['RSCU'])
				  output_line =output_line+"\t"+"\t"+"CGC"+"\t"+str(compare_sequence_array['C_array']['CGC']['number'])+"\t"+str(compare_sequence_array['C_array']['CGC']['RSCU'])+ "\n"
				  writer.write(output_line)
				  
				  output_line = ""+"\t"+"CUA"+"\t"+str(compare_sequence_array['C_array']['CUA']['number'])+"\t"+str(compare_sequence_array['C_array']['CUA']['RSCU'])
				  output_line =output_line+"\t"+"\t"+"CCA"+"\t"+str(compare_sequence_array['C_array']['CCA']['number'])+"\t"+str(compare_sequence_array['C_array']['CCA']['RSCU'])
				  output_line =output_line+"\t"+"Gln"+"\t"+"CAA"+"\t"+str(compare_sequence_array['C_array']['CAA']['number'])+"\t"+str(compare_sequence_array['C_array']['CAA']['RSCU'])
				  output_line =output_line+"\t"+"\t"+"CAG"+"\t"+str(compare_sequence_array['C_array']['CAG']['number'])+"\t"+str(compare_sequence_array['C_array']['CAG']['RSCU'])+ "\n"
				  writer.write(output_line)
				  
				  output_line = ""+"\t"+"CUG"+"\t"+str(compare_sequence_array['C_array']['CUG']['number'])+"\t"+str(compare_sequence_array['C_array']['CUG']['RSCU'])
				  output_line = output_line + "\t"+"\t"+"CCG"+"\t"+str(compare_sequence_array['C_array']['CCG']['number'])+"\t"+str(compare_sequence_array['C_array']['CCG']['RSCU'])
				  output_line = output_line + "\t"+"\t"+"CAG"+"\t"+str(compare_sequence_array['C_array']['CAG']['number'])+"\t"+str(compare_sequence_array['C_array']['CAG']['RSCU'])
				  output_line = output_line + "\t"+"\t"+"CGG"+"\t"+str(compare_sequence_array['C_array']['CGG']['number'])+"\t"+str(compare_sequence_array['C_array']['CGG']['RSCU'])+ "\n"
				  writer.write(output_line)
				  
				  output_line = "Iso"+"\t"+"AUU"+"\t"+str(compare_sequence_array['A_array']['AUU']['number'])+"\t"+str(compare_sequence_array['A_array']['AUU']['RSCU'])
				  output_line =output_line+"\t"+"Thr"+"\t"+"ACU"+"\t"+str(compare_sequence_array['A_array']['ACU']['number'])+"\t"+str(compare_sequence_array['A_array']['ACU']['RSCU'])
				  output_line =output_line+"\t"+"Asn"+"\t"+"AAU"+"\t"+str(compare_sequence_array['A_array']['AAU']['number'])+"\t"+str(compare_sequence_array['A_array']['AAU']['RSCU'])
				  output_line =output_line+"\t"+"Ser"+"\t"+"AGU"+"\t"+str(compare_sequence_array['A_array']['AGU']['number'])+"\t"+str(compare_sequence_array['A_array']['AGU']['RSCU'])+ "\n"
				  writer.write(output_line)
				  
				  output_line = ""+"\t"+"AUC"+"\t"+str(compare_sequence_array['A_array']['AUC']['number'])+"\t"+str(compare_sequence_array['A_array']['AUC']['RSCU'])
				  output_line =output_line+"\t"+"\t"+"ACC"+"\t"+str(compare_sequence_array['A_array']['ACC']['number'])+"\t"+str(compare_sequence_array['A_array']['ACC']['RSCU'])
				  output_line =output_line+"\t"+"\t"+"AAC"+"\t"+str(compare_sequence_array['A_array']['AAC']['number'])+"\t"+str(compare_sequence_array['A_array']['AAC']['RSCU'])
				  output_line =output_line+"\t"+"\t"+"AGC"+"\t"+str(compare_sequence_array['A_array']['AGC']['number'])+"\t"+str(compare_sequence_array['A_array']['AGC']['RSCU'])+ "\n"
				  writer.write(output_line)
				  
				  output_line = ""+"\t"+"AUA"+"\t"+str(compare_sequence_array['A_array']['AUA']['number'])+"\t"+str(compare_sequence_array['A_array']['AUA']['RSCU'])
				  output_line =output_line+"\t"+"\t"+"ACA"+"\t"+str(compare_sequence_array['A_array']['ACA']['number'])+"\t"+str(compare_sequence_array['A_array']['ACA']['RSCU'])
				  output_line =output_line+"\t"+"Lys"+"\t"+"AAA"+"\t"+str(compare_sequence_array['A_array']['AAA']['number'])+"\t"+str(compare_sequence_array['A_array']['AAA']['RSCU'])
				  output_line =output_line+"\t"+"Arg"+"\t"+"AGA"+"\t"+str(compare_sequence_array['A_array']['AGA']['number'])+"\t"+str(compare_sequence_array['A_array']['AGA']['RSCU'])+ "\n"
				  writer.write(output_line)
				  
				  output_line = "Met"+"\t"+"AUG"+"\t"+str(compare_sequence_array['A_array']['AUG']['number'])+"\t"+str(compare_sequence_array['A_array']['AUG']['RSCU'])
				  output_line =output_line+"\t"+"\t"+"ACG"+"\t"+str(compare_sequence_array['A_array']['ACG']['number'])+"\t"+str(compare_sequence_array['A_array']['ACG']['RSCU'])
				  output_line =output_line+"\t"+"\t"+"AAG"+"\t"+str(compare_sequence_array['A_array']['AAG']['number'])+"\t"+str(compare_sequence_array['A_array']['AAG']['RSCU'])
				  output_line =output_line+"\t"+"\t"+"AGG"+"\t"+str(compare_sequence_array['A_array']['AGG']['number'])+"\t"+str(compare_sequence_array['A_array']['AGG']['RSCU'])+ "\n"
				  writer.write(output_line)
				  
				  output_line = "Val"+"\t"+"GUU"+"\t"+str(compare_sequence_array['G_array']['GUU']['number'])+"\t"+str(compare_sequence_array['G_array']['GUU']['RSCU'])
				  output_line =output_line+"\t"+"Ala"+"\t"+"GCU"+"\t"+str(compare_sequence_array['G_array']['GCU']['number'])+"\t"+str(compare_sequence_array['G_array']['GCU']['RSCU'])
				  output_line =output_line+"\t"+"Asp"+"\t"+"GAU"+"\t"+str(compare_sequence_array['G_array']['GAU']['number'])+"\t"+str(compare_sequence_array['G_array']['GAU']['RSCU'])
				  output_line =output_line+"\t"+"Gly"+"\t"+"GGU"+"\t"+str(compare_sequence_array['G_array']['GGU']['number'])+"\t"+str(compare_sequence_array['G_array']['GGU']['RSCU'])+ "\n"
				  writer.write(output_line)
				  
				  output_line = ""+"\t"+"GUC"+"\t"+str(compare_sequence_array['G_array']['GUC']['number'])+"\t"+str(compare_sequence_array['G_array']['GUC']['RSCU'])
				  output_line =output_line+"\t"+"\t"+"GCC"+"\t"+str(compare_sequence_array['G_array']['GCC']['number'])+"\t"+str(compare_sequence_array['G_array']['GCC']['RSCU'])
				  output_line =output_line+"\t"+"\t"+"GAC"+"\t"+str(compare_sequence_array['G_array']['GAC']['number'])+"\t"+str(compare_sequence_array['G_array']['GAC']['RSCU'])
				  output_line =output_line+"\t"+"\t"+"GGC"+"\t"+str(compare_sequence_array['G_array']['GGC']['number'])+"\t"+str(compare_sequence_array['G_array']['GGC']['RSCU']) + "\n"
				  writer.write(output_line)
				  
				  output_line = ""+"\t"+"GUA"+"\t"+str(compare_sequence_array['G_array']['GUA']['number'])+"\t"+str(compare_sequence_array['G_array']['GUA']['RSCU'])
				  output_line =output_line+"\t"+"\t"+"GCA"+"\t"+str(compare_sequence_array['G_array']['GCA']['number'])+"\t"+str(compare_sequence_array['G_array']['GCA']['RSCU'])
				  output_line =output_line+"\t"+"Glu"+"\t"+"GAA"+"\t"+str(compare_sequence_array['G_array']['GAA']['number'])+"\t"+str(compare_sequence_array['G_array']['GAA']['RSCU'])
				  output_line =output_line+"\t"+"\t"+"GGA"+"\t"+str(compare_sequence_array['G_array']['GGA']['number'])+"\t"+str(compare_sequence_array['G_array']['GGA']['RSCU'])+ "\n"
				  writer.write(output_line)
				  
				  output_line = ""+"\t"+"GUG"+"\t"+str(compare_sequence_array['G_array']['GUG']['number'])+"\t"+str(compare_sequence_array['G_array']['GUG']['RSCU'])
				  output_line =output_line+"\t"+"\t"+"GCG"+"\t"+str(compare_sequence_array['G_array']['GCG']['number'])+"\t"+str(compare_sequence_array['G_array']['GCG']['RSCU'])
				  output_line =output_line+"\t"+"\t"+"GAG"+"\t"+str(compare_sequence_array['G_array']['GAG']['number'])+"\t"+str(compare_sequence_array['G_array']['GAG']['RSCU'])
				  output_line =output_line+"\t"+"\t"+"GGG"+"\t"+str(compare_sequence_array['G_array']['GGG']['number'])+"\t"+str(compare_sequence_array['G_array']['GGG']['RSCU'])+ "\n\n"
				  writer.write(output_line)
			  compare_sequence_name=line[1:]
			  compare_sequence_start+=1
		  elif(line[0:1]!='>' and base_sequence_start==1 and compare_sequence_start==0):
			  line=line.upper()
			  error_found=2
			  for z in range(0,len(line)):
				  test_string=line[z:]
				  ch=test_string[0:1]
				  if ch!='A' and ch!='T' and ch!='C' and ch!='G':
					error_found = 1 			  
			  if error_found ==1:
				  success= "Error in Input File" + str(linenum)
				  return success
			  else:
				base_sequence=base_sequence + line
				is_base_sequence_exist= "true"
		  elif(line[0:1]!='>' and compare_sequence_start>=1):
			  line=line.upper()
			  error_found=2
			  for z in range(0,len(line)):
				  test_string=line[z:]
				  ch=test_string[0:1]
				  if ch!='A' and ch!='T' and ch!='C' and ch!='G':
					error_found = 1 			  
			  if error_found ==1:
				  success= "Error in Input File"  + str(linenum)
				  return success
			  else:			  
				  compare_sequence=compare_sequence + line
				  is_compare_sequence_exist = "true"
				  is_base_sequence_exist = "false"
		  elif(line[0:1]!='>' and base_sequence_start==0):
			  success= "Error in Input File line number  " + str(linenum)
			  return success			  
	  
	  if base_sequence != "":
			  base_sequence=re.sub('T', 'U', base_sequence)
			  base_sequence=re.sub(' ', '', base_sequence)
			  base_sequence=base_sequence.upper()			  
			  base_sequence_length=math.floor((len(base_sequence)/3)*3)
			  base_sequence=base_sequence[0:int(base_sequence_length)]
			  base_sequence_start+=1;
			  base_array={
						'U_array':{
									'UUU':{'number' : 0,'RNA' : 'Phe','RSCU' : 0.00},
									'UUC':{'number' : 0,'RNA' : 'Phe','RSCU' : 0.00},
									'UUA':{'number' : 0,'RNA' : 'Leu','RSCU' : 0.00},
									'UUG':{'number' : 0,'RNA' : 'Leu','RSCU' : 0.00},
									'UCU':{'number' : 0,'RNA' : 'Ser','RSCU' : 0.00},
									'UCC':{'number' : 0,'RNA' : 'Ser','RSCU' : 0.00},
									'UCA':{'number' : 0,'RNA' : 'Ser','RSCU' : 0.00},
									'UCG':{'number' : 0,'RNA' : 'Ser','RSCU' : 0.00},
									'UAU':{'number' : 0,'RNA' : 'Tyr','RSCU' : 0.00},
									'UAC':{'number' : 0,'RNA' : 'Tyr','RSCU' : 0.00},
									'UAA':{'number' : 0,'RNA' : 'Stop','RSCU' : 0.00},
									'UAG':{'number' : 0,'RNA' : 'Stop','RSCU' : 0.00},
									'UGU':{'number' : 0,'RNA' : 'Cys','RSCU' : 0.00},
									'UGC':{'number' : 0,'RNA' : 'Cys','RSCU' : 0.00},
									'UGA':{'number' : 0,'RNA' : 'Stop','RSCU' : 0.00},
									'UGG':{'number' : 0,'RNA' : 'Trp','RSCU' : 0.00}
									},
						'G_array':{
									'GUU':{'number' : 0,'RNA' : 'Val','RSCU' :  0.00}, 
									'GUC':{'number' : 0,'RNA' : 'Val','RSCU' :  0.00}, 
									'GUA':{'number' : 0,'RNA' : 'Val','RSCU' :  0.00}, 
									'GUG':{'number' : 0,'RNA' : 'Val','RSCU' :  0.00}, 
									'GCU':{'number' : 0,'RNA' : 'Ala','RSCU' :  0.00}, 
									'GCC':{'number' : 0,'RNA' : 'Ala','RSCU' :  0.00}, 
									'GCA':{'number' : 0,'RNA' : 'Ala','RSCU' :  0.00}, 
									'GCG':{'number' : 0,'RNA' : 'Ala','RSCU' :  0.00}, 
									'GAU':{'number' : 0,'RNA' : 'Asp','RSCU' :  0.00}, 
									'GAC':{'number' : 0,'RNA' : 'Asp','RSCU' :  0.00}, 
									'GAA':{'number' : 0,'RNA' : 'Glu','RSCU' :  0.00}, 
									'GAG':{'number' : 0,'RNA' : 'Glu','RSCU' :  0.00}, 
									'GGU':{'number' : 0,'RNA' : 'Gly','RSCU' :  0.00}, 
									'GGC':{'number' : 0,'RNA' : 'Gly','RSCU' :  0.00}, 
									'GGA':{'number' : 0,'RNA' : 'Gly','RSCU' :  0.00}, 
									'GGG':{'number' : 0,'RNA' : 'Gly','RSCU' :  0.00}
									},
						'C_array':{
									'CUU':{'number' : 0,'RNA' : 'Leu','RSCU' :  0.00},
									'CUC':{'number' : 0,'RNA' : 'Leu','RSCU' :  0.00},
									'CUA':{'number' : 0,'RNA' : 'Leu','RSCU' :  0.00},
									'CUG':{'number' : 0,'RNA' : 'Leu','RSCU' :  0.00},
									'CCU':{'number' : 0,'RNA' : 'Pro','RSCU' :  0.00},
									'CCC':{'number' : 0,'RNA' : 'Pro','RSCU' :  0.00},
									'CCA':{'number' : 0,'RNA' : 'Pro','RSCU' :  0.00},
									'CCG':{'number' : 0,'RNA' : 'Pro','RSCU' :  0.00},
									'CAU':{'number' : 0,'RNA' : 'His','RSCU' :  0.00},
									'CAC':{'number' : 0,'RNA' : 'His','RSCU' :  0.00},
									'CAA':{'number' : 0,'RNA' : 'Gln','RSCU' :  0.00},
									'CAG':{'number' : 0,'RNA' : 'Gln','RSCU' :  0.00},
									'CGU':{'number' : 0,'RNA' : 'Arg','RSCU' :  0.00},
									'CGC':{'number' : 0,'RNA' : 'Arg','RSCU' :  0.00},
									'CGA':{'number' : 0,'RNA' : 'Arg','RSCU' :  0.00},
									'CGG':{'number' : 0,'RNA' : 'Arg','RSCU' :  0.00}
									},
						'A_array':{
									'AUU':{'number' : 0,'RNA' : 'Ile','RSCU' :  0.00} ,
									'AUC':{'number' : 0,'RNA' : 'Ile','RSCU' :  0.00} ,
									'AUA':{'number' : 0,'RNA' : 'Ile','RSCU' :  0.00} ,
									'AUG':{'number' : 0,'RNA' : 'Met','RSCU' :  0.00} ,
									'ACU':{'number' : 0,'RNA' : 'Thr','RSCU' :  0.00} ,
									'ACC':{'number' : 0,'RNA' : 'Thr','RSCU' :  0.00} ,
									'ACA':{'number' : 0,'RNA' : 'Thr','RSCU' :  0.00} ,
									'ACG':{'number' : 0,'RNA' : 'Thr','RSCU' :  0.00} ,
									'AAU':{'number' : 0,'RNA' : 'Asn','RSCU' :  0.00} ,
									'AAC':{'number' : 0,'RNA' : 'Asn','RSCU' :  0.00} ,
									'AAA':{'number' : 0,'RNA' : 'Lys','RSCU' :  0.00} ,
									'AAG':{'number' : 0,'RNA' : 'Lys','RSCU' :  0.00} ,
									'AGU':{'number' : 0,'RNA' : 'Ser','RSCU' :  0.00} ,
									'AGC':{'number' : 0,'RNA' : 'Ser','RSCU' :  0.00} ,
									'AGA':{'number' : 0,'RNA' : 'Arg','RSCU' :  0.00} ,
									'AGG':{'number' : 0,'RNA' : 'Arg','RSCU' :  0.00} 
									},
						'Amino_acid_array':{
											'Phe':{'alt_codon_num' : 2,'total_codon' : 0} ,
											'Leu':{'alt_codon_num' : 6,'total_codon' : 0} ,
											'Ile':{'alt_codon_num' : 3,'total_codon' : 0} ,
											'Met':{'alt_codon_num' : 1,'total_codon' : 0} ,
											'Val':{'alt_codon_num' : 4,'total_codon' : 0} ,
											'Ser':{'alt_codon_num' : 6,'total_codon' : 0} ,
											'Pro':{'alt_codon_num' : 4,'total_codon' : 0} ,
											'Thr':{'alt_codon_num' : 4,'total_codon' : 0} ,
											'Ala':{'alt_codon_num' : 4,'total_codon' : 0} ,
											'Tyr':{'alt_codon_num' : 2,'total_codon' : 0} ,
											'Stop':{'alt_codon_num' : 3,'total_codon' : 0} ,
											'His':{'alt_codon_num' : 2,'total_codon' : 0} ,
											'Gln':{'alt_codon_num' : 2,'total_codon' : 0} ,
											'Asn':{'alt_codon_num' : 2,'total_codon' : 0} ,
											'Lys':{'alt_codon_num' : 2,'total_codon' : 0} ,
											'Asp':{'alt_codon_num' : 2,'total_codon' : 0} ,
											'Glu':{'alt_codon_num' : 2,'total_codon' : 0} ,
											'Cys':{'alt_codon_num' : 2,'total_codon' : 0} ,
											'Trp':{'alt_codon_num' : 1,'total_codon' : 0} ,
											'Arg':{'alt_codon_num' : 6,'total_codon' : 0} ,
											'Gly':{'alt_codon_num' : 4,'total_codon' : 0} 
											}
					
						}
			  i=0;			
			  while int(i)<base_sequence_length :
						substring=base_sequence[i:]
						codon = substring[0:3]
						first_ch=codon[0:1]
						
						if(first_ch=='U'):
							if (base_array['U_array'][codon]):
								base_array['U_array'][codon]['number']=base_array['U_array'][codon]['number']+1
								
						if(first_ch=='G'):
								if (base_array['G_array'][codon]):
									base_array['G_array'][codon]['number']=base_array['G_array'][codon]['number']+1
									
						if(first_ch=='C'):
								if (base_array['C_array'][codon]):
									base_array['C_array'][codon]['number']=base_array['C_array'][codon]['number']+1
									
						if(first_ch=='A'):
								if (base_array['A_array'][codon]):
									base_array['A_array'][codon]['number']=base_array['A_array'][codon]['number']+1
						i=i+3
						
			  base_array['Amino_acid_array']['Phe']['total_codon']=base_array['U_array']['UUU']['number']+base_array['U_array']['UUC']['number'];
			  base_array['Amino_acid_array']['Leu']['total_codon']=base_array['U_array']['UUA']['number']+base_array['U_array']['UUG']['number']+base_array['C_array']['CUU']['number']+base_array['C_array']['CUC']['number']+base_array['C_array']['CUA']['number']+base_array['C_array']['CUG']['number'];
			  base_array['Amino_acid_array']['Ile']['total_codon']=base_array['A_array']['AUU']['number']+base_array['A_array']['AUC']['number']+base_array['A_array']['AUA']['number'];
			  base_array['Amino_acid_array']['Met']['total_codon']=base_array['A_array']['AUG']['number'];
			  base_array['Amino_acid_array']['Val']['total_codon']=base_array['G_array']['GUU']['number']+base_array['G_array']['GUC']['number']+base_array['G_array']['GUA']['number']+base_array['G_array']['GUG']['number'];
			  base_array['Amino_acid_array']['Ser']['total_codon']=base_array['U_array']['UCU']['number']+base_array['U_array']['UCC']['number']+base_array['U_array']['UCA']['number']+base_array['U_array']['UCG']['number']+base_array['A_array']['AGU']['number']+base_array['A_array']['AGC']['number'];
			  base_array['Amino_acid_array']['Pro']['total_codon']=base_array['C_array']['CCU']['number']+base_array['C_array']['CCC']['number']+base_array['C_array']['CCA']['number']+base_array['C_array']['CCG']['number'];
			  base_array['Amino_acid_array']['Thr']['total_codon']=base_array['A_array']['ACU']['number']+base_array['A_array']['ACC']['number']+base_array['A_array']['ACA']['number']+base_array['A_array']['ACG']['number'];
			  base_array['Amino_acid_array']['Ala']['total_codon']=base_array['G_array']['GCU']['number']+base_array['G_array']['GCC']['number']+base_array['G_array']['GCA']['number']+base_array['G_array']['GCG']['number'];
			  base_array['Amino_acid_array']['Tyr']['total_codon']=base_array['U_array']['UAU']['number']+base_array['U_array']['UAC']['number'];
			  base_array['Amino_acid_array']['Stop']['total_codon']=base_array['U_array']['UAA']['number']+base_array['U_array']['UGA']['number']+base_array['U_array']['UAG']['number'];
			  base_array['Amino_acid_array']['His']['total_codon']=base_array['C_array']['CAU']['number']+base_array['C_array']['CAC']['number'];
			  base_array['Amino_acid_array']['Gln']['total_codon']=base_array['C_array']['CAA']['number']+base_array['C_array']['CAG']['number'];
			  base_array['Amino_acid_array']['Asn']['total_codon']=base_array['A_array']['AAU']['number']+base_array['A_array']['AAC']['number'];
			  base_array['Amino_acid_array']['Lys']['total_codon']=base_array['A_array']['AAA']['number']+base_array['A_array']['AAG']['number'];
			  base_array['Amino_acid_array']['Asp']['total_codon']=base_array['G_array']['GAU']['number']+base_array['G_array']['GAC']['number'];
			  base_array['Amino_acid_array']['Glu']['total_codon']=base_array['G_array']['GAA']['number']+base_array['G_array']['GAG']['number'];
			  base_array['Amino_acid_array']['Cys']['total_codon']=base_array['U_array']['UGU']['number']+base_array['U_array']['UGC']['number'];
			  base_array['Amino_acid_array']['Trp']['total_codon']=base_array['U_array']['UGG']['number'];
			  base_array['Amino_acid_array']['Arg']['total_codon']=base_array['C_array']['CGU']['number']+base_array['C_array']['CGC']['number']+base_array['C_array']['CGA']['number']+base_array['C_array']['CGG']['number']+base_array['A_array']['AGA']['number']+base_array['A_array']['AGG']['number'];
			  base_array['Amino_acid_array']['Gly']['total_codon']=base_array['G_array']['GGU']['number']+base_array['G_array']['GGC']['number']+base_array['G_array']['GGA']['number']+base_array['G_array']['GGG']['number'];
			  
			  i=0;			
			  while int(i)<base_sequence_length :
					substring=base_sequence[i:]
					codon = substring[0:3]
					first_ch=codon[0:1]
					
					if(first_ch=='U'):
						if (base_array['U_array'][codon]):
							amio_acid_name=base_array['U_array'][codon]['RNA'];
							if base_array['Amino_acid_array'][amio_acid_name]['total_codon']>0:							
								base_array['U_array'][codon]['RSCU']=(base_array['U_array'][codon]['number']*base_array['Amino_acid_array'][amio_acid_name]['alt_codon_num'])/base_array['Amino_acid_array'][amio_acid_name]['total_codon']
								base_array['U_array'][codon]['RSCU']=round(base_array['U_array'][codon]['RSCU'],3)
							
					if(first_ch=='G'):
						if (base_array['G_array'][codon]):
							amio_acid_name=base_array['G_array'][codon]['RNA'];
							if base_array['Amino_acid_array'][amio_acid_name]['total_codon']>0:							
								base_array['G_array'][codon]['RSCU']=(base_array['G_array'][codon]['number']*base_array['Amino_acid_array'][amio_acid_name]['alt_codon_num'])/base_array['Amino_acid_array'][amio_acid_name]['total_codon'];
								base_array['G_array'][codon]['RSCU']=round(base_array['G_array'][codon]['RSCU'],3)
							
					if(first_ch=='C'):
						if (base_array['C_array'][codon]):
							amio_acid_name=base_array['C_array'][codon]['RNA'];
							if base_array['Amino_acid_array'][amio_acid_name]['total_codon']>0:
								base_array['C_array'][codon]['RSCU']=(base_array['C_array'][codon]['number']*base_array['Amino_acid_array'][amio_acid_name]['alt_codon_num'])/base_array['Amino_acid_array'][amio_acid_name]['total_codon'];
								base_array['C_array'][codon]['RSCU']=round(base_array['C_array'][codon]['RSCU'],3)
							
					if(first_ch=='A'):
						if (base_array['A_array'][codon]):
							amio_acid_name=base_array['A_array'][codon]['RNA'];
							if base_array['Amino_acid_array'][amio_acid_name]['total_codon']>0:							
								base_array['A_array'][codon]['RSCU']=(base_array['A_array'][codon]['number']*base_array['Amino_acid_array'][amio_acid_name]['alt_codon_num'])/base_array['Amino_acid_array'][amio_acid_name]['total_codon'];
								base_array['A_array'][codon]['RSCU']=round(base_array['A_array'][codon]['RSCU'],3)
					i=i+3		  
			  writer = open(ofile_path, 'wb')
			  output_line = "Sequence name :  " + base_sequence_name + "\n\n"
			  writer.write(output_line)
			  
			  output_line = "Phe"+"\t"+"UUU"+"\t"+ str(base_array['U_array']['UUU']['number'])+"\t"+ str(base_array['U_array']['UUU']['RSCU'])
			  output_line =output_line+"\t"+"Ser"+"\t"+"UCU"+"\t"+str(base_array['U_array']['UCU']['number'])+"\t"+str(base_array['U_array']['UCU']['RSCU'])
			  output_line =output_line+"\t"+"Tyr"+"\t"+"UAU"+"\t"+str(base_array['U_array']['UAU']['number'])+"\t"+str(base_array['U_array']['UAU']['RSCU'])
			  output_line =output_line+"\t"+"Cys"+"\t"+"UGU"+"\t"+str(base_array['U_array']['UGU']['number'])+"\t"+str(base_array['U_array']['UGU']['RSCU'])+ "\n"
			  writer.write(output_line)
			  
			  output_line = ""+"\t"+"UUC"+"\t"+str(base_array['U_array']['UUC']['number'])+"\t"+str(base_array['U_array']['UUC']['RSCU'])
			  output_line =output_line+"\t"+"\t"+"UCC"+"\t"+str(base_array['U_array']['UCC']['number'])+"\t"+str(base_array['U_array']['UCC']['RSCU'])
			  output_line =output_line+"\t"+"\t"+"UAC"+"\t"+str(base_array['U_array']['UAC']['number'])+"\t"+str(base_array['U_array']['UAC']['RSCU'])
			  output_line =output_line+"\t"+"\t"+"UGC"+"\t"+str(base_array['U_array']['UGC']['number'])+"\t"+str(base_array['U_array']['UGC']['RSCU'])+ "\n"
			  writer.write(output_line)
			  
			  output_line = "Leu"+"\t"+"UUA"+"\t"+str(base_array['U_array']['UUA']['number'])+"\t"+str(base_array['U_array']['UUA']['RSCU'])
			  output_line =output_line+"\t"+"\t"+"UCA"+"\t"+str(base_array['U_array']['UCA']['number'])+"\t"+str(base_array['U_array']['UCA']['RSCU'])
			  output_line =output_line+"\t"+"Stop"+"\t"+"UAA"+"\t"+str(base_array['U_array']['UAA']['number'])+"\t"+str(base_array['U_array']['UAA']['RSCU'])
			  output_line =output_line+"\t"+"Stop"+"\t"+"UGA"+"\t"+str(base_array['U_array']['UGA']['number'])+"\t"+str(base_array['U_array']['UGA']['RSCU'])+ "\n"
			  writer.write(output_line)
			  
			  output_line = ""+"\t"+"UUG"+"\t"+str(base_array['U_array']['UUG']['number'])+"\t"+str(base_array['U_array']['UUG']['RSCU'])
			  output_line =output_line+"\t"+"\t"+"UCG"+"\t"+str(base_array['U_array']['UCG']['number'])+"\t"+str(base_array['U_array']['UCG']['RSCU'])
			  output_line =output_line+"\t"+"\t"+"UAG"+"\t"+str(base_array['U_array']['UAG']['number'])+"\t"+str(base_array['U_array']['UAG']['RSCU'])
			  output_line =output_line+"\t"+"Trp"+"\t"+"UGG"+"\t"+str(base_array['U_array']['UGG']['number'])+"\t"+str(base_array['U_array']['UGG']['RSCU'])+ "\n"
			  writer.write(output_line)
			  
			  output_line = ""+"\t"+"CUU"+"\t"+str(base_array['C_array']['CUU']['number'])+"\t"+str(base_array['C_array']['CUU']['RSCU'])
			  output_line =output_line+"\t"+"Pro"+"\t"+"CCU"+"\t"+str(base_array['C_array']['CCU']['number'])+"\t"+str(base_array['C_array']['CCU']['RSCU'])
			  output_line =output_line+"\t"+"His"+"\t"+"CAU"+"\t"+str(base_array['C_array']['CAU']['number'])+"\t"+str(base_array['C_array']['CAU']['RSCU'])
			  output_line =output_line+"\t"+"Arg"+"\t"+"CGU"+"\t"+str(base_array['C_array']['CGU']['number'])+"\t"+str(base_array['C_array']['CGU']['RSCU'])+ "\n"
			  writer.write(output_line)
			  
			  output_line = ""+"\t"+"CUC"+"\t"+str(base_array['C_array']['CUC']['number'])+"\t"+str(base_array['C_array']['CUC']['RSCU'])
			  output_line =output_line+"\t"+"\t"+"CCC"+"\t"+str(base_array['C_array']['CCC']['number'])+"\t"+str(base_array['C_array']['CCC']['RSCU'])
			  output_line =output_line+"\t"+"\t"+"CAC"+"\t"+str(base_array['C_array']['CAC']['number'])+"\t"+str(base_array['C_array']['CAC']['RSCU'])
			  output_line =output_line+"\t"+"\t"+"CGC"+"\t"+str(base_array['C_array']['CGC']['number'])+"\t"+str(base_array['C_array']['CGC']['RSCU'])+ "\n"
			  writer.write(output_line)
			  
			  output_line = ""+"\t"+"CUA"+"\t"+str(base_array['C_array']['CUA']['number'])+"\t"+str(base_array['C_array']['CUA']['RSCU'])
			  output_line =output_line+"\t"+"\t"+"CCA"+"\t"+str(base_array['C_array']['CCA']['number'])+"\t"+str(base_array['C_array']['CCA']['RSCU'])
			  output_line =output_line+"\t"+"Gln"+"\t"+"CAA"+"\t"+str(base_array['C_array']['CAA']['number'])+"\t"+str(base_array['C_array']['CAA']['RSCU'])
			  output_line =output_line+"\t"+"\t"+"CAG"+"\t"+str(base_array['C_array']['CAG']['number'])+"\t"+str(base_array['C_array']['CAG']['RSCU'])+ "\n"
			  writer.write(output_line)
			  
			  output_line = ""+"\t"+"CUG"+"\t"+str(base_array['C_array']['CUG']['number'])+"\t"+str(base_array['C_array']['CUG']['RSCU'])
			  output_line = output_line + "\t"+"\t"+"CCG"+"\t"+str(base_array['C_array']['CCG']['number'])+"\t"+str(base_array['C_array']['CCG']['RSCU'])
			  output_line = output_line + "\t"+"\t"+"CAG"+"\t"+str(base_array['C_array']['CAG']['number'])+"\t"+str(base_array['C_array']['CAG']['RSCU'])
			  output_line = output_line + "\t"+"\t"+"CGG"+"\t"+str(base_array['C_array']['CGG']['number'])+"\t"+str(base_array['C_array']['CGG']['RSCU'])+ "\n"
			  writer.write(output_line)
			  
			  output_line = "Iso"+"\t"+"AUU"+"\t"+str(base_array['A_array']['AUU']['number'])+"\t"+str(base_array['A_array']['AUU']['RSCU'])
			  output_line =output_line+"\t"+"Thr"+"\t"+"ACU"+"\t"+str(base_array['A_array']['ACU']['number'])+"\t"+str(base_array['A_array']['ACU']['RSCU'])
			  output_line =output_line+"\t"+"Asn"+"\t"+"AAU"+"\t"+str(base_array['A_array']['AAU']['number'])+"\t"+str(base_array['A_array']['AAU']['RSCU'])
			  output_line =output_line+"\t"+"Ser"+"\t"+"AGU"+"\t"+str(base_array['A_array']['AGU']['number'])+"\t"+str(base_array['A_array']['AGU']['RSCU'])+ "\n"
			  writer.write(output_line)
			  
			  output_line = ""+"\t"+"AUC"+"\t"+str(base_array['A_array']['AUC']['number'])+"\t"+str(base_array['A_array']['AUC']['RSCU'])
			  output_line =output_line+"\t"+"\t"+"ACC"+"\t"+str(base_array['A_array']['ACC']['number'])+"\t"+str(base_array['A_array']['ACC']['RSCU'])
			  output_line =output_line+"\t"+"\t"+"AAC"+"\t"+str(base_array['A_array']['AAC']['number'])+"\t"+str(base_array['A_array']['AAC']['RSCU'])
			  output_line =output_line+"\t"+"\t"+"AGC"+"\t"+str(base_array['A_array']['AGC']['number'])+"\t"+str(base_array['A_array']['AGC']['RSCU'])+ "\n"
			  writer.write(output_line)
			  
			  output_line = ""+"\t"+"AUA"+"\t"+str(base_array['A_array']['AUA']['number'])+"\t"+str(base_array['A_array']['AUA']['RSCU'])
			  output_line =output_line+"\t"+"\t"+"ACA"+"\t"+str(base_array['A_array']['ACA']['number'])+"\t"+str(base_array['A_array']['ACA']['RSCU'])
			  output_line =output_line+"\t"+"Lys"+"\t"+"AAA"+"\t"+str(base_array['A_array']['AAA']['number'])+"\t"+str(base_array['A_array']['AAA']['RSCU'])
			  output_line =output_line+"\t"+"Arg"+"\t"+"AGA"+"\t"+str(base_array['A_array']['AGA']['number'])+"\t"+str(base_array['A_array']['AGA']['RSCU'])+ "\n"
			  writer.write(output_line)
			  
			  output_line = "Met"+"\t"+"AUG"+"\t"+str(base_array['A_array']['AUG']['number'])+"\t"+str(base_array['A_array']['AUG']['RSCU'])
			  output_line =output_line+"\t"+"\t"+"ACG"+"\t"+str(base_array['A_array']['ACG']['number'])+"\t"+str(base_array['A_array']['ACG']['RSCU'])
			  output_line =output_line+"\t"+"\t"+"AAG"+"\t"+str(base_array['A_array']['AAG']['number'])+"\t"+str(base_array['A_array']['AAG']['RSCU'])
			  output_line =output_line+"\t"+"\t"+"AGG"+"\t"+str(base_array['A_array']['AGG']['number'])+"\t"+str(base_array['A_array']['AGG']['RSCU'])+ "\n"
			  writer.write(output_line)
			  
			  output_line = "Val"+"\t"+"GUU"+"\t"+str(base_array['G_array']['GUU']['number'])+"\t"+str(base_array['G_array']['GUU']['RSCU'])
			  output_line =output_line+"\t"+"Ala"+"\t"+"GCU"+"\t"+str(base_array['G_array']['GCU']['number'])+"\t"+str(base_array['G_array']['GCU']['RSCU'])
			  output_line =output_line+"\t"+"Asp"+"\t"+"GAU"+"\t"+str(base_array['G_array']['GAU']['number'])+"\t"+str(base_array['G_array']['GAU']['RSCU'])
			  output_line =output_line+"\t"+"Gly"+"\t"+"GGU"+"\t"+str(base_array['G_array']['GGU']['number'])+"\t"+str(base_array['G_array']['GGU']['RSCU'])+ "\n"
			  writer.write(output_line)
			  
			  output_line = ""+"\t"+"GUC"+"\t"+str(base_array['G_array']['GUC']['number'])+"\t"+str(base_array['G_array']['GUC']['RSCU'])
			  output_line =output_line+"\t"+"\t"+"GCC"+"\t"+str(base_array['G_array']['GCC']['number'])+"\t"+str(base_array['G_array']['GCC']['RSCU'])
			  output_line =output_line+"\t"+"\t"+"GAC"+"\t"+str(base_array['G_array']['GAC']['number'])+"\t"+str(base_array['G_array']['GAC']['RSCU'])
			  output_line =output_line+"\t"+"\t"+"GGC"+"\t"+str(base_array['G_array']['GGC']['number'])+"\t"+str(base_array['G_array']['GGC']['RSCU'])+ "\n"
			  writer.write(output_line)
			  
			  output_line = ""+"\t"+"GUA"+"\t"+str(base_array['G_array']['GUA']['number'])+"\t"+str(base_array['G_array']['GUA']['RSCU'])
			  output_line =output_line+"\t"+"\t"+"GCA"+"\t"+str(base_array['G_array']['GCA']['number'])+"\t"+str(base_array['G_array']['GCA']['RSCU'])
			  output_line =output_line+"\t"+"Glu"+"\t"+"GAA"+"\t"+str(base_array['G_array']['GAA']['number'])+"\t"+str(base_array['G_array']['GAA']['RSCU'])
			  output_line =output_line+"\t"+"\t"+"GGA"+"\t"+str(base_array['G_array']['GGA']['number'])+"\t"+str(base_array['G_array']['GGA']['RSCU'])+ "\n"
			  writer.write(output_line)
			  
			  output_line = ""+"\t"+"GUG"+"\t"+str(base_array['G_array']['GUG']['number'])+"\t"+str(base_array['G_array']['GUG']['RSCU'])
			  output_line =output_line+"\t"+"\t"+"GCG"+"\t"+str(base_array['G_array']['GCG']['number'])+"\t"+str(base_array['G_array']['GCG']['RSCU'])
			  output_line =output_line+"\t"+"\t"+"GAG"+"\t"+str(base_array['G_array']['GAG']['number'])+"\t"+str(base_array['G_array']['GAG']['RSCU'])
			  output_line =output_line+"\t"+"\t"+"GGG"+"\t"+str(base_array['G_array']['GGG']['number'])+"\t"+str(base_array['G_array']['GGG']['RSCU'])+ "\n\n"
			  writer.write(output_line)
			  
			  writer.close()
			  success = "yes"
	  elif compare_sequence != "":	  
		  compare_sequence=re.sub('T', 'U', compare_sequence)
		  compare_sequence=re.sub(' ', '', compare_sequence)
		  compare_sequence=compare_sequence.upper()	  
		  compare_sequence_length=math.floor((len(compare_sequence)/3)*3)
		  compare_sequence=compare_sequence[0:int(compare_sequence_length)]
		  compare_sequence_array={
					'U_array':{
								'UUU':{'number' : 0,'RNA' : 'Phe','RSCU' : 0.00},
								'UUC':{'number' : 0,'RNA' : 'Phe','RSCU' : 0.00},
								'UUA':{'number' : 0,'RNA' : 'Leu','RSCU' : 0.00},
								'UUG':{'number' : 0,'RNA' : 'Leu','RSCU' : 0.00},
								'UCU':{'number' : 0,'RNA' : 'Ser','RSCU' : 0.00},
								'UCC':{'number' : 0,'RNA' : 'Ser','RSCU' : 0.00},
								'UCA':{'number' : 0,'RNA' : 'Ser','RSCU' : 0.00},
								'UCG':{'number' : 0,'RNA' : 'Ser','RSCU' : 0.00},
								'UAU':{'number' : 0,'RNA' : 'Tyr','RSCU' : 0.00},
								'UAC':{'number' : 0,'RNA' : 'Tyr','RSCU' : 0.00},
								'UAA':{'number' : 0,'RNA' : 'Stop','RSCU' : 0.00},
								'UAG':{'number' : 0,'RNA' : 'Stop','RSCU' : 0.00},
								'UGU':{'number' : 0,'RNA' : 'Cys','RSCU' : 0.00},
								'UGC':{'number' : 0,'RNA' : 'Cys','RSCU' : 0.00},
								'UGA':{'number' : 0,'RNA' : 'Stop','RSCU' : 0.00},
								'UGG':{'number' : 0,'RNA' : 'Trp','RSCU' : 0.00}
								},
					'G_array':{
								'GUU':{'number' : 0,'RNA' : 'Val','RSCU' :  0.00}, 
								'GUC':{'number' : 0,'RNA' : 'Val','RSCU' :  0.00}, 
								'GUA':{'number' : 0,'RNA' : 'Val','RSCU' :  0.00}, 
								'GUG':{'number' : 0,'RNA' : 'Val','RSCU' :  0.00}, 
								'GCU':{'number' : 0,'RNA' : 'Ala','RSCU' :  0.00}, 
								'GCC':{'number' : 0,'RNA' : 'Ala','RSCU' :  0.00}, 
								'GCA':{'number' : 0,'RNA' : 'Ala','RSCU' :  0.00}, 
								'GCG':{'number' : 0,'RNA' : 'Ala','RSCU' :  0.00}, 
								'GAU':{'number' : 0,'RNA' : 'Asp','RSCU' :  0.00}, 
								'GAC':{'number' : 0,'RNA' : 'Asp','RSCU' :  0.00}, 
								'GAA':{'number' : 0,'RNA' : 'Glu','RSCU' :  0.00}, 
								'GAG':{'number' : 0,'RNA' : 'Glu','RSCU' :  0.00}, 
								'GGU':{'number' : 0,'RNA' : 'Gly','RSCU' :  0.00}, 
								'GGC':{'number' : 0,'RNA' : 'Gly','RSCU' :  0.00}, 
								'GGA':{'number' : 0,'RNA' : 'Gly','RSCU' :  0.00}, 
								'GGG':{'number' : 0,'RNA' : 'Gly','RSCU' :  0.00}
								},
					'C_array':{
								'CUU':{'number' : 0,'RNA' : 'Leu','RSCU' :  0.00},
								'CUC':{'number' : 0,'RNA' : 'Leu','RSCU' :  0.00},
								'CUA':{'number' : 0,'RNA' : 'Leu','RSCU' :  0.00},
								'CUG':{'number' : 0,'RNA' : 'Leu','RSCU' :  0.00},
								'CCU':{'number' : 0,'RNA' : 'Pro','RSCU' :  0.00},
								'CCC':{'number' : 0,'RNA' : 'Pro','RSCU' :  0.00},
								'CCA':{'number' : 0,'RNA' : 'Pro','RSCU' :  0.00},
								'CCG':{'number' : 0,'RNA' : 'Pro','RSCU' :  0.00},
								'CAU':{'number' : 0,'RNA' : 'His','RSCU' :  0.00},
								'CAC':{'number' : 0,'RNA' : 'His','RSCU' :  0.00},
								'CAA':{'number' : 0,'RNA' : 'Gln','RSCU' :  0.00},
								'CAG':{'number' : 0,'RNA' : 'Gln','RSCU' :  0.00},
								'CGU':{'number' : 0,'RNA' : 'Arg','RSCU' :  0.00},
								'CGC':{'number' : 0,'RNA' : 'Arg','RSCU' :  0.00},
								'CGA':{'number' : 0,'RNA' : 'Arg','RSCU' :  0.00},
								'CGG':{'number' : 0,'RNA' : 'Arg','RSCU' :  0.00}
								},
					'A_array':{
								'AUU':{'number' : 0,'RNA' : 'Ile','RSCU' :  0.00} ,
								'AUC':{'number' : 0,'RNA' : 'Ile','RSCU' :  0.00} ,
								'AUA':{'number' : 0,'RNA' : 'Ile','RSCU' :  0.00} ,
								'AUG':{'number' : 0,'RNA' : 'Met','RSCU' :  0.00} ,
								'ACU':{'number' : 0,'RNA' : 'Thr','RSCU' :  0.00} ,
								'ACC':{'number' : 0,'RNA' : 'Thr','RSCU' :  0.00} ,
								'ACA':{'number' : 0,'RNA' : 'Thr','RSCU' :  0.00} ,
								'ACG':{'number' : 0,'RNA' : 'Thr','RSCU' :  0.00} ,
								'AAU':{'number' : 0,'RNA' : 'Asn','RSCU' :  0.00} ,
								'AAC':{'number' : 0,'RNA' : 'Asn','RSCU' :  0.00} ,
								'AAA':{'number' : 0,'RNA' : 'Lys','RSCU' :  0.00} ,
								'AAG':{'number' : 0,'RNA' : 'Lys','RSCU' :  0.00} ,
								'AGU':{'number' : 0,'RNA' : 'Ser','RSCU' :  0.00} ,
								'AGC':{'number' : 0,'RNA' : 'Ser','RSCU' :  0.00} ,
								'AGA':{'number' : 0,'RNA' : 'Arg','RSCU' :  0.00} ,
								'AGG':{'number' : 0,'RNA' : 'Arg','RSCU' :  0.00} 
								},
					'Amino_acid_array':{
										'Phe':{'alt_codon_num' : 2,'total_codon' : 0} ,
										'Leu':{'alt_codon_num' : 6,'total_codon' : 0} ,
										'Ile':{'alt_codon_num' : 3,'total_codon' : 0} ,
										'Met':{'alt_codon_num' : 1,'total_codon' : 0} ,
										'Val':{'alt_codon_num' : 4,'total_codon' : 0} ,
										'Ser':{'alt_codon_num' : 6,'total_codon' : 0} ,
										'Pro':{'alt_codon_num' : 4,'total_codon' : 0} ,
										'Thr':{'alt_codon_num' : 4,'total_codon' : 0} ,
										'Ala':{'alt_codon_num' : 4,'total_codon' : 0} ,
										'Tyr':{'alt_codon_num' : 2,'total_codon' : 0} ,
										'Stop':{'alt_codon_num' : 3,'total_codon' : 0} ,
										'His':{'alt_codon_num' : 2,'total_codon' : 0} ,
										'Gln':{'alt_codon_num' : 2,'total_codon' : 0} ,
										'Asn':{'alt_codon_num' : 2,'total_codon' : 0} ,
										'Lys':{'alt_codon_num' : 2,'total_codon' : 0} ,
										'Asp':{'alt_codon_num' : 2,'total_codon' : 0} ,
										'Glu':{'alt_codon_num' : 2,'total_codon' : 0} ,
										'Cys':{'alt_codon_num' : 2,'total_codon' : 0} ,
										'Trp':{'alt_codon_num' : 1,'total_codon' : 0} ,
										'Arg':{'alt_codon_num' : 6,'total_codon' : 0} ,
										'Gly':{'alt_codon_num' : 4,'total_codon' : 0} 
										}
				
					}
		  i=0;			
		  while int(i)<compare_sequence_length :
					substring=compare_sequence[i:]
					codon = substring[0:3]
					first_ch=codon[0:1]
					if(first_ch=='U'):
						if (compare_sequence_array['U_array'][codon]):
							compare_sequence_array['U_array'][codon]['number']=compare_sequence_array['U_array'][codon]['number']+1
							
					if(first_ch=='G'):
							if (compare_sequence_array['G_array'][codon]):
								compare_sequence_array['G_array'][codon]['number']=compare_sequence_array['G_array'][codon]['number']+1
								
					if(first_ch=='C'):
							if (compare_sequence_array['C_array'][codon]):
								compare_sequence_array['C_array'][codon]['number']=compare_sequence_array['C_array'][codon]['number']+1
								
					if(first_ch=='A'):
							if (compare_sequence_array['A_array'][codon]):
								compare_sequence_array['A_array'][codon]['number']=compare_sequence_array['A_array'][codon]['number']+1
					i=i+3
					
		  compare_sequence_array['Amino_acid_array']['Phe']['total_codon']=compare_sequence_array['U_array']['UUU']['number']+compare_sequence_array['U_array']['UUC']['number'];
		  compare_sequence_array['Amino_acid_array']['Leu']['total_codon']=compare_sequence_array['U_array']['UUA']['number']+compare_sequence_array['U_array']['UUG']['number']+compare_sequence_array['C_array']['CUU']['number']+compare_sequence_array['C_array']['CUC']['number']+compare_sequence_array['C_array']['CUA']['number']+compare_sequence_array['C_array']['CUG']['number'];
		  compare_sequence_array['Amino_acid_array']['Ile']['total_codon']=compare_sequence_array['A_array']['AUU']['number']+compare_sequence_array['A_array']['AUC']['number']+compare_sequence_array['A_array']['AUA']['number'];
		  compare_sequence_array['Amino_acid_array']['Met']['total_codon']=compare_sequence_array['A_array']['AUG']['number'];
		  compare_sequence_array['Amino_acid_array']['Val']['total_codon']=compare_sequence_array['G_array']['GUU']['number']+compare_sequence_array['G_array']['GUC']['number']+compare_sequence_array['G_array']['GUA']['number']+compare_sequence_array['G_array']['GUG']['number'];
		  compare_sequence_array['Amino_acid_array']['Ser']['total_codon']=compare_sequence_array['U_array']['UCU']['number']+compare_sequence_array['U_array']['UCC']['number']+compare_sequence_array['U_array']['UCA']['number']+compare_sequence_array['U_array']['UCG']['number']+compare_sequence_array['A_array']['AGU']['number']+compare_sequence_array['A_array']['AGC']['number'];
		  compare_sequence_array['Amino_acid_array']['Pro']['total_codon']=compare_sequence_array['C_array']['CCU']['number']+compare_sequence_array['C_array']['CCC']['number']+compare_sequence_array['C_array']['CCA']['number']+compare_sequence_array['C_array']['CCG']['number'];
		  compare_sequence_array['Amino_acid_array']['Thr']['total_codon']=compare_sequence_array['A_array']['ACU']['number']+compare_sequence_array['A_array']['ACC']['number']+compare_sequence_array['A_array']['ACA']['number']+compare_sequence_array['A_array']['ACG']['number'];
		  compare_sequence_array['Amino_acid_array']['Ala']['total_codon']=compare_sequence_array['G_array']['GCU']['number']+compare_sequence_array['G_array']['GCC']['number']+compare_sequence_array['G_array']['GCA']['number']+compare_sequence_array['G_array']['GCG']['number'];
		  compare_sequence_array['Amino_acid_array']['Tyr']['total_codon']=compare_sequence_array['U_array']['UAU']['number']+compare_sequence_array['U_array']['UAC']['number'];
		  compare_sequence_array['Amino_acid_array']['Stop']['total_codon']=compare_sequence_array['U_array']['UAA']['number']+compare_sequence_array['U_array']['UGA']['number']+compare_sequence_array['U_array']['UAG']['number'];
		  compare_sequence_array['Amino_acid_array']['His']['total_codon']=compare_sequence_array['C_array']['CAU']['number']+compare_sequence_array['C_array']['CAC']['number'];
		  compare_sequence_array['Amino_acid_array']['Gln']['total_codon']=compare_sequence_array['C_array']['CAA']['number']+compare_sequence_array['C_array']['CAG']['number'];
		  compare_sequence_array['Amino_acid_array']['Asn']['total_codon']=compare_sequence_array['A_array']['AAU']['number']+compare_sequence_array['A_array']['AAC']['number'];
		  compare_sequence_array['Amino_acid_array']['Lys']['total_codon']=compare_sequence_array['A_array']['AAA']['number']+compare_sequence_array['A_array']['AAG']['number'];
		  compare_sequence_array['Amino_acid_array']['Asp']['total_codon']=compare_sequence_array['G_array']['GAU']['number']+compare_sequence_array['G_array']['GAC']['number'];
		  compare_sequence_array['Amino_acid_array']['Glu']['total_codon']=compare_sequence_array['G_array']['GAA']['number']+compare_sequence_array['G_array']['GAG']['number'];
		  compare_sequence_array['Amino_acid_array']['Cys']['total_codon']=compare_sequence_array['U_array']['UGU']['number']+compare_sequence_array['U_array']['UGC']['number'];
		  compare_sequence_array['Amino_acid_array']['Trp']['total_codon']=compare_sequence_array['U_array']['UGG']['number'];
		  compare_sequence_array['Amino_acid_array']['Arg']['total_codon']=compare_sequence_array['C_array']['CGU']['number']+compare_sequence_array['C_array']['CGC']['number']+compare_sequence_array['C_array']['CGA']['number']+compare_sequence_array['C_array']['CGG']['number']+compare_sequence_array['A_array']['AGA']['number']+compare_sequence_array['A_array']['AGG']['number'];
		  compare_sequence_array['Amino_acid_array']['Gly']['total_codon']=compare_sequence_array['G_array']['GGU']['number']+compare_sequence_array['G_array']['GGC']['number']+compare_sequence_array['G_array']['GGA']['number']+compare_sequence_array['G_array']['GGG']['number'];
		  
		  i=0			
		  while int(i)<compare_sequence_length :
				substring=compare_sequence[i:]
				codon = substring[0:3]
				first_ch=codon[0:1]
				
				if(first_ch=='U'):
					if (compare_sequence_array['U_array'][codon]):
						amio_acid_name=compare_sequence_array['U_array'][codon]['RNA'];
						if compare_sequence_array['Amino_acid_array'][amio_acid_name]['total_codon']>0:						
							compare_sequence_array['U_array'][codon]['RSCU']=(compare_sequence_array['U_array'][codon]['number']*compare_sequence_array['Amino_acid_array'][amio_acid_name]['alt_codon_num'])/compare_sequence_array['Amino_acid_array'][amio_acid_name]['total_codon']
							compare_sequence_array['U_array'][codon]['RSCU']=round(compare_sequence_array['U_array'][codon]['RSCU'],3)
				if(first_ch=='G'):
					if (compare_sequence_array['G_array'][codon]):
						amio_acid_name=compare_sequence_array['G_array'][codon]['RNA'];
						if compare_sequence_array['Amino_acid_array'][amio_acid_name]['total_codon']>0:						
							compare_sequence_array['G_array'][codon]['RSCU']=(compare_sequence_array['G_array'][codon]['number']*compare_sequence_array['Amino_acid_array'][amio_acid_name]['alt_codon_num'])/compare_sequence_array['Amino_acid_array'][amio_acid_name]['total_codon'];
							compare_sequence_array['G_array'][codon]['RSCU']=round(compare_sequence_array['G_array'][codon]['RSCU'],3)
				if(first_ch=='C'):
					if (compare_sequence_array['C_array'][codon]):
						amio_acid_name=compare_sequence_array['C_array'][codon]['RNA'];
						if compare_sequence_array['Amino_acid_array'][amio_acid_name]['total_codon']>0:						
							compare_sequence_array['C_array'][codon]['RSCU']=(compare_sequence_array['C_array'][codon]['number']*compare_sequence_array['Amino_acid_array'][amio_acid_name]['alt_codon_num'])/compare_sequence_array['Amino_acid_array'][amio_acid_name]['total_codon'];
							compare_sequence_array['C_array'][codon]['RSCU']=round(compare_sequence_array['C_array'][codon]['RSCU'],3)
				if(first_ch=='A'):
					if (compare_sequence_array['A_array'][codon]):
						amio_acid_name=compare_sequence_array['A_array'][codon]['RNA'];
						if compare_sequence_array['Amino_acid_array'][amio_acid_name]['total_codon']>0:						
							compare_sequence_array['A_array'][codon]['RSCU']=(compare_sequence_array['A_array'][codon]['number']*compare_sequence_array['Amino_acid_array'][amio_acid_name]['alt_codon_num'])/compare_sequence_array['Amino_acid_array'][amio_acid_name]['total_codon'];
							compare_sequence_array['A_array'][codon]['RSCU']=round(compare_sequence_array['A_array'][codon]['RSCU'],3)
				i=i+3
		  output_line = "Sequence name :  " + compare_sequence_name + "\n\n"
		  writer.write(output_line)
		  
		  output_line = "Phe"+"\t"+"UUU"+"\t"+str(compare_sequence_array['U_array']['UUU']['number'])+"\t"+str(compare_sequence_array['U_array']['UUU']['RSCU'])
		  output_line =output_line+"\t"+"Ser"+"\t"+"UCU"+"\t"+str(compare_sequence_array['U_array']['UCU']['number'])+"\t"+str(compare_sequence_array['U_array']['UCU']['RSCU'])
		  output_line =output_line+"\t"+"Tyr"+"\t"+"UAU"+"\t"+str(compare_sequence_array['U_array']['UAU']['number'])+"\t"+str(compare_sequence_array['U_array']['UAU']['RSCU'])
		  output_line =output_line+"\t"+"Cys"+"\t"+"UGU"+"\t"+str(compare_sequence_array['U_array']['UGU']['number'])+"\t"+str(compare_sequence_array['U_array']['UGU']['RSCU'])+ "\n"
		  writer.write(output_line)
		  
		  output_line = ""+"\t"+"UUC"+"\t"+str(compare_sequence_array['U_array']['UUC']['number'])+"\t"+str(compare_sequence_array['U_array']['UUC']['RSCU'])
		  output_line =output_line+"\t"+"\t"+"UCC"+"\t"+str(compare_sequence_array['U_array']['UCC']['number'])+"\t"+str(compare_sequence_array['U_array']['UCC']['RSCU'])
		  output_line =output_line+"\t"+"\t"+"UAC"+"\t"+str(compare_sequence_array['U_array']['UAC']['number'])+"\t"+str(compare_sequence_array['U_array']['UAC']['RSCU'])
		  output_line =output_line+"\t"+"\t"+"UGC"+"\t"+str(compare_sequence_array['U_array']['UGC']['number'])+"\t"+str(compare_sequence_array['U_array']['UGC']['RSCU'])+ "\n"
		  writer.write(output_line)
		  
		  output_line = "Leu"+"\t"+"UUA"+"\t"+str(compare_sequence_array['U_array']['UUA']['number'])+"\t"+str(compare_sequence_array['U_array']['UUA']['RSCU'])
		  output_line =output_line+"\t"+"\t"+"UCA"+"\t"+str(compare_sequence_array['U_array']['UCA']['number'])+"\t"+str(compare_sequence_array['U_array']['UCA']['RSCU'])
		  output_line =output_line+"\t"+"Stop"+"\t"+"UAA"+"\t"+str(compare_sequence_array['U_array']['UAA']['number'])+"\t"+str(compare_sequence_array['U_array']['UAA']['RSCU'])
		  output_line =output_line+"\t"+"Stop"+"\t"+"UGA"+"\t"+str(compare_sequence_array['U_array']['UGA']['number'])+"\t"+str(compare_sequence_array['U_array']['UGA']['RSCU'])+ "\n"
		  writer.write(output_line)
		  
		  output_line = ""+"\t"+"UUG"+"\t"+str(compare_sequence_array['U_array']['UUG']['number'])+"\t"+str(compare_sequence_array['U_array']['UUG']['RSCU'])
		  output_line =output_line+"\t"+"\t"+"UCG"+"\t"+str(compare_sequence_array['U_array']['UCG']['number'])+"\t"+str(compare_sequence_array['U_array']['UCG']['RSCU'])
		  output_line =output_line+"\t"+"\t"+"UAG"+"\t"+str(compare_sequence_array['U_array']['UAG']['number'])+"\t"+str(compare_sequence_array['U_array']['UAG']['RSCU'])
		  output_line =output_line+"\t"+"Trp"+"\t"+"UGG"+"\t"+str(compare_sequence_array['U_array']['UGG']['number'])+"\t"+str(compare_sequence_array['U_array']['UGG']['RSCU'])+ "\n"
		  writer.write(output_line)
		  
		  output_line = ""+"\t"+"CUU"+"\t"+str(compare_sequence_array['C_array']['CUU']['number'])+"\t"+str(compare_sequence_array['C_array']['CUU']['RSCU'])
		  output_line =output_line+"\t"+"Pro"+"\t"+"CCU"+"\t"+str(compare_sequence_array['C_array']['CCU']['number'])+"\t"+str(compare_sequence_array['C_array']['CCU']['RSCU'])
		  output_line =output_line+"\t"+"His"+"\t"+"CAU"+"\t"+str(compare_sequence_array['C_array']['CAU']['number'])+"\t"+str(compare_sequence_array['C_array']['CAU']['RSCU'])
		  output_line =output_line+"\t"+"Arg"+"\t"+"CGU"+"\t"+str(compare_sequence_array['C_array']['CGU']['number'])+"\t"+str(compare_sequence_array['C_array']['CGU']['RSCU'])+ "\n"
		  writer.write(output_line)
		  
		  output_line = ""+"\t"+"CUC"+"\t"+str(compare_sequence_array['C_array']['CUC']['number'])+"\t"+str(compare_sequence_array['C_array']['CUC']['RSCU'])
		  output_line =output_line+"\t"+"\t"+"CCC"+"\t"+str(compare_sequence_array['C_array']['CCC']['number'])+"\t"+str(compare_sequence_array['C_array']['CCC']['RSCU'])
		  output_line =output_line+"\t"+"\t"+"CAC"+"\t"+str(compare_sequence_array['C_array']['CAC']['number'])+"\t"+str(compare_sequence_array['C_array']['CAC']['RSCU'])
		  output_line =output_line+"\t"+"\t"+"CGC"+"\t"+str(compare_sequence_array['C_array']['CGC']['number'])+"\t"+str(compare_sequence_array['C_array']['CGC']['RSCU'])+ "\n"
		  writer.write(output_line)
		  
		  output_line = ""+"\t"+"CUA"+"\t"+str(compare_sequence_array['C_array']['CUA']['number'])+"\t"+str(compare_sequence_array['C_array']['CUA']['RSCU'])
		  output_line =output_line+"\t"+"\t"+"CCA"+"\t"+str(compare_sequence_array['C_array']['CCA']['number'])+"\t"+str(compare_sequence_array['C_array']['CCA']['RSCU'])
		  output_line =output_line+"\t"+"Gln"+"\t"+"CAA"+"\t"+str(compare_sequence_array['C_array']['CAA']['number'])+"\t"+str(compare_sequence_array['C_array']['CAA']['RSCU'])
		  output_line =output_line+"\t"+"\t"+"CAG"+"\t"+str(compare_sequence_array['C_array']['CAG']['number'])+"\t"+str(compare_sequence_array['C_array']['CAG']['RSCU'])+ "\n"
		  writer.write(output_line)
		  
		  output_line = ""+"\t"+"CUG"+"\t"+str(compare_sequence_array['C_array']['CUG']['number'])+"\t"+str(compare_sequence_array['C_array']['CUG']['RSCU'])
		  output_line = output_line + "\t"+"\t"+"CCG"+"\t"+str(compare_sequence_array['C_array']['CCG']['number'])+"\t"+str(compare_sequence_array['C_array']['CCG']['RSCU'])
		  output_line = output_line + "\t"+"\t"+"CAG"+"\t"+str(compare_sequence_array['C_array']['CAG']['number'])+"\t"+str(compare_sequence_array['C_array']['CAG']['RSCU'])
		  output_line = output_line + "\t"+"\t"+"CGG"+"\t"+str(compare_sequence_array['C_array']['CGG']['number'])+"\t"+str(compare_sequence_array['C_array']['CGG']['RSCU'])+ "\n"
		  writer.write(output_line)
		  
		  output_line = "Iso"+"\t"+"AUU"+"\t"+str(compare_sequence_array['A_array']['AUU']['number'])+"\t"+str(compare_sequence_array['A_array']['AUU']['RSCU'])
		  output_line =output_line+"\t"+"Thr"+"\t"+"ACU"+"\t"+str(compare_sequence_array['A_array']['ACU']['number'])+"\t"+str(compare_sequence_array['A_array']['ACU']['RSCU'])
		  output_line =output_line+"\t"+"Asn"+"\t"+"AAU"+"\t"+str(compare_sequence_array['A_array']['AAU']['number'])+"\t"+str(compare_sequence_array['A_array']['AAU']['RSCU'])
		  output_line =output_line+"\t"+"Ser"+"\t"+"AGU"+"\t"+str(compare_sequence_array['A_array']['AGU']['number'])+"\t"+str(compare_sequence_array['A_array']['AGU']['RSCU'])+ "\n"
		  writer.write(output_line)
		  
		  output_line = ""+"\t"+"AUC"+"\t"+str(compare_sequence_array['A_array']['AUC']['number'])+"\t"+str(compare_sequence_array['A_array']['AUC']['RSCU'])
		  output_line =output_line+"\t"+"\t"+"ACC"+"\t"+str(compare_sequence_array['A_array']['ACC']['number'])+"\t"+str(compare_sequence_array['A_array']['ACC']['RSCU'])
		  output_line =output_line+"\t"+"\t"+"AAC"+"\t"+str(compare_sequence_array['A_array']['AAC']['number'])+"\t"+str(compare_sequence_array['A_array']['AAC']['RSCU'])
		  output_line =output_line+"\t"+"\t"+"AGC"+"\t"+str(compare_sequence_array['A_array']['AGC']['number'])+"\t"+str(compare_sequence_array['A_array']['AGC']['RSCU'])+ "\n"
		  writer.write(output_line)
		  
		  output_line = ""+"\t"+"AUA"+"\t"+str(compare_sequence_array['A_array']['AUA']['number'])+"\t"+str(compare_sequence_array['A_array']['AUA']['RSCU'])
		  output_line =output_line+"\t"+"\t"+"ACA"+"\t"+str(compare_sequence_array['A_array']['ACA']['number'])+"\t"+str(compare_sequence_array['A_array']['ACA']['RSCU'])
		  output_line =output_line+"\t"+"Lys"+"\t"+"AAA"+"\t"+str(compare_sequence_array['A_array']['AAA']['number'])+"\t"+str(compare_sequence_array['A_array']['AAA']['RSCU'])
		  output_line =output_line+"\t"+"Arg"+"\t"+"AGA"+"\t"+str(compare_sequence_array['A_array']['AGA']['number'])+"\t"+str(compare_sequence_array['A_array']['AGA']['RSCU'])+ "\n"
		  writer.write(output_line)
		  
		  output_line = "Met"+"\t"+"AUG"+"\t"+str(compare_sequence_array['A_array']['AUG']['number'])+"\t"+str(compare_sequence_array['A_array']['AUG']['RSCU'])
		  output_line =output_line+"\t"+"\t"+"ACG"+"\t"+str(compare_sequence_array['A_array']['ACG']['number'])+"\t"+str(compare_sequence_array['A_array']['ACG']['RSCU'])
		  output_line =output_line+"\t"+"\t"+"AAG"+"\t"+str(compare_sequence_array['A_array']['AAG']['number'])+"\t"+str(compare_sequence_array['A_array']['AAG']['RSCU'])
		  output_line =output_line+"\t"+"\t"+"AGG"+"\t"+str(compare_sequence_array['A_array']['AGG']['number'])+"\t"+str(compare_sequence_array['A_array']['AGG']['RSCU'])+ "\n"
		  writer.write(output_line)
		  
		  output_line = "Val"+"\t"+"GUU"+"\t"+str(compare_sequence_array['G_array']['GUU']['number'])+"\t"+str(compare_sequence_array['G_array']['GUU']['RSCU'])
		  output_line =output_line+"\t"+"Ala"+"\t"+"GCU"+"\t"+str(compare_sequence_array['G_array']['GCU']['number'])+"\t"+str(compare_sequence_array['G_array']['GCU']['RSCU'])
		  output_line =output_line+"\t"+"Asp"+"\t"+"GAU"+"\t"+str(compare_sequence_array['G_array']['GAU']['number'])+"\t"+str(compare_sequence_array['G_array']['GAU']['RSCU'])
		  output_line =output_line+"\t"+"Gly"+"\t"+"GGU"+"\t"+str(compare_sequence_array['G_array']['GGU']['number'])+"\t"+str(compare_sequence_array['G_array']['GGU']['RSCU'])+ "\n"
		  writer.write(output_line)
		  
		  output_line = ""+"\t"+"GUC"+"\t"+str(compare_sequence_array['G_array']['GUC']['number'])+"\t"+str(compare_sequence_array['G_array']['GUC']['RSCU'])
		  output_line =output_line+"\t"+"\t"+"GCC"+"\t"+str(compare_sequence_array['G_array']['GCC']['number'])+"\t"+str(compare_sequence_array['G_array']['GCC']['RSCU'])
		  output_line =output_line+"\t"+"\t"+"GAC"+"\t"+str(compare_sequence_array['G_array']['GAC']['number'])+"\t"+str(compare_sequence_array['G_array']['GAC']['RSCU'])
		  output_line =output_line+"\t"+"\t"+"GGC"+"\t"+str(compare_sequence_array['G_array']['GGC']['number'])+"\t"+str(compare_sequence_array['G_array']['GGC']['RSCU']) + "\n"
		  writer.write(output_line)
		  
		  output_line = ""+"\t"+"GUA"+"\t"+str(compare_sequence_array['G_array']['GUA']['number'])+"\t"+str(compare_sequence_array['G_array']['GUA']['RSCU'])
		  output_line =output_line+"\t"+"\t"+"GCA"+"\t"+str(compare_sequence_array['G_array']['GCA']['number'])+"\t"+str(compare_sequence_array['G_array']['GCA']['RSCU'])
		  output_line =output_line+"\t"+"Glu"+"\t"+"GAA"+"\t"+str(compare_sequence_array['G_array']['GAA']['number'])+"\t"+str(compare_sequence_array['G_array']['GAA']['RSCU'])
		  output_line =output_line+"\t"+"\t"+"GGA"+"\t"+str(compare_sequence_array['G_array']['GGA']['number'])+"\t"+str(compare_sequence_array['G_array']['GGA']['RSCU'])+ "\n"
		  writer.write(output_line)
		  
		  output_line = ""+"\t"+"GUG"+"\t"+str(compare_sequence_array['G_array']['GUG']['number'])+"\t"+str(compare_sequence_array['G_array']['GUG']['RSCU'])
		  output_line =output_line+"\t"+"\t"+"GCG"+"\t"+str(compare_sequence_array['G_array']['GCG']['number'])+"\t"+str(compare_sequence_array['G_array']['GCG']['RSCU'])
		  output_line =output_line+"\t"+"\t"+"GAG"+"\t"+str(compare_sequence_array['G_array']['GAG']['number'])+"\t"+str(compare_sequence_array['G_array']['GAG']['RSCU'])
		  output_line =output_line+"\t"+"\t"+"GGG"+"\t"+str(compare_sequence_array['G_array']['GGG']['number'])+"\t"+str(compare_sequence_array['G_array']['GGG']['RSCU'])+ "\n\n"
		  writer.write(output_line)
		  writer.close()
		  success = "yes"
	  elif base_sequence == "" and base_sequence_name != "" :
		  success = "Base sequence is not found in line number " + str(linenum)
	  elif compare_sequence == "" and compare_sequence_name != "" :
		  success = "Compare sequence is not found in line number  " + str(linenum)
	  return success
		  
