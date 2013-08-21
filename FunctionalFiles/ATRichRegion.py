"""
Tools Name : Base Analyzer Tool
Original Developer: Afsana Rahman Snigdha
Original Development Date : 01-01-2012
Version : 1.0
Purpose : 
		Input File: A text file containing one fasta sequence.
		Output File : A html file indecating the AT rich region.
"""
import pdb
import os 
import sys
import fileinput
import math
import re
import csv
from PyQt4 import QtCore, QtGui

class ATRichRegion(QtGui.QMainWindow):

  def __init__(self, parent=None):
    super(ATRichRegion, self).__init__()
  
  def ExactCheck(self,array,index,first_position,last_position,max_AT_array_counter):
	  duplicacy_check=1
	  for k in range (0, index):
		  if int(array[k,3])==first_position and int(array[k,4])==last_position:
			  duplicacy_check=0
			  #print str(max_AT_array_counter)+"\t\t"+str(first_position)+"\t\t"+str(last_position) +"\t\t"+str(duplicacy_check)+"\t\tfirst if\n"
			  #print str(first_position)+"\t\t"+str(array[k,3])+"\t\t"+str(last_position) +"\t\t"+str(array[k,4])+"\t\t"+str(duplicacy_check)+"\t\tfirst if\n"
		  if first_position>int(array[k,3]) and first_position<int(array[k,4]) and last_position<=int(array[k,4]):
			  duplicacy_check=0
			  #print str(max_AT_array_counter)+"\t\t"+str(first_position)+"\t\t"+str(array[k,3])+"\t\t"+str(last_position) +"\t\t"+str(array[k,4])+"\t\t"+str(duplicacy_check)+"\t\tsecond if\n"
			  #print str(max_AT_array_counter)+"\t\t"+str(first_position)+"\t\t"+str(last_position) +"\t\t"+str(duplicacy_check)+"\t\tsecond if\n"
	  #if duplicacy_check==1:
	  #print str(first_position)+"\t\t"+str(last_position) +"\t\t"+str(duplicacy_check)+"\n"
		#writer.write(output_line_help)
	  return duplicacy_check
  
  def GetRichRegion(self,subset_array,max_AT_array_counter,sequence,threshold):
	  """
	  writer = open("/home/afsana/Documents/help_function.txt", 'a+')
	  """
	  prev_counter = int(subset_array[max_AT_array_counter,2])-1
	  next_counter = int(subset_array[max_AT_array_counter,3])+1
	  """
	  output_line_help ="prev_counter====" + str(prev_counter) + "\n"
	  writer.write(output_line_help)
	  output_line_help ="next_counter====" + str(next_counter) + "\n"
	  writer.write(output_line_help)
	  """
	  flag_prev=1
	  flag_next=1
	  	  	  
	  if prev_counter<0 :
		  prev_counter=int(subset_array[max_AT_array_counter,2])
		  flag_prev=2
	    
	  if next_counter==int(len(sequence)):  
		next_counter=int(subset_array[max_AT_array_counter,3])
		flag_next=3
		
	  variable_sequence=str(subset_array[max_AT_array_counter,0])
	  variable_sequence_First_Position=int(subset_array[max_AT_array_counter,2])
	  variable_sequence_Last_Position=int(subset_array[max_AT_array_counter,3])
	  

	  AT=0
	  GC=0
	  temp_array= {}
	  found=0
	  while prev_counter>=0 or next_counter<int(len(sequence)):
		  """
		  output_line_help ="prev_counter====" + str(prev_counter) + "\n"
		  writer.write(output_line_help)
		  output_line_help ="next_counter====" + str(next_counter) + "\n"
		  writer.write(output_line_help)
		  output_line_help ="flag_prev=====" + str(flag_prev) + "\n"
		  writer.write(output_line_help)
		  output_line_help ="flag_next=====" + str(flag_next) + "\n"
		  writer.write(output_line_help)
		  output_line_help ="variable_sequence=====" + str(variable_sequence) + "\n"
		  writer.write(output_line_help)
		  """
		  if prev_counter>0 and flag_prev==1 :
			  prev_counter = prev_counter-1
			  
			  if prev_counter == 0:
				  flag_prev=2
			  
			  temp = sequence[prev_counter:]
			  ch=temp[0:1]
			  variable_sequence = ch + variable_sequence
			  
			  if ch=="A" or ch=="T":
				  AT=AT+1
			  if ch=="G" or ch=="C":
				  GC=GC+1
				  
			  neumator=int(AT)+int(subset_array[max_AT_array_counter,4])
			  denometor=len(variable_sequence)
			  AT_percentage=float((neumator*100)/denometor)
			  variable_sequence_First_Position=int(prev_counter)
			  variable_sequence_Last_Position=int(prev_counter)+len(variable_sequence)
			  
			  if AT_percentage>=threshold:
				  found=1
				  temp_array[0,0] = variable_sequence
				  temp_array[0,1] = len(temp_array[0,0])
				  temp_array[0,2] = str(prev_counter)
				  temp_array[0,3] = str(next_counter)
				  temp_array[0,4] = str(AT_percentage)
				  
		  if next_counter<int(len(sequence)) and flag_next==1:
			  next_counter = next_counter+1
			  
			  if next_counter == int(len(sequence))-1:
				 flag_next=3
			  
			  temp = sequence[next_counter:]
			  ch=temp[0:1]
			  variable_sequence = variable_sequence + ch
			  
			  if ch=="A" or ch=="T":
				  AT=AT+1
			  if ch=="G" or ch=="C":
				  GC=GC+1
				  
			  neumator=int(AT)+int(subset_array[max_AT_array_counter,4])
			  denometor=len(variable_sequence)
			  AT_percentage=float((neumator*100)/denometor)
			  variable_sequence_First_Position=int(subset_array[max_AT_array_counter,1])
			  variable_sequence_Last_Position=int(next_counter)	
			  			  
			  if AT_percentage>=threshold:
				  found=1
				  temp_array[0,0] = variable_sequence
				  temp_array[0,1] = str(len(temp_array[0,0]))
				  temp_array[0,2] = str(prev_counter)
				  temp_array[0,3] = str(next_counter)
				  temp_array[0,4] = str(AT_percentage)
				  
		  if flag_prev==2 and flag_next==3:
			  break;
	  
	  if found ==0:
		  temp_array[0,0] = variable_sequence
		  temp_array[0,1] = len(variable_sequence)
		  temp_array[0,2] = str(variable_sequence_First_Position)
		  temp_array[0,3] = str(variable_sequence_Last_Position)
		  temp_array[0,4] = "100.00"		  
	  return temp_array
	  		  			  
  def ATRichRegion(self,sourcepath,output_file,percentage):
	  sequence_start =0
	  new_character = ""
	  prev_character = ""
	  temp = ""
	  junk_character = ""
	  counter = 0
	  max_AT_length = 0
	  max_GC_length = 0
	  Temp_AT_Counter = 0
	  Temp_GC_Counter = 0
	  Junk_AT_Counter = 0
	  Junk_GC_Counter = 0	
	  sequence = ""
	  
	  subset_array= {}	    
	  Rich_Region_Index_Array = []
	  result_array= {}
	  final_array = {}
	  
	  result_sequence="";
	  result_sequence_length=0;
	  result_sequence_AT_percentage=0.00
	  result_sequence_First_Position=0
	  result_sequence_Last_Position=0	  
	  linenum=0
	  for line in open(sourcepath,'r'):
		  line=line.rstrip('\n') 
		  linenum=linenum+1
		  if(line[0:1]=='>' and sequence_start==0):
			  sequence_start+=1;
			  sequence_name=line[1:]
			  
		  elif line[0:1]=='>' and sequence_start==1:
			  sequence_start+=1;
			  success= "Error in Input File line number " + str(linenum)
			  return success;
			  
		  elif(line[0:1]!='>' and sequence_start==1):
			  #print line
			  error_found=2
			  for z in range(0,len(line)):
				  test_string=line[z:]
				  ch=test_string[0:1]
				  if ch!='A' and ch!='T' and ch!='C' and ch!='G':
					error_found = 1 			  
			  if error_found ==1:
				  success= "Error in Input File line number " + str(linenum)
				  return success
			  else :			  
				sequence=sequence + line 
	  if sequence_start ==0:
		  success= "Error in Input File line number " + str(linenum)
		  return success;

	  sequence=re.sub(' ', '', sequence)	  
	  sequence_length=len(sequence)
	  i=0
	  while int(i)<sequence_length :
		  substring=sequence[i:]
		  new_character=substring[0:1];
		  ## #print "new_character====" + new_character
		  if i==0:
			  if new_character=='A' or new_character=='T':
				  temp=temp + str(new_character)
				  # #print "new_character====" + new_character + "====prev_character====" + prev_character
				  prev_character = new_character
				  
				  if i==int(sequence_length-1):
					 subset_array[counter,0]=temp;# sub sequence.
					 subset_array[counter,1]=int(len(temp)) # length of sub sequence.
					 subset_array[counter,2]=i-int(len(temp)) # to set first postion of this subsequence formula is current position - length of sub sequence.
					 subset_array[counter,3]=i;
					 subset_array[counter,4]=int(len(temp)) 
					 subset_array[counter,5]=0;
					 subset_array[counter,6]=1;
					 Rich_Region_Index_Array.append(counter)
					 #max_AT_length=subset_array[counter,1];
					 #max_AT_array_counter=counter;
					 
					 counter = counter +1  # increment counter by one to set new subsequence
			  elif new_character=='C' or new_character=='G':
				  temp=temp + str(new_character)
				  # #print "new_character====" + new_character + "====prev_character====" + prev_character
				  prev_character = new_character
				  
				  if i== int(sequence_length-1):
					 subset_array[counter,0]=temp;# sub sequence.
					 subset_array[counter,1]=len(temp); # length of sub sequence.
					 subset_array[counter,2]=i-int(len(temp)); # to set first postion of this subsequence formula is current position - length of sub sequence.
					 subset_array[counter,3]=i;
					 subset_array[counter,4]=0;
					 subset_array[counter,5]=int(len(temp))
					 subset_array[counter,6]=2;
					
					 #max_GC_length=subset_array[counter,1];
					 #max_GC_array_counter=counter;
					 counter = counter +1  # increment counter by one to set new subsequence 
		  elif i>0:
			  if (new_character=='A' or new_character=='T') and prev_character!='G' and prev_character != 'C':
				  temp=temp + str(new_character)
				  # #print "new_character====" + new_character + "====prev_character====" + prev_character
				  prev_character = new_character
				  
				  if i==int(sequence_length-1) and int(len(temp)>=4):
					  if junk_character!="":
						 subset_array[counter,0]=junk_character;# sub sequence.
						 subset_array[counter,1]=int(len(junk_character)) # length of sub sequence.
						 subset_array[counter,2]=i-int(len(temp))-int(len(junk_character))+1; # to set first postion of this subsequence formula is current position - length of sub sequence.
						 subset_array[counter,3]=i-int(len(temp))
						 subset_array[counter,4]=Junk_AT_Counter;
						 subset_array[counter,5]=Junk_GC_Counter;
						 subset_array[counter,6]=3;# set sequnce type is GC
						
						 Junk_AT_Counter=0;
						 Junk_GC_Counter=0;
						 junk_character="";
						 counter = counter +1
					  
					  subset_array[counter,0]=temp# sub sequence.
					  subset_array[counter,1]=int(len(temp)) # length of sub sequence.
					  subset_array[counter,2]=i-int(len(temp))+1 # to set first postion of this subsequence formula is current position - length of sub sequence.
					  subset_array[counter,3]=i;
					  subset_array[counter,4]=int(len(temp))
					  subset_array[counter,5]=0;
					  subset_array[counter,6]=1
					  Rich_Region_Index_Array.append(counter)
					  
					  #if subset_array[counter,1]>int(max_AT_length):
						  #max_AT_length=subset_array[counter,1]
						  #max_AT_array_counter=counter
					  counter = counter +1 
				  
				  elif i==int(sequence_length-1) and int(len(temp)<4) :
					   junk_character=junk_character +temp;
					   Junk_AT_Counter=Junk_AT_Counter+int(len(temp))
					   
					   subset_array[counter,0]=junk_character# sub sequence.
					   subset_array[counter,1]=int(len(junk_character)) # length of sub sequence.
					   subset_array[counter,2]=i-int(len(junk_character))+1; # to set first postion of this subsequence formula is current position - length of sub sequence.
					   subset_array[counter,3]=i;
					   subset_array[counter,4]=Junk_AT_Counter;
					   subset_array[counter,5]=Junk_GC_Counter;
					   subset_array[counter,6]=3;
					   counter=counter +1
					   
					   Junk_AT_Counter=0;
					   Junk_GC_Counter=0;
					   junk_character="";
			  elif(new_character=='G' or new_character=='C') and (prev_character=='A' or prev_character == 'T'):
				  # #print "new_character====" + new_character + "====prev_character====" + prev_character
				  prev_character = new_character
				  if i<int(sequence_length)-1 and int(len(temp))>=4:
					  if junk_character!="":
						 # #print "Helooo=====junk_character====" + junk_character
						 subset_array[counter,0]=junk_character # sub sequence.
						 subset_array[counter,1]=str(int(len(junk_character))) # length of sub sequence.
						 subset_array[counter,2]=str(i-int(len(temp))-int(len(junk_character))) # to set first postion of this subsequence formula is current position - length of sub sequence.
						 subset_array[counter,3]=str(i-int(len(temp))-1)
						 subset_array[counter,4]=str(Junk_AT_Counter)
						 subset_array[counter,5]=str(Junk_GC_Counter)
						 subset_array[counter,6]=str(3) # set sequnce type is GC
						
						 Junk_AT_Counter=0
						 Junk_GC_Counter=0
						 junk_character=""
						 counter = counter +1  # increment counter by one to set new subsequence
					 
					  subset_array[counter,0]=temp # sub sequence.
					  subset_array[counter,1]=int(len(temp)) # length of sub sequence.
					  subset_array[counter,2]=i-int(len(temp)) # to set first postion of this subsequence formula is current position - length of sub sequence.
					  subset_array[counter,3]=i-1 # to set last postion of this subsequence formula is current position - 1.
					  subset_array[counter,4]=int(len(temp))
					  subset_array[counter,5]=0
					  subset_array[counter,6]=1
					  Rich_Region_Index_Array.append(counter)
					  
					  #if subset_array[counter,1]>int(max_AT_length):
						  #max_AT_length=subset_array[counter,1]
						  #max_AT_array_counter=counter
						  
					  counter = counter +1  # increment counter by one to set new subsequence
					  
				  elif i<int(sequence_length)-1 and  int(len(temp)<4):
					  junk_character=junk_character + temp
					  Junk_AT_Counter=Junk_AT_Counter+int(len(temp))
				  
				  temp=new_character
				  prev_character=new_character
				  
				  if i==int(sequence_length)-1:
					  Junk_GC_Counter=Junk_GC_Counter+int(len(temp))
					  subset_array[counter,0]=junk_character# sub sequence.
					  subset_array[counter,1]=int(len(junk_character)) # length of sub sequence.
					  subset_array[counter,2]=i # to set first postion of this subsequence formula is current position - length of sub sequence.
					  subset_array[counter,3]=i
					  subset_array[counter,4]=Junk_AT_Counter
					  subset_array[counter,5]=Junk_GC_Counter
					  subset_array[counter,6]=3# set sequnce type is GC
					
					  counter = counter +1  # increment counter by one to set new subsequence
					
					  Junk_AT_Counter=0
					  Junk_GC_Counter=0
					  junk_character=""
			  elif((new_character=='G' or new_character=='C') and prev_character!='A' and prev_character != 'T'):
				  temp=temp + new_character
				  # #print "new_character====" + new_character + "====prev_character====" + prev_character
				  prev_character=new_character
				  
				  if i==int(sequence_length)-1 and (int(len(temp)>=4)):
					  if junk_character!="":
						  
						  subset_array[counter,0]=junk_character # sub sequence.
						  subset_array[counter,1]=int(len(junk_character)) # length of sub sequence.
						  subset_array[counter,2]=i-int(len(temp))-int(len(junk_character))+1; # to set first postion of this subsequence formula is current position - length of sub sequence.
						  subset_array[counter,3]=i-int(len(temp))
						  subset_array[counter,4]=Junk_AT_Counter
						  subset_array[counter,5]=Junk_GC_Counter
						  subset_array[counter,6]=3 # set sequnce type is GC
						
						  Junk_AT_Counter=0;
						  Junk_GC_Counter=0;
						  junk_character="";
						  counter = counter +1;  # increment counter by one to set new subsequence
					  
					  subset_array[counter,0]=temp;# sub sequence.
					  subset_array[counter,1]=int(len(temp)) # length of sub sequence.
					  subset_array[counter,2]=i-int(len(temp))+1 # to set first postion of this subsequence formula is current position - length of sub sequence.
					  subset_array[counter,3]=i
					  subset_array[counter,4]=0
					  subset_array[counter,5]=int(len(temp))
					  subset_array[counter,6]=2;
					  
					  #if subset_array[counter,1]>int(max_GC_length):
						 #max_GC_length=subset_array[counter,1]
						  #max_GC_array_counter=counter				  
					  counter = counter +1
				  elif i==sequence_length-1 and int(len(temp)<4):
						junk_character=junk_character + temp
						
						subset_array[counter,0]=junk_character# sub sequence.
						subset_array[counter,1]=int(len(junk_character)) # length of sub sequence.
						subset_array[counter,2]=i-int(len(junk_character)) # to set first postion of this subsequence formula is current position - length of sub sequence.
						subset_array[counter,3]=i
						subset_array[counter,4]=Junk_AT_Counter
						subset_array[counter,5]=Junk_GC_Counter
						subset_array[counter,6]=3
						counter = counter +1
						
						Junk_AT_Counter=0
						Junk_GC_Counter=0
						junk_character=""					  
			  elif((new_character=='A' or new_character=='T') and (prev_character=='C' or prev_character == 'G')):
				  # #print "new_character====" + newresult_sequence_character + "====prev_character====" + prev_character
				  prev_character = new_character
				  
				  if i<int(sequence_length)-1 and int(len(temp)>=4):
					  if junk_character!="":
						  # #print "junk_character====" + junk_character + "counter==" + str(counter)
						  
						  subset_array[counter,0]=junk_character # sub sequence.
						  subset_array[counter,1]=int(len(junk_character)) # length of sub sequence.
						  subset_array[counter,2]=i-int(len(temp))-int(len(junk_character)) # to set first postion of this subsequence formula is current position - length of sub sequence.
						  subset_array[counter,3]=i-int(len(temp))-1
						  subset_array[counter,4]=Junk_AT_Counter
						  subset_array[counter,5]=Junk_GC_Counter
						  subset_array[counter,6]= 3 # set sequnce type is GC
						  # #print "subset_array[counter,0]====" + subset_array[counter,0] + "counter==" + str(counter)
						  Junk_AT_Counter=0
						  Junk_GC_Counter=0
						  junk_character=""
						  counter = counter +1  # increment counter by one to set new subsequence						  	
					  # #print "temp====" + temp + "counter==" + str(counter)
					 
					  subset_array[counter,0]=temp # sub sequence.
					  subset_array[counter,1]=int(len(temp)) # length of sub sequence.
					  subset_array[counter,2]=i-int(len(temp)) # to set first postion of this subsequence formula is current position - length of sub sequence.
					  subset_array[counter,3]=i-1 # to set last postion of this subsequence formula is current position - 1.
					  subset_array[counter,4]=0
					  subset_array[counter,5]=int(len(temp))
					  subset_array[counter,6]=2;# Set the type tag value is AT
					  #if subset_array[counter,1]>max_GC_length :
						  #max_GC_length=subset_array[counter,1]
						  #max_GC_array_counter=counter						  
					  counter = counter +1
				  elif(i<(sequence_length)-1 and int(len(temp)<4)):
					    junk_character=junk_character + temp
					    Junk_GC_Counter=Junk_GC_Counter+ int(len(temp))					  
				  temp=new_character
				  prev_character=new_character
				  if i==int(sequence_length)-1 :
					   Junk_AT_Counter=Junk_AT_Counter+int(len(temp))
					   subset_array[counter,0]=junk_character# sub sequence.
					   subset_array[counter,1]=int(len(junk_character)) # length of sub sequence.
					   subset_array[counter,2]=i-int(len(junk_character)) # to set first postion of this subsequence formula is current position - length of sub sequence.
					   subset_array[counter,3]=i
					   subset_array[counter,4]=Junk_AT_Counter
					   subset_array[counter,5]=Junk_GC_Counter
					   subset_array[counter,6]=3 # set sequnce type is GC
					
					   counter = counter +1  # increment counter by one to set new subsequence
					
					   Junk_AT_Counter = 0
					   Junk_GC_Counter = 0
					   junk_character = ""					  
		  i=i+1
	  
	  array_size=counter
	  #for i in range (0,counter) :
		  #print "array_posaition====" + str(i)
		  #for j in range (0,7):
			  #print str(subset_array[i,j]) + "====="
		  #print "\n"
	  #print "array_size======"+str(array_size) +"\n"
	  counter=0	  
	  #writer = open("/home/afsana/Documents/help.txt", 'wb')
	  #print Rich_Region_Index_Array
	 
	  for i in range (0,len(Rich_Region_Index_Array)):
		  flag=1
		  max_AT_array_counter=Rich_Region_Index_Array[i]
		  """
		  #print "max_AT_array_counter===" + str(max_AT_array_counter)+"\n"
		  output_line_help= "max_AT_array_counter===" +str(max_AT_array_counter) + "***************\n"
		  writer.write(output_line_help)
		  output_line_help= "First Position===" +str(subset_array[max_AT_array_counter,2]) + "\n"
		  writer.write(output_line_help)
		  output_line_help= "Last Position===" +str(subset_array[max_AT_array_counter,3]) + "\n"
		  writer.write(output_line_help)
		  """
		  temp_array=self.GetRichRegion(subset_array,max_AT_array_counter,sequence,int(percentage))
		  result_sequence=temp_array[0,0]
		  result_sequence_length=temp_array[0,1]
		  result_sequence_First_Position = temp_array[0,2]
		  result_sequence_Last_Position = temp_array[0,3]
		  result_sequence_AT_percentage = temp_array[0,4]
		  check=self.ExactCheck(result_array,counter,result_sequence_First_Position,result_sequence_Last_Position,max_AT_array_counter)
		  if check==1:
			result_array[counter,0]=str(result_sequence)
			result_array[counter,1]=str(result_sequence_length)
			result_array[counter,2]=str(result_sequence_AT_percentage)
			result_array[counter,3]=str(result_sequence_First_Position)
			result_array[counter,4]=str(result_sequence_Last_Position)
			result_array[counter,5]=str(max_AT_array_counter)
			result_array[counter,6]="if max_start_position==1 :"
			counter=counter+1

		  
	  """	
	  for k in range(0,counter):
		  output_line_help="befor sorting===Rich_Region_Index_Array[i]====" + str(result_array[k,5]) +"***************\n"
		  writer.write(output_line_help)
		  output_line_help="k======" + str(k) +"***************\n"
		  writer.write(output_line_help)		  		  
		  #output_line_help="result_sequence====" + str(temp_array[k,0]) +"\n"
		  #writer.write(output_line_help)
		  #output_line_help="result_sequence_length====" + str(temp_array[k,1]) +"\n"
		  #writer.write(output_line_help)
		  #output_line_help="result_sequence_AT_percentage====" + str(temp_array[k,2]) +"\n"
		  #writer.write(output_line_help)
		  output_line_help=str(result_array[k,3]) +"\t\t\t\t"+str(result_array[k,4])+"\t\t\t\t"+str(result_array[k,6])+"\n"
		  writer.write(output_line_help)
		  #output_line_help="result_sequence_Last_Position====" + str(result_array[k,4]) +"***************\n"
		  #writer.write(output_line_help)
	  """
	  
	  #sort the result_array accoording first position
	  for i in range(0,counter):
		  #output_line_help="first loop iteration====="+str(i)+"\n************\n"
		  #writer.write(output_line_help)
		  for j in range(0,counter-1-i):
			  #output_line_help="second loop iteration====="+str(j)+"\n************\n"
			  #writer.write(output_line_help)
			  if int(result_array[j+1,3]) < int(result_array[j,3]):
				  #output_line_help="before swaped =====result_array[j,3]==="+str(result_array[j,3])+"******result_array[j+1,3]==="+str(result_array[j+1,3])+"\n*************\n" 
				  #writer.write(output_line_help)
				  result_array[j+1,0], result_array[j,0] = result_array[j,0], result_array[j+1,0]
				  result_array[j+1,1], result_array[j,1] = result_array[j,1], result_array[j+1,1]
				  result_array[j+1,2], result_array[j,2] = result_array[j,2], result_array[j+1,2]
				  result_array[j+1,3], result_array[j,3] = result_array[j,3], result_array[j+1,3]
				  result_array[j+1,4], result_array[j,4] = result_array[j,4], result_array[j+1,4]
				  result_array[j+1,5], result_array[j,5] = result_array[j,5], result_array[j+1,5]
				  result_array[j+1,6], result_array[j,6] = result_array[j,6], result_array[j+1,6]
				  #output_line_help="after swaped =====result_array[j,3]==="+str(result_array[j,3])+"******result_array[j+1,3]==="+str(result_array[j+1,3])+"\n*************\n\n"
				  #writer.write(output_line_help)
			
	  #print result_array	
	  #output_line_help=str(result_array[counter-1,3]) +"\t\t\t\t"+str(result_array[counter-1,4])+"\t\t\t\t"+str(result_array[counter-1,6])+"\t\t\t\t"+str(counter)+"\n"
	  #writer.write(output_line_help)
	  """
	  for k in range(0,counter):
		  output_line_help="after sorting===Rich_Region_Index_Array[i]====" + str(result_array[k,5]) +"***************\n"
		  writer.write(output_line_help)
		  output_line_help="k======" + str(k) +"***************\n"
		  writer.write(output_line_help)		  		  
		  #output_line_help="result_sequence====" + str(temp_array[k,0]) +"\n"
		  #writer.write(output_line_help)
		  #output_line_help="result_sequence_length====" + str(temp_array[k,1]) +"\n"
		  #writer.write(output_line_help)
		  #output_line_help="result_sequence_AT_percentage====" + str(temp_array[k,2]) +"\n"
		  #writer.write(output_line_help)
		  output_line_help=str(result_array[k,3]) +"\t\t\t\t"+str(result_array[k,4])+"\t\t\t\t"+str(result_array[k,6])+"\n"
		  writer.write(output_line_help)
		  #output_line_help="result_sequence_Last_Position====" + str(result_array[k,4]) +"***************\n"
		  #writer.write(output_line_help)		  
	  """	  
	  #coppy the result_array to a temporary array.
	  temp_array=result_array
	  #blank the result_array to get the unique data from temp_array.
	  result_array= {}
	  count=0
	  for k in range(0,counter):
		   if k==0:
			  result_array[k,0] =temp_array[k,0]
			  result_array[k,1] =temp_array[k,1]
			  result_array[k,2] =temp_array[k,2]
			  result_array[k,3] =temp_array[k,3]
			  result_array[k,4] =temp_array[k,4]
			  result_array[k,5] =temp_array[k,5]
			  result_array[k,6] =temp_array[k,6]
		   elif k>0:
			   insert_index=-1
			   for j in range (0,len(result_array)/7):
				   if int(temp_array[k,3])==int(result_array[j,3]) and int(temp_array[k,4])>int(result_array[j,4]):
					  insert_index=j
					  #output_line_help="insert_index" + str(insert_index) +"\n"
					  #writer.write(output_line_help)
				   elif int(temp_array[k,3])>int(result_array[j,4]):
					  insert_index=j+1
					  #output_line_help="insert_index" + str(insert_index) +"\n"
					  #writer.write(output_line_help)
				   elif int(temp_array[k,3])>int(result_array[j,3]) and int(temp_array[k,3])<int(result_array[j,4])and int(temp_array[k,4])<=int(result_array[j,4]):
					   insert_index=-2
					   #output_line_help="insert_index" + str(insert_index) +"\n"
					   #writer.write(output_line_help)
			   if insert_index>=0 :
				   result_array[insert_index,0] =temp_array[k,0]
				   result_array[insert_index,1] =temp_array[k,1]
				   result_array[insert_index,2] =temp_array[k,2]
				   result_array[insert_index,3] =temp_array[k,3]
				   result_array[insert_index,4] =temp_array[k,4]
				   result_array[insert_index,5] =temp_array[k,5]
				   result_array[insert_index,6] =temp_array[k,6]	  
	  """
	  for k in range(0,len(result_array)/7):
		  output_line_help="2nd time filtering=====" + str(result_array[k,5]) +"***************\n"
		  writer.write(output_line_help)
		  output_line_help="k======" + str(k) +"***************\n"
		  writer.write(output_line_help)		  		  
		  #output_line_help="result_sequence====" + str(temp_array[k,0]) +"\n"
		  #writer.write(output_line_help)
		  #output_line_help="result_sequence_length====" + str(temp_array[k,1]) +"\n"
		  #writer.write(output_line_help)
		  #output_line_help="result_sequence_AT_percentage====" + str(temp_array[k,2]) +"\n"
		  #writer.write(output_line_help)
		  output_line_help=str(result_array[k,3]) +"\t\t\t\t"+str(result_array[k,4])+"\t\t\t\t"+str(result_array[k,6])+"\n"
		  writer.write(output_line_help)
		  #output_line_help="result_sequence_Last_Position====" + str(result_array[k,4]) +"***************\n"
		  #writer.write(output_line_help)
	  """
	  j=0
	  counter=len(result_array)/7
	  #output_line_help="counter====" + str(counter) +"sequence_length====" + str(sequence_length)+"\n"
	  #writer.write(output_line_help)
	  if counter==0:
		  final_array[j,0]=sequence
		  final_array[j,1]=sequence_length
		  final_array[j,2]=0.00
		  final_array[j,3]=0
		  final_array[j,4]=sequence_length-1
		  final_array[j,5]="None"
		  final_array[j,6]="black"
		  final_array[j,7]="if counter==0:"
		  j=j+1
	  elif 	counter>0:	  	  	  
		  for i in range (0,counter) :
			  if int(result_array[i,3])==0 and int(result_array[i,4]) + 1 ==sequence_length and i == counter-1:
				  final_array[j,0]=result_array[i,0]
				  final_array[j,1]=result_array[i,1]
				  final_array[j,2]=result_array[i,2]
				  final_array[j,3]=result_array[i,3]
				  final_array[j,4]=result_array[i,4]
				  final_array[j,5]=str(result_array[i,5])
				  final_array[j,6]="green"
				  final_array[j,7]="if result_array[i,3]==0 and result_array[i,4] + 1 ==sequence_length and i == counter:"
				  j=j+1
			  elif int(result_array[i,3])==0 and int(result_array[i,4]) + 1 < sequence_length:		  	 
				  final_array[j,0]=result_array[i,0]
				  final_array[j,1]=result_array[i,1]
				  final_array[j,2]=result_array[i,2]
				  final_array[j,3]=result_array[i,3]
				  final_array[j,4]=result_array[i,4]
				  final_array[j,5]=result_array[i,5]
				  final_array[j,6]="green"
				  final_array[j,7]="elif result_array[i,3]==0 and result_array[i,4] + 1 < sequence_length:--1st part"
				  j=j+1
				  if i==len(result_array)/7:
					  temp_sequence= sequence[int(result_array[i,4]) +1:]
					  final_array[j,0]=temp_sequence[0 : sequence_length - int(result_array[i,4])-1]
					  final_array[j,1]=len(temp_sequence[0 : sequence_length - int(result_array[i,4])-1])
					  final_array[j,2]=-1.0
					  final_array[j,3]=result_array[i,4]
					  final_array[j,4]=sequence_length-1
					  final_array[j,5]="none"
					  final_array[j,6]="black"
					  final_array[j,7]="elif result_array[i,3]==0 and result_array[i,4] + 1 < sequence_length:---2nd part"
					  j=j+1		  			  
			  elif int(result_array[i,3])>0 and int(result_array[i,4]) + 1 == sequence_length:
				  if i==0:
					  temp_sequence= sequence[0:]
					  final_array[j,0]=temp_sequence[0 : int(result_array[i,3])]
					  final_array[j,1]=len(temp_sequence[0 : int(result_array[i,3])])
					  final_array[j,2]=-1.0
					  final_array[j,3]=0
					  final_array[j,4]=int(result_array[i,3])-1
					  final_array[j,5]="none"
					  final_array[j,6]="black"
					  final_array[j,7]="elif result_array[i,3]>0 and result_array[i,4] + 1 == sequence_length:if i==0:---1st part"
					  j=j+1
					  final_array[j,0]=result_array[i,0]
					  final_array[j,1]=result_array[i,1]
					  final_array[j,2]=result_array[i,2]
					  final_array[j,3]=result_array[i,3]
					  final_array[j,4]=result_array[i,4]
					  final_array[j,5]=result_array[i,5]
					  final_array[j,6]="green"
					  final_array[j,7]="elif result_array[i,3]>0 and result_array[i,4] + 1 == sequence_length:---2nd part"
					  j=j+1
				  elif i>0:
					  #print "seqlen==="+str(result_array[i,3]-result_array[i-1,4]-1)
					  #print "cutting position===" +str(result_array[i-1,4]+1)
					  #print "whol sequence ===" + str(sequence)
					  #print "sequence=====" +sequence[135:182]
					  temp_sequence= sequence[int(result_array[i-1,4]) +1:]
					  final_array[j,0]=temp_sequence[0: int(result_array[i,3])-int(result_array[i-1,4])-1]
					  final_array[j,1]=len(temp_sequence[0: int(result_array[i,3])-int(result_array[i-1,4])-1])
					  final_array[j,2]=-1.0
					  final_array[j,3]=int(result_array[i-1,4])+1
					  final_array[j,4]=int(result_array[i,3])-1
					  final_array[j,5]="none"
					  final_array[j,6]="black"
					  final_array[j,7]="elif result_array[i,3]>0 and result_array[i,4] + 1 == sequence_length:i>0:---1st part"
					  j=j+1
					  final_array[j,0]=result_array[i,0]
					  final_array[j,1]=result_array[i,1]
					  final_array[j,2]=result_array[i,2]
					  final_array[j,3]=result_array[i,3]
					  final_array[j,4]=result_array[i,4]
					  final_array[j,5]=result_array[i,5]
					  final_array[j,6]="green"
					  final_array[j,7]="elif result_array[i,3]>0 and result_array[i,4] + 1 == sequence_length:---2nd part"
					  j=j+1				  			 
			  elif int(result_array[i,3])>0 and int(result_array[i,4]) + 1<sequence_length :		  
				  if i==0:
					  temp_sequence= sequence[0:]
					  final_array[j,0]=temp_sequence[0 : int(result_array[i,3])]
					  final_array[j,1]=len(temp_sequence[0 : int(result_array[i,3])])
					  final_array[j,2]=-1.0
					  final_array[j,3]=0
					  final_array[j,4]=int(result_array[i,3])-1
					  final_array[j,5]="none"
					  final_array[j,6]="black"
					  final_array[j,7]="elif result_array[i,3]>0 and result_array[i,4] + 1<sequence_length :---1st part"
					  j=j+1
					  final_array[j,0]=result_array[i,0]
					  final_array[j,1]=result_array[i,1]
					  final_array[j,2]=result_array[i,2]
					  final_array[j,3]=result_array[i,3]
					  final_array[j,4]=result_array[i,4]
					  final_array[j,5]=result_array[i,5]
					  final_array[j,6]="green"
					  final_array[j,7]="elif result_array[i,3]>0 and result_array[i,4] + 1<sequence_length :---2nd part"
					  j=j+1
					  if int(result_array[i,4])<sequence_length	:
						  temp_sequence= sequence[int(result_array[i,4])+1:]
						  final_array[j,0]=temp_sequence[0 : sequence_length-int(result_array[i,4])]
						  final_array[j,1]=len(temp_sequence[0 : sequence_length-int(result_array[i,4])])
						  final_array[j,2]=-1.0
						  final_array[j,3]=int(result_array[i,4])+1
						  final_array[j,4]=sequence_length-1
						  final_array[j,5]="none"
						  final_array[j,6]="black"
						  final_array[j,7]="elif result_array[i,3]>0 and result_array[i,4] + 1<sequence_length : if int(result_array[i,4])<sequence_length	:---1st part"
						  j=j+1						    
				  elif i>0  and i < counter:
					  temp_sequence= sequence[int(result_array[i-1,4]) +1:]
					  final_array[j,0]=temp_sequence[0: int(result_array[i,3])-int(result_array[i-1,4])-1]
					  final_array[j,1]=len(temp_sequence[0 : int(result_array[i,3])-int(result_array[i-1,4])-1])
					  final_array[j,2]=-1.0
					  final_array[j,3]=int(result_array[i-1,4])+1
					  final_array[j,4]=int(result_array[i,3])-1
					  final_array[j,5]="none"
					  final_array[j,6]="black"
					  final_array[j,7]="elif i>0  and i < counter:---1st part"
					  j=j+1
					  final_array[j,0]=result_array[i,0]
					  final_array[j,1]=result_array[i,1]
					  final_array[j,2]=result_array[i,2]
					  final_array[j,3]=result_array[i,3]
					  final_array[j,4]=result_array[i,4]
					  final_array[j,5]=result_array[i,5]
					  final_array[j,6]="green"
					  final_array[j,7]="elif i>0  and i < counter:---2nd part"
					  j=j+1				  			  
				  elif i == counter:
					  final_array[j,0]=result_array[i,0]
					  final_array[j,1]=result_array[i,1]
					  final_array[j,2]=result_array[i,2]
					  final_array[j,3]=result_array[i,3]
					  final_array[j,4]=result_array[i,4]
					  final_array[j,5]=result_array[i,5]
					  final_array[j,6]="green"
					  final_array[j,7]="elif i == counter:---1st part"
					  j=j+1
					  temp_sequence= sequence[int(result_array[i-1,4]) +1:]				  
					  final_array[j,0]=temp_sequence[0 : sequence_length - int(result_array[i,4])-1]
					  final_array[j,1]=len(temp_sequence[0 : sequence_length - int(result_array[i,4])-1])
					  final_array[j,2]=-1.0
					  final_array[j,3]=result_array[i,4]
					  final_array[j,4]=sequence_length-1
					  final_array[j,5]="none"
					  final_array[j,6]="black"
					  final_array[j,7]="elif i == counter:---2nd part"
					  j=j+1				  						 
	  
	  final_array_size=j
	  """
	  for i in range (0,j) :
		  output_line_help= "array_posaition====" + str(i)+ "\n"
		  writer.write(output_line_help)
		  output_line_help= "sequence ====" +str(final_array[i,0]) + "\n"
		  writer.write(output_line_help)
		  output_line_help= "first_position ====" +str(final_array[i,3])+ "\n"
		  writer.write(output_line_help)
		  output_line_help= "last_position ====" +str(final_array[i,4])+ "\n"
		  writer.write(output_line_help)
		  output_line_help= "color ====" +str(final_array[i,6])+ "\n"
		  writer.write(output_line_help)
		  output_line_help= "condition ====" +str(final_array[i,7]) + "\n"
		  writer.write(output_line_help)
	  """
	  ofile_path = 	 output_file
	   
	  writer = open(ofile_path, 'wb')
	  output_line = "<html><body>"
	  writer.write(output_line)
	  output_line = "<table><tr><td><h5>Sequence name :  " + sequence_name + "</h5></td></td><tr><td>"
	  writer.write(output_line)
	  output_line=""
	  letter_count=0
	  string=""
	  
	  for i in range(0,final_array_size):
		 temp_seq= final_array[i,0]
		 output_line=output_line +"<b><font size='2' face='verdana' color='"+final_array[i,6]+"'>"
		 for k in range(0,int(final_array[i,1])):
			 substring=temp_seq[k:]
			 first_ch=substring[0:1]
			 if letter_count<69: 
				output_line=output_line +first_ch
				letter_count=letter_count+1
			 elif letter_count==69:
				output_line=output_line+first_ch+"</font></b></td></tr><tr><td><b><font size='2' face='verdana' color='" + final_array[i,6]+"'>"				
				writer.write(output_line)
				output_line=""
				letter_count=0				
		 output_line=output_line + "</font>"
	  output_line = "</b></td></tr></body></html>"
	  writer.write(output_line)	  
	  writer.close()
	  	  
	  success="yes"
	  return success	  


		
