#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import time,os
def p(n,k,p0,d):
	return (p0+(d*k*(n-k)/n**2))

def Psum(W,n,k,m):
	result=0
	for i in xrange(n+1,m.N+1):
		for j in xrange(max(k,i-m.N2),min(k+i-n,m.N1)+1):
			result = result + (2*p(i,j,m.P0,m.d)*W[i][j])/((j+1)*(i-j+1)-2)
	return result		

def Qsum(W,n,k,m):
	result=0
	if n==1:
		return result
	else:
		for i in xrange(1,n):
			for j in xrange(max(i+k-n,0),min(i,k)+1):
				result=result+(m.Pm*W[i][j]*W[n-i][k-j])
	return result			 	
 
class MergeSplit:
	P0=None
	d=None
	Pm=None
	N=None
	N1=None
	N2=None
	
	def __init__(self,N,N1,P0,d,Pm):
		self.P0=P0 #Base Split rate
		self.N=N #Total Population 
		self.N1=N1 #Type-I population
		self.d=d #delta
		self.Pm=Pm #Merge rate
		self.N2=N-N1


def main():
	for delta in [10,14,16,20]: #Edit delta in this list
		print (delta)
		start=time.time()
		tmax=100
		m = MergeSplit(100,25,1,delta,5) #Enter the values of N,N1,P0,delta,Pm here
		pop = m.N
		Wi = [[0 for j in xrange(0,i+1)] for i in xrange(0,pop+1)]
		Wf = [[0 for j in xrange(0,i+1)] for i in xrange(0,pop+1)]
		f = open("input100_f","r")
		
		i=0
		j=0
		for line in f:
			line=line.rstrip('\n')
			line=line.split('\t') 
			for j in xrange(0,len(line)):
				#print(i,j,float(line[j]))
				try:
					Wi[i][j]=float(line[j])
				except (KeyboardInterrupt, SystemExit):
					raise
				except:
					pass	
			i+=1
		#print Wi		
		
		for t in xrange(0,tmax): # FPI
			print t
			Z0=0
			for i in xrange(2,m.N+1): # Mean Field
				for j in xrange(max(0,i-m.N2),min(i,m.N1)+1):
					Z0 = Z0 + (p(i,j,m.P0,m.d)*Wi[i][j])
			Z0=Z0/m.Pm
			
			for n in xrange(1,pop+1):
				for k in xrange(0,n+1):
					Wf[n][k]=(Qsum(Wi,n,k,m)+(Psum(Wi,n,k,m)/Z0))/((p(n,k,m.P0,m.d)/Z0)+2*m.Pm)
			total=0
			for item in Wf:
				total=total+sum(item)
			for n in xrange(1,pop+1):
				for k in xrange(0,n+1):
					Wi[n][k]=Wf[n][k]/total 	
		
		filename = "FPI_N=%s_N1=%s_P0=%s_d=%s_Pm=%s.csv" % (m.N,m.N1,m.P0,m.d,m.Pm)						 
		g=open(filename,"w")
		for item in Wi:
			for subitem in item:
				g.write(str(subitem)+',')
			g.write('\n')
		g.close()	
		print time.time()-start
		
		foldername="N=%d_N1=%d_P0=%f_Pm=%f_d=%f" % (m.N,m.N1,m.P0,m.Pm,m.d)
		os.system('mkdir %s' % foldername)
		for i in xrange(0,len(Wi)):
			plt.figure(i)
			plt.plot(Wi[i])
			plt.savefig(foldername+'/'+str(i)+'.png')
			plt.close()
			
				
if __name__== "__main__": main()	
