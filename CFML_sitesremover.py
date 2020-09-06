#!/usr/bin/env python3
# coding: utf-8

import sys,os,logging,time
import argparse
parser = argparse.ArgumentParser(prog='CFML_siteremover',
				 usage='./CFML_siteremover.py [-h] [-f RECOMB_SITES] [-i INPUT_FASTA] [-o OUTPUT_FASTA] [-l SEQ_LENGTH] \n python3 CFML_siteremover.py [-h] [-f RECOMB_SITES] [-i INPUT_FASTA] [-o OUTPUT_FASTA] [-l SEQ_LENGTH]',
				 description='''\ Generates a fasta file that excludes sites of recombination from an alignment and records the list of sites excluded''',
				 epilog="Please send your suggestions and comments to Aditya <lankapalli@shh.mpg.de>")

parser.add_argument("-f", "--Recomb",dest="recomb_sites", help="Recombination Sites")
parser.add_argument("-i", "--input",dest="input_fasta", help="Input multifasta")
parser.add_argument("-o", "--output",dest="output_fasta", help="Output multifasta",default='CFML_recombsitesrm.fasta')
parser.add_argument("-l", "--length", dest="seq_length", help="Number of bases to be printed per line in output fasta [100]", type=int,default=100)
args = parser.parse_args()
logger=logging.getLogger()
logging.basicConfig(filename='CFML_siteremover.log',level=logging.DEBUG,format='%(asctime)s %(levelname)-8s %(message)s',filemode='w')

g=open(args.recomb_sites).read().strip().split('\n')
m=[list(map(int,i.split('\t')[1:3])) for i in g[1:]]
p=sorted(m,key=lambda x:x[0])
h2=list(filter(None,open(args.input_fasta).read().strip().split('>')))

h3={i.split('\n')[0]:''.join(i.split('\n')[1:]) for i in h2}
h5={i.split('\n')[0]:'' for i in h2}

chrstart=0
chrend=set([len(h3[numb]) for numb in h3])).pop()



def fu1(x,y):
    '''Identifies if segment1 is overlapping segment2
       and provides start and end of merged longer segment'''
    if x[1]<y[0] or x[0]>y[1]:
        return(x,y)
    else:
        return([min(x[0],y[0]),max(x[1],y[1])])

def frun1(a,j):
	'''checks the result from fu1 and stores them as a list '''
	q1=[]
	x,y=fu1(a,p[j])
	if type(x)==list:
		q.append(x)
		q1.append(j-1)
		a=y
		j+=1
		return(a,j)
	elif type(x)==int:
		a=[x,y]
		j+=1
		return(a,j)


def frun2(chrstart,chrend,q):
	''' re-orders list to find sites that are not excluded from the alignment'''
	r=[]
	if len(q)==1:
		r.append([chrstart,q[0][0]])
		r.append([q[0][1],chrend])
		return r
	else:
		for n,i in enumerate(q):
			if n==0:
				r.append([chrstart,q[n][0]])
				r.append([i[1],q[n+1][0]])
			elif n==(len(q)-1):
				r.append([q[n][1],chrend])
			else:
				r.append([i[1],q[n+1][0]])
		return r


def frun3(indict,outdict):
	'''Creates dictionary from the multi-alignment fasta dictionary '''
	outdict={}
	for n,i in enumerate(r):
		if n==0:
			outdict={j: indict[j][i[0]:i[1]] for j in indict}
		elif n!=0:
			for j in outdict:
				outdict[j]+=indict[j][i[0]:i[1]]
	return outdict


if __name__=='__main__':
	q=[];q1=[];j=1
	a=p[0]
	while j<(len(p)):
		a,j=frun1(a,j)

	q.append(a)
	q1.append(j-1)
	r=frun2(chrstart,chrend,q)
	hout=frun3(h3,h5)

	excluded_list=open('CFML_siteremover.excludedsites.txt','w')
	excluded_list.write('start\tend\n')
	for i in q:
		excluded_list.write(str(i[0])+'\t'+str(i[1])+'\n')

	excluded_list.close()
	output=open(args.output_fasta,'w')
	'''output.write('\n'.join(wrap(hout[out],args.seq_length))+'\n')'''
	for out in hout:
		output.write('>'+out+'\n')
		sequence=hout[out]
		logging.info("{0} Sequence length of {1}bp after {2} recombination sites were excluded".format(out,len(sequence),(len(h3[out])-len(sequence))))
		while len(sequence)>0:
			output.write(sequence[:args.seq_length]+'\n')
			sequence=sequence[args.seq_length:]

	output.close()

	sys.exit("Program completed successfully !!")
