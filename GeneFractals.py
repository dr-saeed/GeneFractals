### GENEFRACTALS - Software for Fractal Analysis of Genetic Sequences 
### Copyright (C) GPL-3.0 (2020) Dr. Mohammad Saeed

## Citation: Mohammad Saeed (2020). Fractal Genomics of SOD1 Evolution. doi: https://doi.org/10.1101/2020.05.14.097022.
## Link: https://www.biorxiv.org/content/10.1101/2020.05.14.097022v1.full

## This program is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.

## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.

## The GNU General Public License can be accessed at https://opensource.org/licenses/GPL-3.0.


# Import libraries

from __future__ import print_function
from __future__ import division
import os
import shutil
import time

import math
import numpy as np
import turtle as tu



# Store the list and order of files containing sequences from various organisms to be compared with reference file

filelist = []

def seqlist():
    
    global filelist
    
    fl = open('filelist.txt', 'r')
    for i in fl:
        i = i.rstrip()
        filelist.append(i)
        
    return filelist



# k-mer Analysis: Parse gene sequence into 10mer

def parseq():

    global genefile, seqlen

    seqlen = 0
    
    f1=open(genefile+'.fasta')

    i=f1.read()

    seqlen = len(i)

    f2=open('seqs.txt', 'w')
    
    seq = 10 # seq=input('Please set the length of nucleotide string sequence: ')

    iseq=int(seq)

    n=1

    seql=len(i)

    while seql>0:

        f2.write(i[n:n+iseq]+'\n')
        
        n=n+1
        
        seql=seql-1

    # print('Step 1 - complete: Sequences annotated. Analyzing...')

    f1.close()

    f2.close()

    

    f3=open('seqs.txt')

    statfile = os.path.join(path, 'Stats', genefile+'_Stats.txt')

    f4=open(statfile, 'w')

       

    minf = 0 # minf=input('Choose minimum frequency of occurence for sequence: ')
    
    minseq=9 # minseq=input('Minimum length of sequence to be considered: ')
    

    mf=int(minf)

    msq=int(minseq)


    counts=dict()

    for line in f3:

        words = line.split()

        for word in words:

            if word not in counts:

                counts[word] = 1

            else:

                counts[word] += 1

    lst = list()
    keyf = [] # Add for freq of k-mer
    
    for key, val in counts.items():

        if val>mf and len(key)>msq:

            lst.append( (val, key) )
            keyf.append(val) # Add for freq of k-mer

            
    lst.sort(reverse=True)
    kmers = sum(keyf)

    freq = {x:keyf.count(x) for x in keyf} # Add for freq of k-mer
    print (genefile, kmers, freq, file=f5) # Add for freq of k-mer

    for key, val in lst :

        ll=len(val)

        print(val, ',', key, file=f4)

    # print('Sliding Window Analysis Complete: ', genefile)
    print(genefile, seqlen)

    
    f4.close()

    f3.close()

    return seqlen



## Function to compare gene sequences of a list of organisms with a ref gene seq

def Seqcomp():

    # Read the reference parsed-sequence file and store the sequences as a list (maintain order)

    global genefile

    try:
        os.chdir('Stats')
        file1 = open (genefile+'_Stats.txt')
        os.chdir('..')
        fout = open ('SeqTab.txt', 'w')

        print ("", file=fstat)
        

    except:
        print ('File not found:', genefile)
        exit()

# Store ref sequence in a list for comparison:

    refseq = []
    seq2 = {}
    
    for line in file1:
        var = line.split(',')
        refseq.append(var[0])

    print ('')

    print (genefile, '\t\t', '# Matched Seq', '\t\t', '% Match')
    print (genefile, '\t', '# Matched Seq', '\t', '% Match', file=fstat)
    
    for s in refseq:
        print (s, end=' ', file=fout)

# Open the files in the filelist one by one and store their sequences and their frequencies in a dictionary
# Run the function chk_list to compare these with the reference sequence

    global filelist, H

    os.chdir('Stats')

    
    try:
        for fname in filelist:
            seq2 = {}
            file2 = open(fname+'_Stats.txt')
            print ("", file=fout)
            print (fname, end=' ', file=fout)
           
            for l in file2:
                l = l.rstrip()
                v = l.split(',')
                seq2[v[0]] = v[1]
                
            # Chk_list():

            pos = 0
            neg = 0
                        
            for key in refseq:
                if key in seq2:
                    print (seq2[key], end=' ', file=fout)
                    pos += 1
                else:
                    print (' 0', end=' ', file=fout)
                    neg += 1
                    
            perhomol = float(pos*100)/(pos+neg)
            H.append(round(perhomol, 2))

            
            print ('\t', fname, '\t\t', pos, '\t\t', round(perhomol, 2))
            
            print (fname, '\t', pos, '\t', round(perhomol, 2), file=fstat)
            
       
            
    except:
        print ('File not found:', fname)
        file2.close()
        exit()

    file2.close()

    
    os.chdir('..')
    fout.close()
    


## Detrended Fluctuation Analysis (DFA) for Genomic Sequences

def dfa():
    
    global genefile
    
    fi=open(genefile+'.fasta')

    fout=open('dfaseq.txt', 'w')
     
    while True:
        n=fi.read(1)
                
        if n=='G':
            print ('3', file=fout)
            
        if n=='C':
            print ('4', file=fout)
            
        if n=='A':
            print ('1', file=fout)
            
        if n=='T':
            print ('2', file=fout)
            
        if not n: break

    
    fi.close()
    fout.close()

    os.startfile('dfabat')




## Relative Dispersion Analysis (RDA) for Genomic Sequences

# Calculate RDA

def rda():

    global gene, geneseq, xy

    seql = len(geneseq)

    logseq = math.log(seql, 2)

    box = int(logseq)

    # print (gene, 'Seq length = ', seql, 'Max Box =', box) # Uncomment if RDA processing needs to be visualized
    print (gene, 'Seq length = ', seql, 'Max Box =', box, file=rdaout)
    # print ('')
    print ('', file=rdaout)

    # print ('Box', '\t', 'Mean', '\t', 'SD', '\t', 'RD%', '\t', 'ln(2^Box)', '\t', 'ln(RD%)') # Uncomment if RDA processing needs to be visualized
    print ('Box', '\t', 'Mean', '\t', 'SD', '\t', 'RD%', '\t', 'ln(2^Box)', '\t', 'ln(RD%)', file=rdaout)

    
    m = range(1, box)
    
    for n in m:

        i = 0
        r = 0
        ave = []
       
        while r <= 2**box:

            r += 2**n
                            
            mseq = np.mean(geneseq[i:r])                            
            ave.append(mseq)

            i += 2**n

        
        sc = np.std(ave)
        ac = np.mean(ave)
        ratio = (sc * 100)/ ac
        logbox = math.log(2**n)
        logratio = math.log(ratio)

        xy[n] = logbox, logratio

        # print (n, '\t', round(ac, 2), '\t', round(sc, 2), '\t', round(ratio, 2), '\t', round(logbox, 2), '\t', round(logratio, 2)) # Uncomment if RDA processing needs to be visualized
        print (n, '\t', round(ac, 2), '\t', round(sc, 2), '\t', round(ratio, 2), '\t', round(logbox, 2), '\t', round(logratio, 2), file=rdaout)

    return xy       



# Linear regression of RDA / DFA output

g = 0
c = 0
r = 0


def linreg():

    global gene, g, c, r, xy
    
       
    sx=0
    sy=0
    sxx=0
    syy=0
    sxy=0
    n=0

    for k, v in sorted(xy.items()):
        x = v[0]
        y = v[1]

            
        sx=sx+x
        sy=sy+y
        x2 = x**2
        sxx=sxx+x2
        y2=y**2
        syy=syy+y2
        xy1 = x*y
        sxy=sxy+xy1
        n +=1
        
    
    # s2xy=(sxy/n)-((sx*sy)/(n**2))
    a1=sxy/n
    a2=sx*sy
    a3=n**2
    a4 = a2/a3
    s2xy = a1-a4   

    # s2x=(sxx/n)-((sx**2)/(n**2))
    b1=sxx/n
    b2=sx**2
    b4 = b2/a3
    s2x = b1-b4

    # s2y=(syy/n)-((sy**2)/(n**2))
    c1=syy/n
    c2=sy**2
    c4 = c2/a3
    s2y = c1-c4
        
    # r=(s2xy) /((s2x*s2y)**.5)
    d1 = s2x*s2y
    d2 = (d1)**0.5
    r = s2xy/d2
    
    # g=(sx*sy-(n*sxy))/((sx**2)-(n*sxx))
    e1=sx*sy
    e2=n*sxy
    e3=sx**2
    e4=n*sxx
    e5=e1-e2
    e6=e3-e4
    g=e5/e6

    # c=((sx*sxy)-(sy*sxx))/((sx**2)-(n*sxx))
    h1=sx*sxy
    h2=sy*sxx
    h3=sx**2
    h4=n*sxx
    h5= h1-h2
    h6=h3-h4
    c= h5/h6

    
    return g, c, r





## Fractal Diagram for Genomic Sequences


# Create list of reference gene sequence with nucleotides coded as numbers
# Ref: {ORIGINAL GW-BASIC CODE A=1 (N); T=2 (S); G=3 (E); C=4 (W)}
# A, T, G, C = 1, 2, 3, 4


def seqnuc():

    global gene, geneseq
    
    fgene=open(gene+'.fasta', 'r')

    cnt = 0
        
    while True:
        nuc=fgene.read(1)
        cnt +=1

        if nuc=='A':
            n=1
            
        elif nuc=='T':
            n=2
            
        elif nuc=='G':
            n=3
            
        elif nuc=='C':
            n=4

            
        geneseq.append(n)
                
        if not nuc: break

    fgene.close()

    return geneseq



def fracgene():

    global fracseq, box
    
    stpi = 4 # open diagram = 1

    # Start canvas and get ready to draw sequence
        
    # Imp Note: tu.mode('logo') - offsets screen. Therefore do not use.

    tu.hideturtle()
    tu.speed('fastest')
    tu.tracer(False)


    ## Set box size

    wni = int((8*(2**box))**0.5)
    # wni = 200 # for closed diagram
    # wni=5000 # for open diagram

    tu.screensize(wni, wni) # for closed diagram

    ##tu.pu() # open diagram
    tu.setpos(0, 0)
    ##tu.pd() # open diagram

    tu.color('green')
    tu.begin_fill()
    tu.circle(5)
    tu.end_fill()

    # For open diagram comment out this region
    tu.color('grey')
    tu.pu()
    tu.setpos(-wni-stpi, wni+stpi)

    tu.seth(0)
    tu.pd()
    for border in range(1, 5):
        tu.fd(2*wni+(2*stpi))
        tu.rt(90)
    tu.pu()
    tu.home()
    tu.pd()
    # Till here

    tu.color('blue')


    # Draw Fractal Gene Sequence

    for nuc in fracseq:
        # For open diagram comment out this region    
        tpx = tu.xcor()
        tpy = tu.ycor()

        if tpx >=wni:
            tu.pu()
            tu.setx(-wni)
            tu.pd()
        elif tpx <=-wni:
            tu.pu()
            tu.setx(wni)
            tu.pd()
        elif tpy >=wni:
            tu.pu()
            tu.sety(-wni)
            tu.pd()
        elif tpy <=-wni:
            tu.pu()
            tu.sety(wni)
            tu.pd()
        # Till here
            
        # ORIGINAL GW-BASIC CODE A=1 (N); T=2 (S); G=3 (E); C=4 (W)
        # In Py mode: N -> 90; S -> 270; E -> 0; W ->180

        if nuc==3:    # 'G'
            tu.seth(0)
            tu.forward(stpi)
            
        elif nuc==4:    # 'C'
            tu.seth(180)
            tu.forward(stpi)
            
        elif nuc==1:    # 'A'
            tu.seth(90)
            tu.forward(stpi)
            
        elif nuc==2:    # 'T'
            tu.seth(270)
            tu.forward(stpi)
            
        

    # End drawing on canvas and save image

    tu.color('red')
    tu.begin_fill()
    tu.circle(5)  
    tu.end_fill()
    tu.hideturtle()
    tu.pu()

    ts = tu.getscreen()
    diagfile = os.path.join(path, 'Diag', gene+'_'+str(box)+"_image.eps")

    ts.getcanvas().postscript(file=diagfile)

    

    tu.clear()
    # print (gene)
    
    # print ("Use EPS converter to convert your .eps image to JPEG: https://www.epsconverter.com/")
    # print ("OR Download: https://epsviewer.org/downloadfinal.aspx")

    # tu.done() # uncomment if you want to visualize the image on screen but it will load memory and program can get stuck




## Main

# Create path

path = os.getcwd()
print (path)




# k-mer Analysis

if not os.path.exists('Stats'):
    os.mkdir('Stats')
else:
    print ('Stats folder already exists')
    
seqlist()

print ()
print ()
print ('Results: Parsed Sequences')
print ()


f5 = open ('kmer_stats.txt', 'w')
sqlen = []

for genefile in filelist:
    parseq()
    sqlen.append(seqlen)

f5.close()

# Uncomment this block if you want only 1 seq compared and not all seq with each other

print ()
print ()
print ('Results: Homolog')

fstat = open ('Homolog.txt', 'w')
ftab = open ('Tab.txt', 'w')

spp = 0
for genefile in filelist:
    H = []
    Seqcomp()
    print (genefile, '\t', H, '\t', sqlen[spp], file=ftab)
    spp +=1
     
fstat.close()
ftab.close()

os.remove ("seqs.txt")
os.remove ("SeqTab.txt")


# Seqcomp() # comment out this line for comparing all seqs




## Excute DFA, RDA, Fractal Diagram code


# RDA - Statistically evaluate Self Similarity and Long range order

rdaout = open ('RDA.txt', 'w')

print ('')

Dseq = {}

print ('')
print ('Results: RDA')
print ('')
print ('Gene', '\t ', 'Gene_length', '\t', 'Max Box', '\t', 'D', '\t', 'c', '\t', 'r^2', '\t', 'AT%', '\t', 'GC%')

    
for gene in filelist:
    geneseq = []
    seqnuc()
    xy = {}
    rda()
    linreg()
    
    print (gene, '\t ', round(1+(g*-1), 3), '\t', round(c, 2), '\t', round(r**2, 2), file=rdaout)
    print ('', file=rdaout)
    
    # freq = np.unique(geneseq, return_counts=True) # this will also give A,T,G,C frequencies (as arrays) 
    
    freq = {x:geneseq.count(x) for x in geneseq} # Obtaining freq of A,T,G,C using a dictionary
    # print ('')
    # print ('A:', freq[1], '\t', 'T:', freq[2], '\t', 'G:', freq[3], '\t', 'C:', freq[4])

    print ('', file=rdaout)
    print ('A:', freq[1], '\t', 'T:', freq[2], '\t', 'G:', freq[3], '\t', 'C:', freq[4], file=rdaout)

    AT = int(freq[1])+int(freq[2])
    GC = int(freq[3])+int(freq[4])
    AT_p = AT*100/len(geneseq)
    GC_p = GC*100/len(geneseq)

    # print ('A+T:', AT, '\t', 'G+C:', GC, '\t', 'AT%:', round(AT_p, 2), '\t', 'GC%', round(GC_p, 2))
    print (gene, '\t ', len(geneseq), '\t', int(math.log(len(geneseq), 2)), '\t', round(1+(g*-1), 3), '\t', round(c, 2), '\t', round(r**2, 2), round(AT_p, 2), '\t', round(GC_p, 2))

    print ('A+T:', AT, '\t', 'G+C:', GC, '\t', 'AT%:', round(AT_p, 2), '\t', 'GC%', round(GC_p, 2), file=rdaout)
    print ('', file=rdaout)
    print ('', file=rdaout)

    
##    This block of code can change name of keys from 1, 2, 3, 4 to A, T, G, C
##    ntd = ['A', 'T', 'G', 'C']
##    freq=dict(zip(ntd,list(freq.values())))
##    print (freq)
##    print ('')
    

    Dseq[gene] = len(geneseq), 1+(g*-1), AT_p, GC_p # Dictionary for later statistics in Excel

print ('Gene', '\t ', 'Gene_length', '\t', 'D', '\t', 'AT%', '\t', 'GC%', file=rdaout)    
for k, v in sorted(Dseq.items()):
    # print (k, '\t', v[0], '\t', round(v[1], 2), '\t', round(v[2], 2), '\t', round(v[3], 2))
    print (k, '\t', v[0], '\t', round(v[1], 2), '\t', round(v[2], 2), '\t', round(v[3], 2), file=rdaout)

rdaout.close()


# Calculate DFA
    
print ()
print ()
path = os.getcwd()

if not os.path.exists('DFA'):
    os.mkdir('DFA')
else:
    print ('DFA folder already exists')
    shutil.rmtree(os.path.join(path, 'DFA'))
    os.mkdir('DFA')
 
    
print ()
print ('Results: DFA')
print ()
print ('Species', '\t ', 'alpha', '\t', 'c', '\t', 'r2')
print ()

dfaout = open('DFA_out.txt', 'w')
print ('Species', '\t ', 'alpha', '\t', 'c', '\t', 'r2', file=dfaout)

for genefile in filelist:
    xy = {}
    dfa()
    time.sleep(3)
    
    fdfa = open('dfa.txt', 'r')

    line = 0

    for i in fdfa:
        i = i.rstrip()
        v = i.split(' ')
        var1 = float(v[0])
        var2 = float(v[1])
        xy[line] = var1, var2
        line +=1

    fdfa.close()
    shutil.move('dfa.txt', os.path.join(path, 'DFA'))
    os.rename (os.path.join(path, 'DFA', 'dfa.txt'), os.path.join(path, 'DFA', genefile+'_DFA.txt'))
    

    linreg()
    
    print (genefile, '\t ', round(g, 3), '\t', round(c, 2), '\t', round(r**2, 2))
    print (genefile, '\t ', round(g, 3), '\t', round(c, 2), '\t', round(r**2, 2), file=dfaout)
    # print (genefile, ': ', 'alpha = ', g, 'c = ', c, 'r2 = ', r**2, file=dfaout)

    
dfaout.close()
os.remove ("dfaseq.txt")


# Gene Fractal Diagram - Graphically evaluate Self Similarity

print ()
path = os.getcwd()
if not os.path.exists('Diag'):
    os.mkdir('Diag')
else:
    print ('Diag folder already exists')
    
for gene in filelist:
    geneseq = []
    b=0
    seqnuc()
        
    b = int(math.log(len(geneseq), 2))
    if b < 13:
        box = b
    else:
        box = 13
        
    fracseq = geneseq
    fracgene()

    box = 8
        
    while box <= b or box == 13:
        fracseq = geneseq[0:2**box]
        fracgene()
        box += 1


print ('')
print ('Gene Fractal Diagrams - completed')        

tu.done()
