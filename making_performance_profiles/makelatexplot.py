#########################
#  Given: a table with two columns of values,
#         perhaps runtimes for two different algorithms
#         on a collection of instances.
#  Note: a first row of titles is necessary, which is omitted
#        from the output.
#
#  Output: a latex file using tikz with three figures:
#          (a) a scatter plot of the values
#          (b) a performance profile
#          (c) a performance profile with logs in base 10.
#
#  Assume: all values are nonnegative.
#
#  Output choices:
#           (a) all values greater than cutoff are replaced by cutoff
#           (b) all ratios greater than 100 are omitted from both
#               performance profiles.
#########################
s=input("what is the file containing the data?\n")
f=open(s,'r')
print (f)
outfile=input("where do you want the output file?\n")
g=open(outfile,'w')
print(g)
import numpy
import math
a,b=numpy.loadtxt(f,skiprows=1,unpack=True)
#print(a)
#print(b)
f.close()
#  shrink numbers larger than cutoff
#increase zeroes
cutoff=3600
bada=numpy.zeros(a.size)
for i in range(0,a.size):
    g.write("%% {}\n".format(a[i]))
    if a[i] >= cutoff:
        bada[i]=1
for i in range(0,a.size):
    a[i]=min(a[i],cutoff)
    a[i]=max(a[i],0.05)
g.write("%%\n")
badb=numpy.zeros(b.size)
for i in range(0,b.size):
    g.write("%% {}\n".format(b[i]))
    if b[i] >= cutoff:
        badb[i]=1
for i in range(0,b.size):
    b[i]=min(b[i],cutoff)
    b[i]=max(b[i],0.05)
print(a.ndim)
print(a.itemsize)
print(a.size)
#print(bada)
#print(badb)
#print(a)
#print(b)
amax=max(a)
bmax=max(b)
allmax=max(amax,bmax)
allmaxf=float("{0:.1f}".format(allmax))
print(amax)
print(bmax)
print(allmax)
scale=12/allmax
print(scale)
# latex preamble stuff
g.write("\\documentclass{article}\n\n")
g.write("\\usepackage{graphicx,pgf,tikz}\n")
g.write("\\usetikzlibrary{fadings}\n\n")
g.write("\\begin{document}\n\n")
# scatter plot
g.write("\\begin{center}\n")
g.write("\\begin{tikzpicture}\n\n")
g.write("\\draw[thick] (0,0) -- (0,12) -- (12,12) -- (12,0) -- (0,0) -- (12,12);\n\n")
g.write("\\draw (0,0) node [below left] { 0 };\n")
g.write("\\draw (12,0) node [below] { ")
g.write("{}".format(allmaxf))
g.write(" };\n\n")
g.write("\\draw (0,12) node [left] { ")
g.write("{}".format(allmaxf))
g.write(" };\n\n")
g.write("\\draw (-0.1,11) node [left] {first};\n")
g.write("\\draw (-0.1,10.5) node [left] {column};\n")
g.write("\\draw (12,-0.4) node [below left] {second column};\n")
for i in range(0,a.size):
    aitem=float("{0:.2f}".format(scale*a[i]))
    bitem=float("{0:.2f}".format(scale*b[i]))
    g.write("\\fill ({},{}) circle(0.1);\n".format(bitem,aitem))
#  finish off the scatter plot
g.write("\n")
g.write("\\end{tikzpicture}\n")
g.write("\\end{center}\n\n")
# performance profile
yscale=float("{0:.3f}".format(10/(1+a.size)))
print(yscale)
c=b/a
c2=a/b
#print(c)
for i in range(0,a.size):
    c[i]=max(c[i],1)
    c2[i]=max(c2[i],1)
    if bada[i] >= 0.5:
        c2[i]=-1
    if badb[i] >= 0.5:
        c[i]=-1
#  don't plot any ratios larger than 100
cmax=max(c)
c2max=max(c2)
cc2max=max(cmax,c2max)
print(cc2max)
if cc2max > 100:
    for i in range(0,c.size):
        if c[i] >= 100:
            c[i]=-1
        if c2[i] >= 100:
            c2[i]=-1
    cc2max=100
print(cc2max)
xscale=float("{0:.3f}".format(11/(cc2max-1)))
print(xscale)
d=numpy.sort(c)
#print(d)
for i in range(0,d.size):
    d[i]=float("{0:.2f}".format(d[i]))
#print(c)
#print("d:\n")
#print(d)
d2=numpy.sort(c2)
for i in range(0,d2.size):
    d2[i]=float("{0:.2f}".format(d2[i]))
#print(c2)
#print("d2:\n")
#print(d2)
g.write("\\begin{center}\n")
g.write("\\begin{tikzpicture}\n\n")
g.write("\\draw[thick,->] (0,0) -- (0,10.5);\n\n")
g.write("\\draw[thick,->] (0,0) -- (12,0);\n\n")
cc2f=float("{0:.1f}".format(cc2max))
g.write("\\draw (11,0) node [below] { ")
g.write("{}".format(cc2f))
g.write(" };\n\n")
g.write("\\draw (11,-0.05) -- (11,0.05);\n")
g.write("\\draw (0,0) node [below] { 1 };\n")
g.write("\\draw (0,0) node [left] { 0 };\n")
g.write("\\draw (0,{}) node [left] ".format(a.size*yscale))
g.write(" { ")
g.write("{}".format(a.size))
g.write(" };\n\n")
g.write("\\draw (-0.1,8.5) node [left] {\\# instances};\n")
g.write("\\draw (12,-0.4) node [below left] {ratio};\n")
g.write("\\draw[thick, blue] (0,0) -- ")
negd=0
for i in range(0,d.size):
    if d[i] > 1.001:
        g.write("({},{}) -- ({},{}) -- ".format(xscale*(d[i]-1),(i-negd)*yscale,xscale*(d[i]-1),(i+1-negd)*yscale))
    else:
        if d[i] > 0:
            g.write("({},{}) -- ".format(xscale*(d[i]-1),(i+1-negd)*yscale))
        else:
            negd = negd+1
g.write("(11.1,{});\n\n".format((a.size-negd)*yscale))
negd=0
g.write("\\draw[thick, red, dotted] (0,0) -- ")
for i in range(0,d.size):
    if d2[i] > 1.001:
        g.write("({},{}) -- ({},{}) -- ".format(xscale*(d2[i]-1),(i-negd)*yscale,xscale*(d2[i]-1),(i+1-negd)*yscale))
    else:
        if d2[i] > 0:
            g.write("({},{}) -- ".format(xscale*(d2[i]-1),(i+1-negd)*yscale))
        else:
            negd=negd+1
g.write("(11.1,{});\n\n".format((a.size-negd)*yscale))
#  finish off performance profile
g.write("\n")
g.write("\\end{tikzpicture}\n")
g.write("\\end{center}\n\n")
#  performance profile with logs
for i in range(0,d.size):
    if d[i] > 0:
        d[i] = math.log10(d[i])
    if d2[i] > 0:
        d2[i] = math.log10(d2[i])
ldmax=max(d)
ld2max=max(d2)
lddmax=max(ldmax,ld2max)
lddmax=min(lddmax,2)
lddmax=math.log10(cc2max)
print(cc2max)
print(lddmax)
xscale=float("{0:.3f}".format(11/lddmax))
#print("log d:\n")
#print(d)
#print("log d2:\n")
#print(d2)
g.write("\\begin{center}\n")
g.write("\\begin{tikzpicture}\n\n")
g.write("\\draw[thick,->] (0,0) -- (0,10.5);\n\n")
g.write("\\draw[thick,->] (0,0) -- (12,0);\n\n")
cc2f=float("{0:.1f}".format(lddmax))
g.write("\\draw (11,0) node [below] { ")
g.write("{}".format(cc2f))
g.write(" };\n\n")
g.write("\\draw (11,-0.05) -- (11,0.05);\n")
g.write("\\draw (0,0) node [below left] { 0 };\n")
g.write("\\draw (0,{}) node [left] ".format(a.size*yscale))
g.write(" { ")
g.write("{}".format(a.size))
g.write(" };\n\n")
g.write("\\draw (-0.1,8.5) node [left] {\\# instances};\n")
g.write("\\draw (12,-0.4) node [below left] {$log_{10}$(ratio)};\n")
g.write("\\draw[thick, blue] (0,0) -- ")
negd=0
for i in range(0,d.size):
    if d[i] > 0.001:
        g.write("({},{}) -- ({},{}) -- ".format(xscale*d[i],(i-negd)*yscale,xscale*d[i],(i+1-negd)*yscale))
    else:
        if d[i] > -0.1:
            g.write("({},{}) -- ".format(xscale*d[i],(i+1-negd)*yscale))
        else:
            negd = negd+1
print(negd)
g.write("(11.1,{});\n\n".format((a.size-negd)*yscale))
negd=0
g.write("\\draw[thick, red, dotted] (0,0) -- ")
for i in range(0,d.size):
    if d2[i] > 0.001:
        g.write("({},{}) -- ({},{}) -- ".format(xscale*d2[i],(i-negd)*yscale,xscale*d2[i],(i+1-negd)*yscale))
    else:
        if d2[i] > -0.1:
            g.write("({},{}) -- ".format(xscale*d2[i],(i+1-negd)*yscale))
        else:
            negd=negd+1
print(negd)
g.write("(11.1,{});\n\n".format((a.size-negd)*yscale))
#  finish off performance profile
g.write("\n")
g.write("\\end{tikzpicture}\n")
g.write("\\end{center}\n\n")
#  finish off latex document
g.write("\\end{document}\n")
g.close()
