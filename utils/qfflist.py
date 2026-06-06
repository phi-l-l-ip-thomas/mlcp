import sys
import json

tol = 1E-8

def readthespectrofile(fnm):
    with open(fnm, 'r') as file:
       f = json.load(file)
    f2 = f["freq"][0]
    f3 = f["f3qcm"]
    f4 = f["f4qcm"]
    return f2,f3,f4

def harmoniclist(nmode):
    "Generates index list for harmonic terms"
    fulllist=[]
    for ii in range(nmode):
       ans = [ii+1 for jj in range(2)]
       fulllist.append(ans)
    return fulllist

def qffindxlist(nmode,ncoup):
    "Generates QFF index list"
    fulllist=[]
    ans = [0 for ii in range(ncoup)]
    reset = False
    while not reset:
       ans_prt = [ans[jj]+1 for jj in range(len(ans))]
       fulllist.append(ans_prt)
       ans,reset = nextqffindx(ans,nmode)
    return fulllist

def nextqffindx(idx,nmode):
    "Find next QFF index"
    ncoup = len(idx)
    reset = False
    ans = [idx[ii] for ii in range(ncoup)]
    for ii in range(ncoup-1,-1,-1):
        if ii == 0:
            if ans[ii] == nmode-1:
               ans = [0 for ii in range(ncoup)]
               reset = True
            else: 
               ans[ii] += 1
               for jj in range(ncoup-1,ii,-1):
                   ans[jj] = 0
        elif ans[ii] < ans[ii-1]:
            ans[ii] += 1
            for jj in range(ncoup-1,ii,-1):
                ans[jj] = 0
            break
    return ans,reset

def writefcfiles(tag,inds,fc):
    "Writes force constants to file"
    if len(inds) != len(fc):
       print("ERROR: force constant and index lengths differ!")
       return False
    ncoup = len(inds[0])
    fnm="f"+str(ncoup)+tag+".dat"
    f = open(fnm,'w')
    for ii in range(len(fc)):
        val = float(fc[ii])
        if abs(val) > tol:
           for jj in range(ncoup):
               f.write(" %3d " %(inds[ii][jj]))
           f.write(" %16.8f \n" %(val))
    f.close()
    return True

def CreateQffFileMain(nm):
    "Main routine for creating QFF files"
    f2,f3,f4 = readthespectrofile('finish.spectro')
    nmode = len(f2)
    i2 = harmoniclist(nmode)
    i3 = qffindxlist(nmode,3)
    i4 = qffindxlist(nmode,4)
    ok = writefcfiles(nm,i2,f2)
    ok = writefcfiles(nm,i3,f3)
    ok = writefcfiles(nm,i4,f4)
    return

if len(sys.argv) == 2:
   CreateQffFileMain(sys.argv[1])
else:
   print("syntax: '$> %s <system name>'" %(sys.argv[0]))

