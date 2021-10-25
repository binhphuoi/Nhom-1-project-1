import math
import numpy as np

np.set_printoptions(6, suppress=True)

tn = 3 #total nodes
te = 3 #total elements
xco = [0, 4472.135955, 0] #x co ordinate of nodes
yco = [0,0, 4000] #y co ordinate of nodes

A = 2300 ##mm2
E = 200000  ##Pa

snofel = [1,1,3] #start node of elements
enofel = [2,3,2] #end node of elements
lenofel = [] #length of the element
elcon = [] #constant of the element
cosofel = [] #cos of element
sinofel = [] #sin of element

for i in range(te):  
    a = snofel[i]
    b = enofel[i]
    x1 = float(xco[a-1])
    y1 = float(yco[a-1])
    x2 = float(xco[b-1])
    y2 = float(yco[b-1])
    l = math.sqrt((x2-x1)**2+(y2-y1)**2)
    con = A*E/l
    cos = (x2-x1)/l
    sin = (y2-y1)/l
    
    lenofel.append(l)
    elcon.append(con)
    cosofel.append(cos)
    sinofel.append(sin)
    
##print(lenofel)
##print(elcon)
##print(cosofel)
##print(sinofel)

elstmat = [] #element stiffness matrix

for i in range(te):
    cc = float(cosofel[i])**2
    ss = float(sinofel[i])**2
    cs = float(cosofel[i])*float(sinofel[i])
    
    mat = elcon[i]*np.array([[cc, cs, -cc, -cs],
                             [cs, ss, -cs, -ss],
                             [-cc, -cs, cc, cs],
                             [-cs, -ss, cs, ss]])
    elstmat.append(mat)
##print(elstmat)

gstmatmap = []                          ## Global stiffness matrix mapping, gstmatmap will be the sqare matrix of tn*
for i in range(te):                     ## do this for each elements
    m = snofel[i]*2                     ## taking the start node of element(i) and multiply by 2
    n = enofel[i]*2                     ## taking the end node of element(i) and multiply by 2
    add = [m-1, m, n-1, n]              ## Address of columns and rows of gstmatmap for elemet(i)
                                            # if startnode is 1 and end node is 2 then add=[1,2,3,4]
                                            # if startnode is 1 and end node is 3 then add=[1,2,5,6]
    gmat = np.zeros((tn*2, tn*2))       ## global stiffness matrix loaded with zeros for element(i)
    elmat = elstmat[i]                  ## taking the element stiffness matrix of element(i)
    for j in range(4):                  
        for k in range(4):              
            a = add[j]-1                ## addressing row of GST matrix for element(i)
            b = add[k]-1                ## addressing column of GST matrix for element(i)
            gmat[a,b] = elmat[j,k]      ## updating the values in GST matrix with EST matrix of element(i)
    gstmatmap.append(gmat)              ## storing the resultant matrix in gstmatmap list
##    print(np.around(gmat, 6))

GSM = np.zeros((tn*2, tn*2))            ## creating an empyty GSM matrix
for mat in gstmatmap:
    GSM = GSM+mat                       ## adding all the matrix in the gstmatmap list

print('\nGlobal Stiffness Matrix of the Truss\n', np.around(GSM, 6))

#-----------------------Boundry condition and Loading---------------------#

displist = []
forcelist = []
for i in range(tn):
    a = str('u')+str(i+1)
    displist.append(a)
    b = str('v')+str(i+1)
    displist.append(b)
    c = str('fx')+str(i+1)
    forcelist.append(c)
    d = str('fy')+str(i+1)
    forcelist.append(d)

##print(displist)
##print(forcelist)

dispmat = np.ones((tn*2,1))
dispmat[1*2-2, 0] = 0
dispmat[1*2-1, 0] = 0  
dispmat[3*2-2, 0] = 0
##print(dispmat)

##  _________________Loading____________________
forcemat = np.zeros((tn*2,1))
fx = 0
fy = -10000
forcemat[2*2-2, 0] = fx
forcemat[2*2-1, 0] = fy
##print(forcemat)    

###_________________Matrix Reduction_________________###

rcdlist = []
for i in range(tn*2):
    if dispmat[i,0] == 0:
        rcdlist.append(i)

rrgsm = np.delete(GSM, rcdlist, 0)              #row reduction
crgsm = np.delete(rrgsm, rcdlist, 1)            #column reduction
rgsm = crgsm                                    #reduced global stiffness matrix
rforcemat = np.delete(forcemat, rcdlist, 0)     #reduced force mat
rdispmat = np.delete(dispmat, rcdlist, 0)       #reduced disp mat

###_______________Solving____________________###

dispresult = np.matmul(np.linalg.inv(rgsm), rforcemat)
rin = 0
for i in range(tn*2):
    if dispmat[i,0] == 1:
        dispmat[i,0] = dispresult[rin,0]
        rin = rin+1
##print(dispmat)

forceresult = np.matmul(GSM, dispmat)
##print(forceresult)
print('\n***Positive is Tensile\nNegetive is Compressive***\n')
print('\nDisplacement matrix of nodes (mm)\n')
print(dispmat)
##print('\n\nForce matrix of nodes\n', forceresult)

##____________________new co ordinates of nodes____________####

newxco = []
newyco = []
count = 0
for i in range(tn):
    k = xco[i]+dispmat[count,0]
    newxco.append(k)
    count = count+1
    l = yco[i]+dispmat[count,0]
    newyco.append(l)
    count = count+1

###____________________new length of memebers______________####
    
newlenofel = []
for i in range(te):
    a, b = snofel[i], enofel[i]
    x1 = float(newxco[a-1])
    y1 = float(newyco[a-1])
    x2 = float(newxco[b-1])
    y2 = float(newyco[b-1])
    l = math.sqrt((x2-x1)**2+(y2-y1)**2)
    newlenofel.append(l)

##print(newlenofel)
##print(lenofel)

###______________strain in elements_______________________###
    
np.set_printoptions(6, suppress=False)

elstrain = np.zeros((te,1))
for i in range(te):
    elstrain[i,0] = (newlenofel[i]-lenofel[i])/(lenofel[i])

print('\n\nStrain in the elements')
print(elstrain)
np.set_printoptions(6, suppress=True)

###__________________stress in elements______________________###

elstress = np.zeros((te,1))
for i in range(te):
    elstress[i,0] = E*elstrain[i,0]
    
print('\nStress in the elements\n',elstress)

##np.delete(mat, row, 0)
##np.delete(mat, col, 1)
##np.matmul(mat1, mat2)
##np.linalg.inv(mat)
##np.arange(start, end, step)
##array.reshape(row, col)
##np.zeros((row,col))
##np.ones((row,col))
