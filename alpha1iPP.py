import numpy as np
import makesemicrystallineLmpsfEMC as msc
import math

def linearconnections(la,lb,la1, lb1, la2, lb2,start,end,mid,spacing):
    npoints = 1000
    xpos1 = np.linspace(start[0], mid[0], npoints)
    xpos2 = np.linspace(mid[0], end[0], npoints)
    ypos1 = la*xpos1 + lb*np.ones(npoints)
    ypos2 = la*xpos2 + lb*np.ones(npoints)
    zpos1 = la1*ypos1 + lb1 
    zpos2 = la2*ypos2 + lb2
    dists1 = np.zeros(npoints)
    dists2 = np.zeros(npoints)
    conx = []
    cony = []
    conz = []
    dists1[0] = 0
    
    nspacing = spacing*0.8
    for i in range(npoints-1):
        dists1[i+1] = dists1[i] + ((xpos1[i] - xpos1[i+1])**2 + (ypos1[i] - ypos1[i+1])**2 + (zpos1[i] - zpos1[i+1])**2)**0.5
        if dists1[i] > nspacing:
            conx.append(xpos1[i+1])
            cony.append(ypos1[i+1])
            conz.append(zpos1[i+1])
            nspacing = nspacing + spacing
    dists2[0] = dists1[-1]
    for i in range(npoints-1):
        dists2[i+1] = dists2[i] + ((xpos2[i] - xpos2[i+1])**2 + (ypos2[i] - ypos2[i+1])**2 + (zpos2[i] - zpos2[i+1])**2)**0.5
        print("ypos2[i+1]", ypos2[i+1], "end[1]-2", end[1]-2)
        if dists2[i] > nspacing and ypos2[i+1] <= end[1] - spacing/5:
            conx.append(xpos2[i+1])
            cony.append(ypos2[i+1])
            conz.append(zpos2[i+1])
            nspacing = nspacing + spacing
            
    return conx, cony, conz
    


def writexyz(xyzfile, atomsinfo):
    Natoms = len(atomsinfo)
    fout = open(xyzfile, 'w')
    fout.write(str(Natoms) + "\nAtoms\n")
    for i in range(Natoms):
        fout.write(str(int(atomsinfo[i,0])) + " " + str(int(atomsinfo[i,2])) + " " + str(atomsinfo[i,3]) + " " + str(atomsinfo[i,4]) + " " + str(atomsinfo[i,5]) + "\n")
    fout.close()

def setupInfiniteSystem(xyzfile, Rx, Ry, Rz):
    atomsinfo = []
    a = 6.63 
    b = 20.78
    c = 6.50
    atype = [1, 1, 1, 1, 1, 1, 1, 1, 1]
    beta = 99.5
    diffangle = 0
    #diffangle = 99.5 - 90
    #find unit vector of new c axis
    zoffset = np.cos(math.radians(diffangle))
    xoffset = -np.sin(math.radians(diffangle))
    print("xoffset", xoffset, "zoffset", zoffset)
    rcoords = np.zeros((9,3))
    rcoords[:,0] = [float(i) for i in "-0.0727 -0.0765 -0.1021 -0.3087 -0.1146 -0.1044 0.2775 0.0872 0.1026".split()]
    rcoords[:,1] = [float(i) for i in "0.2291 0.1592 0.1602 0.0589 0.0928 0.0854 0.0797 0.1156 0.1221".split()]
    rcoords[:,2] = [float(i) for i in "0.2004 0.2788 0.5098 0.4941 0.6057 0.8428 0.9260 0.9730 1.2109".split()]
    #rcoords = np.zeros((8,3))
    widthy = max(rcoords[:,1]) - min(rcoords[:,1])
    widthx = max(rcoords[:,0]) - min(rcoords[:,0])
    print("widthx", widthx)
    print("widthy",widthy)
    spacing = widthy/2
    print("widthy*4",widthy*4)

    #rcoords[:,0] = [0, 0, 0.5, 0.5, 0, 0, 0.5, 0.5]
    #rcoords[:,1] = [0, 0.5, 0, 0.5, 0, 0.5, 0, 0.5]
    #rcoords[:,2] = [0, 0, 0, 0, 0.5, 0.5, 0.5, 0.5 ]
    #rcoords = np.array([[0,0,0]])
    xnegmin = min(-rcoords[:,0])
    xposmin = min(rcoords[:,0])
    diffxmin = xposmin - xnegmin
    print("xnegmin", xnegmin, "xposmin", xposmin, "diffxmin", diffxmin)
    symmetries = np.array([[1,1,1], [-1,-1,-1], [1,-1+0.5, 1+0.5], [-1, 1+0.5, -1+0.5]])
    symmetriesM = np.array([[1,1,1], [-1, -1, -1], [1, 1, 1], [-1, -1, -1]])
    symmetriesA = np.array([[1,spacing,0], [diffxmin, 1-spacing, 1], [0, 1+spacing, 0], [1+diffxmin, 2- spacing, 1]])
    #symmetriesNewM = np.array([[1,1,1],[-1,-1,-1],[1,-1,1],[-1,1,-1]])
    #symmetriesNewA = np.array([[]])
    unitcell = np.zeros((len(symmetries)*len(rcoords),3))
    ucount = 0
    for s in range(len(symmetries)):
        print("symmetries", symmetries[s])
        for ci in range(len(rcoords)):
            #unitcell[ucount,0] = rcoords[ci,0]*np.sign(symmetries[s,0]) + np.sign(symmetries[s,0])*abs(abs(symmetries[s,0])-1) 
            #unitcell[ucount, 0] = rcoords[ci,0]*symmetriesM[s,0] + symmetriesA[s,0]
            unitcell[ucount, 1] = rcoords[ci,1]*symmetriesM[s,1] + symmetriesA[s,1]
            unitcell[ucount, 2] = (rcoords[ci,2]*symmetriesM[s,2] + symmetriesA[s,2])
            unitcell[ucount, 0] = rcoords[ci,0]*symmetriesM[s,0] + symmetriesA[s,0] 
            #unitcell[ucount,1] = rcoords[ci,1]*np.sign(symmetries[s,1]) + np.sign(symmetries[s,1])*abs(abs(symmetries[s,1])-1)
            #unitcell[ucount,2] = rcoords[ci,2]*np.sign(symmetries[s,2]) + np.sign(symmetries[s,2])*abs(abs(symmetries[s,2])-1)
            print("unitcell[ucount]", unitcell[ucount])
            ucount +=1
    print("unitcell", unitcell) 
    #calculate distances between ends:
    #pairs in unit cells
    ends1s = [[unitcell[8][0]*a, unitcell[8][1]*b, unitcell[8][2]*c], [unitcell[17][0]*a, unitcell[17][1]*b, unitcell[17][2]*c], [unitcell[26][0]*a, unitcell[26][1]*b, unitcell[26][2]*c], [unitcell[35][0]*a, unitcell[35][1]*b, unitcell[35][2]*c]]
    ends2s = [[unitcell[10][0]*a, unitcell[10][1]*b+0.05*b, unitcell[10][2]*c+0.1*c], [unitcell[19][0]*a, unitcell[19][1]*b+0.05*b, unitcell[19][2]*c+0.2*c], [unitcell[28][0]*a, unitcell[28][1]*b+0.05*b, unitcell[28][2]*c], [unitcell[1,0]*a, (unitcell[1,1]+2)*b+0.05*b, unitcell[1,2]*c+0.2*c]]
    midpoints = []
    for i in range(len(ends1s)):
        #print("end1", ends1s[i], "end2", ends2s[i])
        xdist = (ends1s[i][0] - ends2s[i][0])*a
        ydist = (ends1s[i][1] - ends2s[i][1])*b
        zdist = (ends1s[i][2] - ends2s[i][2])*c
        print("end1", ends1s[i][0]*a, ends1s[i][1]*b, ends1s[i][2]*c)
        print("end2", ends2s[i][0]*a, ends2s[i][1]*b, ends2s[i][2]*c)
        print("xdist", xdist, "ydist", ydist, "zdist", zdist)
        dist = (xdist**2 + ydist**2 + zdist**2)**0.5
        print("connection", i+1, "dist", dist)
    for i in range(len(ends1s)):
        print("end1", ends1s[i][0], ends1s[i][1], ends1s[i][2])
        print("end2", ends2s[i][0], ends2s[i][1], ends2s[i][2])
        if i%2 == 0:
            print("midpoint", (ends1s[i][0]+ ends2s[i][0])/2, (ends1s[i][1] + ends2s[i][1])/2 , (ends1s[i][2] + ends2s[i][2])/2)
            midpoints.append([(ends1s[i][0] + ends2s[i][0])/2, (ends1s[i][1] + ends2s[i][1])/2 , (ends1s[i][2] + ends2s[i][2])/2 + c])
            print("midpoints", midpoints)
        else:
            print("midpoint", (ends1s[i][0] + ends2s[i][0])/2, (ends1s[i][1] + ends2s[i][1])/2 , (ends1s[i][2] + ends2s[i][2])/2)
            midpoints.append([(ends1s[i][0] + ends2s[i][0])/2, (ends1s[i][1] + ends2s[i][1])/2 , (ends1s[i][2] + ends2s[i][2])/2 - c])
            print("midpoints", midpoints)

    print("unitcell", unitcell)
    #connection 1 
    #parabola in z y
    # linear in x
    Aparabola = []
    Hparabola = []
    Kparabola = []
    Aline = []
    Bline = []

    YZline1A = []
    YZline1B = []

    YZline2A = []
    YZline2B = []
    print('ends1s', ends1s, 'ends1s', ends2s, 'midpoints', midpoints)
    for i in range(len(ends1s)):
        print(ends1s[i][1], midpoints[i][1], ends2s[i][1], ends1s[i][2], midpoints[i][2], ends2s[i][2])
        ph  = midpoints[i][1]
        pk  = midpoints[i][2]
        pa = (ends1s[i][2] - pk)/((ends1s[i][1] - ph)**2)
        Aparabola.append(pa)
        Hparabola.append(ph)
        Kparabola.append(pk)
        la = (ends2s[i][1] - ends1s[i][1])/(ends2s[i][0] - ends1s[i][0])
        print("slope", la)
        lb = ends2s[i][1] - la*ends2s[i][0]
        #[la,lb] = np.polyfit([ends1s[i][1], ends2s[i][1]],[ends1s[i][0], ends2s[i][0]],1)
        Aline.append(la)
        Bline.append(lb) 
        
        yzla1 = (ends1s[i][2] - midpoints[i][2])/(ends1s[i][1] - midpoints[i][1])
        yzla2 = (ends2s[i][2] - midpoints[i][2])/(ends2s[i][1] - midpoints[i][1])
        yzlb1 = ends1s[i][2] -yzla1*ends1s[i][1]
        yzlb2 = ends2s[i][2] - yzla2*ends2s[i][1]
        YZline1A.append(yzla1)
        YZline1B.append(yzlb1)
        YZline2A.append(yzla2)
        YZline2B.append(yzlb2)


    connections1 = []
    connections2 = []
    connections3 = []
    connections4 = []
    sidegroups1 = []
    sidegroups2 = []
    sidegroups3 = []
    sidegroups4 = []
    npoints = [8,6,8,6]
    spacing = [1.7,1.9,1.7,1.9]
    for i in range(4):
        xs1 = np.linspace(ends1s[i][0], midpoints[i][0], npoints[i])
        ys1 = Aline[i]*xs1 + Bline[i]
        zs1 = YZline1A[i]*ys1 + YZline1B[i]
        xs2 = np.linspace(midpoints[i][0], ends2s[i][0], npoints[i])
        ys2 = Aline[i]*xs2 + Bline[i]
        zs2 = YZline2A[i]*ys2 + YZline2B[i]
        conx, cony, conz = linearconnections(Aline[i],Bline[i],YZline1A[i], YZline1B[i], YZline2A[i], YZline2B[i], ends1s[i], ends2s[i],midpoints[i],spacing[i])
        #print("old xs", xs, "old ys", ys, "old zs", zs)
        #xs,ys,zs = paraboliclinear(Aline[i],Bline[i],Aparabola[i],Hparabola[i],Kparabola[i],ends1s[i],ends2s[i],spacing[i]) 
        
        #print("xs", xs, "ys", ys, "zs", zs)
        for j in range(len(conx)):
            vars()["connections"+str(i+1)].append([conx[j], cony[j], conz[j]])
        #add side groups



    
    print("unitcell", unitcell)
    



    count = 1
    mol = 1
    for k in range(Rz):
        basez =k*c*zoffset
        #basez = k*c
        for j in range(Ry):
            basey = j*b*2
            for i in range(Rx):
                basex = i*a*2 + xoffset*k*c*a 
                #basex  = i*a*2
                for ci in range(len(unitcell)):
                        mol = ci//9 + 1
                        #x = basex + (symmetries[s][0]*a + rcoords[ci][0]*a) + xoffset*(symmetries[s][2]*c + rcoords[ci][2]*c)*a
                        x = basex + unitcell[ci,0]*a + unitcell[ci,2]*xoffset*c 
                        #print("x", x,"basex", basex, "symmetries[s][0]*a", symmetries[s][0]*a, "symmetries[s][0]", symmetries[s][0])
                        y = basey + unitcell[ci,1]*b
                        #print("y", y, "basey", basey, "symmetries[s][1]*b", symmetries[s][1]*b, "symmetries[s][1]", symmetries[s][1])
                        #z = basez + (symmetries[s][2]*c + rcoords[ci][2]*c)*zoffset 
                        z = basez + unitcell[ci,2]*c*zoffset
                        atomsinfo.append([count, mol, 1, x, y, z])
                        count = count + 1
                #atomsinfo.append([count, mol, 2, basex+ center[0]*a, basey + center[1]*b, basez + center[2]*c])
                #count = count + 1
                #atomsinfo.append([count, mol, 3, basex, basey, basez])
                #count = count + 1
                #atomsinfo.append([count, mol, 3, basex, basey+2*b, basez])
                #count = count + 1
                #atomsinfo.append([count, mol, 3, basex+2*a, basey, basez])
                #count = count + 1
    #for m in midpoints:
        #atomsinfo.append([count, 5, 2, m[0], m[1], m[2]])
        #count = count + 1
    for j in range(4):
        for i in range(len(vars()["connections"+str(j+1)])):
            atomsinfo.append([count, 5+j, 2+j, vars()["connections" + str(j+1)][i][0], vars()["connections" + str(j+1)][i][1], vars()["connections" + str(j+1)][i][2]])
            count = count + 1
                    

    atomsinfo = np.array(atomsinfo)
    writexyz(xyzfile, atomsinfo)
    xlo = min(atomsinfo[:,3]) - 5
    xhi = max(atomsinfo[:,3]) + 5
    ylo = min(atomsinfo[:,4]) - 5
    yhi = max(atomsinfo[:,4]) + 5
    zlo = min(atomsinfo[:,5]) - 5
    zhi = max(atomsinfo[:,5]) + 5
    #xlo = -(Rx+1)*a
    #ylo = -(Ry+1)*b
    #zlo = -(Rz+1)*c
    #xhi =  (Rx+1)*a
    #yhi = (Ry+1)*b*(1.5)
    #zhi = (Rz+1)*c*(1.5)
    #xlo = -2*a + Rx*c*zoffset
    #xhi = Rx*a 
    #ylo = -2*b
    #yhi = Rx*b
    #zlo = -2*c*zoffset
    #zhi = Rz*c*zoffset
    return atomsinfo, xlo, xhi, ylo, yhi, zlo, zhi

def readlammpsbondsPPctypes(filename,newfile):
    Cmasstypes = [13.019, 14.027, 15.035]
    boxcoords, masstypes, atoms, bonds, angles, dihedrals, atomsinfo, bondsinfo, anglesinfo, dihedralsinfo = msc.readlammpsdata(filename)
    #change atomtypes in based on number of bonds
    #loop through bondsinfo
    numBonds = np.zeros(len(atomsinfo))
    for i in range(len(bondsinfo)):
        bond = bondsinfo[i,:]
        atom1 = bond[2]-1
        atom2 = bond[3]-1
        numBonds[int(atom1)] += 1
        numBonds[int(atom2)] += 1

    #update atomsinfo
    for i in range(len(atomsinfo)):
        atomsinfo[i,2] = numBonds[i]

    #write new lammps data with just atoms 
    msc.writelammpsdatajustatoms(newfile, boxcoords, Cmasstypes, len(atomsinfo), atomsinfo)
        

def readxyz(filename):
    count = 0 
    with open(filename, 'r') as f:
        for i,line in enumerate(f):
            sstring = line.split()
            if i==0:
                N = int(sstring[0])
                atomsinfo = np.ones((N,6))
            elif i > 1:
                
                atomsinfo[count,0] = int(sstring[0])
                atomsinfo[count,1] = 1
                atomsinfo[count,2] = int(sstring[1])
                atomsinfo[count,3] = float(sstring[2])
                atomsinfo[count,4] = float(sstring[3])
                atomsinfo[count,5] = float(sstring[4])
                count = count + 1

    xlo = min(atomsinfo[:,3]) - 5
    xhi = max(atomsinfo[:,3]) + 5
    ylo = min(atomsinfo[:,4]) - 5
    yhi = max(atomsinfo[:,4]) + 5
    zlo = min(atomsinfo[:,5]) - 5
    zhi = max(atomsinfo[:,5]) + 5

    return atomsinfo, xlo, xhi, ylo, yhi, zlo, zhi

def main():
    readlammpsbondsPPctypes("TrialInfa1iPPbonds.data", "TrialInfa1iPPC1type.data")
    atomsinfo, xlo, xhi, ylo, yhi, zlo, zhi = setupInfiniteSystem("TrialInfa1iPP.xyz", 1, 1, 1)
    
    msc.writelammpsdatajustatoms("TrialInfa1iPP.data",[xlo,xhi,ylo,yhi,zlo,zhi], [15], len(atomsinfo), atomsinfo)
    #atomsinfo, xlo, xhi, ylo, yhi, zlo, zhi = readxyz("custompp-iso.xyz")
    #msc.writelammpsdatajustatoms("custompp-iso.data", [xlo, xhi, ylo, yhi, zlo, zhi], [14, 1], len(atomsinfo), atomsinfo)
main()
