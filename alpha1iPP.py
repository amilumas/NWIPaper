import numpy as np
import makesemicrystallineLmpsfEMC as msc
import math


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
    diffangle = 99.5 - 90
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
    symmetries = np.array([[1,1,1], [-1,-1,-1], [1,-1+0.5, 1+0.5], [-1, 1+0.5, -1+0.5]])
    symmetriesM = np.array([[1,1,1], [-1,-1,-1], [1, 1, 1], [-1, -1,-1]])
    symmetriesA = np.array([[1,spacing,0], [1,2 - spacing,1], [0, 1+spacing , 0], [0, 1-spacing , 1]])
    #symmetriesNewM = np.array([[1,1,1],[-1,-1,-1],[1,-1,1],[-1,1,-1]])
    #symmetriesNewA = np.array([[]])
    unitcell = np.zeros((len(symmetries)*len(rcoords),3))
    ucount = 0
    for s in range(len(symmetries)):
        print("symmetries", symmetries[s])
        for ci in range(len(rcoords)):
            #unitcell[ucount,0] = rcoords[ci,0]*np.sign(symmetries[s,0]) + np.sign(symmetries[s,0])*abs(abs(symmetries[s,0])-1) 
            unitcell[ucount, 0] = rcoords[ci,0]*symmetriesM[s,0] + symmetriesA[s,0]
            unitcell[ucount, 1] = rcoords[ci,1]*symmetriesM[s,1] + symmetriesA[s,1]
            unitcell[ucount, 2] = rcoords[ci,2]*symmetriesM[s,2] + symmetriesA[s,2]
            #unitcell[ucount,1] = rcoords[ci,1]*np.sign(symmetries[s,1]) + np.sign(symmetries[s,1])*abs(abs(symmetries[s,1])-1)
            #unitcell[ucount,2] = rcoords[ci,2]*np.sign(symmetries[s,2]) + np.sign(symmetries[s,2])*abs(abs(symmetries[s,2])-1)
            print("unitcell[ucount]", unitcell[ucount])
            ucount +=1
    print("unitcell", unitcell) 
    center= np.zeros(3)
    center[0] = min(unitcell[:,0]) + (max(unitcell[:,0]) - min(unitcell[:,0]))/2
    center[1] = min(unitcell[:,1]) + (max(unitcell[:,1]) - min(unitcell[:,1]))/2
    print('min y', min(unitcell[:,1])*b, 'max y', max(unitcell[:,1])*b, 'center y', center[1]*b)
    center[2] = min(unitcell[:,2]) + (max(unitcell[:,2]) - min(unitcell[:,2]))/2
    print("center", center)
    count = 1
    mol = 1
    for k in range(Rz):
        #basez =k*c*zoffset
        basez = k*c
        for j in range(Ry):
            basey = j*b*2
            for i in range(Rx):
                #basex = i*a + xoffset*k*c*a 
                basex  = i*a*2
                for ci in range(len(unitcell)):
                        mol = ci//9 + 1
                        #x = basex + (symmetries[s][0]*a + rcoords[ci][0]*a) + xoffset*(symmetries[s][2]*c + rcoords[ci][2]*c)*a
                        x = basex + unitcell[ci,0]*a
                        #print("x", x,"basex", basex, "symmetries[s][0]*a", symmetries[s][0]*a, "symmetries[s][0]", symmetries[s][0])
                        y = basey + unitcell[ci,1]*b
                        #print("y", y, "basey", basey, "symmetries[s][1]*b", symmetries[s][1]*b, "symmetries[s][1]", symmetries[s][1])
                        #z = basez + (symmetries[s][2]*c + rcoords[ci][2]*c)*zoffset 
                        z = basez + unitcell[ci,2]*c
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
                    

    atomsinfo = np.array(atomsinfo)
    writexyz(xyzfile, atomsinfo)
    #xlo = -(Rx+1)*a
    #ylo = -(Ry+1)*b
    #zlo = -(Rz+1)*c
    #xhi =  (Rx+1)*a
    #yhi = (Ry+1)*b*(1.5)
    #zhi = (Rz+1)*c*(1.5)
    xlo = min(atomsinfo[:,3]) - 0.1*a
    xhi = max(atomsinfo[:,3]) + 0.1*a
    ylo = min(atomsinfo[:,4]) - 0.1*b
    yhi = max(atomsinfo[:,4]) + 0.1*b
    zlo = min(atomsinfo[:,5]) - 0.1*c
    zhi = max(atomsinfo[:,5]) + 0.1*c
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
    #readlammpsbondsPPctypes("TrialInfa1iPPbonds.data", "TrialInfa1iPPCtype.data")
    atomsinfo, xlo, xhi, ylo, yhi, zlo, zhi = setupInfiniteSystem("TrialInfa1iPP.xyz", 10, 10, 1)
    
    msc.writelammpsdatajustatoms("TrialInfa1iPP.data",[xlo,xhi,ylo,yhi,zlo,zhi], [15], len(atomsinfo), atomsinfo)
    #atomsinfo, xlo, xhi, ylo, yhi, zlo, zhi = readxyz("custompp-iso.xyz")
    #msc.writelammpsdatajustatoms("custompp-iso.data", [xlo, xhi, ylo, yhi, zlo, zhi], [14, 1], len(atomsinfo), atomsinfo)
main()
