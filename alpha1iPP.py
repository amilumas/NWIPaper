import numpy as np
import makesemicrystallineLmpsfEMC as msc
import math


def linear2parabolaXZConnect(la, lb, pa1, ph1, pk1, pa2, ph2, pk2, start, end, mid, spacing):
    npoints = 1000
    xpos1 = np.linspace(start[0], mid[0], npoints)
    xpos2 = np.linspace(mid[0], end[0], npoints)
    ypos1 = la*xpos1 + lb*np.ones(npoints)
    ypos2 = la*xpos2 + lb*np.ones(npoints)
    zpos1 = pa1*(xpos1 - ph1)**2 + pk1
    zpos2 = pa2*(xpos2 - ph2)**2 + pk2
    #print("xpos1", xpos1)
    #print("zpos1", zpos1)
    dists1 = np.zeros(npoints)
    dists2 = np.zeros(npoints)
    conx = []
    cony = []
    conz = []
    sidex = []
    sidey = []
    sidez = []
    dists1[0] = 0
    pad = 1.5
    extra = 0
    nspacing = spacing*0.8
    for i in range(npoints-1):
        dists1[i+1] = dists1[i] + ((xpos1[i] - xpos1[i+1])**2 + (ypos1[i] - ypos1[i+1])**2 + (zpos1[i] - zpos1[i+1])**2)**0.5
        if dists1[i] > nspacing:
            conx.append(xpos1[i+1])
            cony.append(ypos1[i+1])
            conz.append(zpos1[i+1])
            nspacing = nspacing + spacing
    dists2[0] = dists1[-1]
    N1 = len(conx)
    for i in range(N1):
        if i%2 == 0:
            sy = cony[i] - pad
            slope = -1*(pa1*2*(conx[i] - ph1))
            inter = conz[i] - slope*conx[i]
            vec = np.array([slope, 1])
            distvec = sum(vec**2)**0.5
            normvec = vec/distvec
            print("i", i,"normvec", normvec)
            sx = conx[i] + np.sign(normvec[0])*(spacing-pad + extra)*(normvec[0])
            sz = conz[i] + np.sign(normvec[0])*(spacing-pad + extra)*normvec[1]
            sidex.append(sx)
            sidey.append(sy)
            sidez.append(sz)

    for i in range(npoints-1):
        dists2[i+1] = dists2[i] + ((xpos2[i] - xpos2[i+1])**2 + (ypos2[i] - ypos2[i+1])**2 + (zpos2[i] - zpos2[i+1])**2)**0.5
        distend = ((xpos2[i+1] - end[0])**2 + (ypos2[i+1] - end[1])**2 + (zpos2[i+1] - end[2])**2)**2
        if dists2[i] > nspacing and distend > spacing:
            conx.append(xpos2[i+1])
            cony.append(ypos2[i+1])
            conz.append(zpos2[i+1])
            nspacing = nspacing + spacing
    N = len(conx)
    for i in range(N1, N):
        if i%2 == 0:
            sy = cony[i] + pad
            slope = -1*(pa2*2*(conx[i] - ph2))
            inter = conz[i] - slope*conx[i]
            vec = np.array([slope, 1])
            distvec = sum(vec**2)**0.5
            normvec = vec/distvec
            sx = conx[i] -  np.sign(normvec[0])*(spacing-pad + extra)*normvec[0]
            sz = conz[i] -  np.sign(normvec[0])*(spacing-pad + extra)*normvec[1]
            sidex.append(sx)
            sidey.append(sy)
            sidez.append(sz)
    return conx, cony, conz, sidex, sidey, sidez


def linear2parabolaConnect(la,lb,pa1, ph1, pk1, pa2, ph2, pk2, start, end,mid,spacing):
    npoints = 1000
    xpos1 = np.linspace(start[0], mid[0], npoints)
    xpos2 = np.linspace(mid[0], end[0], npoints)
    ypos1 = la*xpos1 + lb*np.ones(npoints)
    ypos2 = la*xpos2 + lb*np.ones(npoints)
    zpos1 = pa1*(ypos1 - ph1)**2 + pk1
    zpos2 = pa2*(ypos2 - ph2)**2 + pk2
    dists1 = np.zeros(npoints)
    dists2 = np.zeros(npoints)
    conx = []
    cony = []
    conz = []
    sidex = []
    sidey = []
    sidez = []
    dists1[0] = 0
    pad = 0
    extra = 0
    nspacing = spacing*0.8
    for i in range(npoints-1):
        dists1[i+1] = dists1[i] + ((xpos1[i] - xpos1[i+1])**2 + (ypos1[i] - ypos1[i+1])**2 + (zpos1[i] - zpos1[i+1])**2)**0.5
        if dists1[i] > nspacing:
            conx.append(xpos1[i+1])
            cony.append(ypos1[i+1])
            conz.append(zpos1[i+1])
            nspacing = nspacing + spacing
    dists2[0] = dists1[-1]
    N1 = len(conx)
    #put side groups for first parabola
    for i in range(N1):
        if i%2 == 0:

            sx = conx[i] 

            slope = -1*(pa1*2*(cony[i] - ph1))
            inter = conz[i] - slope*cony[i]
            vec = np.array([slope, 1])
            distvec = sum(vec**2)**0.5
            normvec = vec/distvec
            print("i", i,"normvec", normvec)
            

            sy = cony[i] - np.sign(normvec[0])*(spacing-pad + extra)*(normvec[0])
            sz = conz[i] - np.sign(normvec[0])*(spacing-pad + extra)*normvec[1]
            sidex.append(sx)
            sidey.append(sy)
            sidez.append(sz)

    for i in range(npoints-1):
        dists2[i+1] = dists2[i] + ((xpos2[i] - xpos2[i+1])**2 + (ypos2[i] - ypos2[i+1])**2 + (zpos2[i] - zpos2[i+1])**2)**0.5
        distend = ((xpos2[i+1] - end[0])**2 + (ypos2[i+1] - end[1])**2 + (zpos2[i+1] - end[2])**2)*2
        if dists2[i] > nspacing and distend > spacing:
            conx.append(xpos2[i+1])
            cony.append(ypos2[i+1])
            conz.append(zpos2[i+1])
            nspacing = nspacing + spacing
    N = len(conx)
    for i in range(N1, N):
        if i%2 == 0:
            sx = conx[i]
            slope = -1*(pa2*2*(cony[i] - ph2))
            inter = conz[i] - slope*cony[i]
            vec = np.array([slope, 1])
            distvec = sum(vec**2)**0.5
            normvec = vec/distvec
            sy = cony[i] +  np.sign(normvec[0])*(spacing-pad + extra)*normvec[0]
            sz = conz[i] + np.sign(normvec[0])*(spacing-pad + extra)*normvec[1]
            sidex.append(sx)
            sidey.append(sy)
            sidez.append(sz)

            
    return conx, cony, conz, sidex, sidey, sidez
    


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
    #diffangle = 0
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
    #calculate distances between ends:
    #pairs in unit cells
    ends1s = [[unitcell[8][0]*a, unitcell[8][1]*b, unitcell[8][2]*c], [unitcell[17][0]*a, unitcell[17][1]*b, unitcell[17][2]*c], [unitcell[26][0]*a, unitcell[26][1]*b, unitcell[26][2]*c], [unitcell[35][0]*a, unitcell[35][1]*b, unitcell[35][2]*c]]
    ends2s = [[unitcell[10][0]*a, unitcell[10][1]*b, unitcell[10][2]*c], [unitcell[19][0], unitcell[19][1]*b, unitcell[19][2]*c], [unitcell[28][0]*a, unitcell[28][1]*b, unitcell[28][2]*c], [unitcell[1,0]*a, (unitcell[1,1]+2)*b, unitcell[1,2]*c]]
    midpoints = []
    for i in range(len(ends1s)):
        #print("end1", ends1s[i], "end2", ends2s[i])
        xdist = (ends1s[i][0] - ends2s[i][0])*a
        ydist = (ends1s[i][1] - ends2s[i][1])*b
        zdist = (ends1s[i][2] - ends2s[i][2])*c
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

    #connection 1 
    #parabola in z y
    # linear in x
    Aparabola1 = []
    Hparabola1 = []
    Kparabola1 = []
    Aparabola2 = []
    Hparabola2 = []
    Kparabola2 = []
    Aline = []
    Bline = []

    print('ends1s', ends1s, 'ends1s', ends2s, 'midpoints', midpoints)
    for i in range(len(ends1s)):
        print(ends1s[i][1], midpoints[i][1], ends2s[i][1], ends1s[i][2], midpoints[i][2], ends2s[i][2])
        ph1  = midpoints[i][1]
        pk1  = midpoints[i][2]
        pa1 = (ends1s[i][2] - pk1)/((ends1s[i][1] - ph1)**2)
        Aparabola1.append(pa1)
        Hparabola1.append(ph1)
        Kparabola1.append(pk1)
        ph2 = midpoints[i][1]
        pk2 = midpoints[i][2]
        pa2 = (ends2s[i][2] - pk2)/((ends2s[i][1] - ph2)**2)
        Aparabola2.append(pa2)
        Hparabola2.append(ph2)
        Kparabola2.append(pk2)
        la = (ends2s[i][1] - ends1s[i][1])/(ends2s[i][0] - ends1s[i][0])
        print("slope", la)
        lb = ends2s[i][1] - la*ends2s[i][0]
        #[la,lb] = np.polyfit([ends1s[i][1], ends2s[i][1]],[ends1s[i][0], ends2s[i][0]],1)
        Aline.append(la)
        Bline.append(lb) 
        


    connections1 = []
    connections2 = []
    connections3 = []
    connections4 = []
    npoints = [8,6,8,6]
    spacing = [1.54,1.54,1.54,1.54]
    for i in range(4):
        #xs1 = np.linspace(ends1s[i][0], midpoints[i][0], npoints[i])
        #ys1 = Aline[i]*xs1 + Bline[i]
        #zs1 = Aparabola1[i]*(ys1 - Hparabola1[i])**2 + Kparabola1[i]
        #xs2 = np.linspace(midpoints[i][0], ends2s[i][0], npoints[i])
        #ys2 = Aline[i]*xs2 + Bline[i]
        #zs2 = Aparabola2[i]*(ys2 - Hparabola2[i])**2 + Kparabola2[i]
        conx, cony, conz, sidex, sidey, sidez = linear2parabolaConnect(Aline[i],Bline[i],Aparabola1[i], Hparabola1[i], Kparabola1[i], Aparabola2[i], Hparabola2[i], Kparabola2[i], ends1s[i], ends2s[i],midpoints[i],spacing[i])
        #print("old xs", xs, "old ys", ys, "old zs", zs)
        #xs,ys,zs = paraboliclinear(Aline[i],Bline[i],Aparabola[i],Hparabola[i],Kparabola[i],ends1s[i],ends2s[i],spacing[i]) 
        
        #print("xs", xs, "ys", ys, "zs", zs)
        sidec = 0
        for j in range(len(conx)):
            vars()["connections"+str(i+1)].append([conx[j], cony[j], conz[j]])
            if j%2 == 0:
                vars()["connections"+str(i+1)].append([sidex[sidec], sidey[sidec], sidez[sidec]])
                sidec += 1
        #add side groups


    #planes connection
    
    pend1 = ends1s[3]
    pend2 = [ends1s[3][0]+2*a, ends1s[3][1], ends1s[3][2]]
    pla = (pend2[1] - pend1[1])/(pend2[0] - pend1[0])
    plb = pend2[1] - pla*pend2[0] 
    pmid = [(pend1[0] + pend2[0])/2, (pend1[1] + pend2[1])/2 , (pend1[2] + pend2[2])/2 - c]
    pendh1 = pmid[0]
    pendh2 = pmid[0]
    pendk1 = pmid[2]
    pendk2 = pmid[2]

    penda1 = (pend1[2] - pendk1)/((pend1[0] - pendh1)**2)
    penda2 = (pend2[2] - pendk2)/((pend2[0] - pendh2)**2)
    pconx, pcony, pconz, psidex, psidey, psidez = linear2parabolaXZConnect(pla,plb,penda1, pendh1, pendk1, penda2, pendh2, pendk2, pend1, pend2, pmid,1.6)
    print("pla", pla, "plb", plb)
    print("pend1", pend1)
    print("pend2", pend2)
    print("pmid", pmid)
    print("penda1", penda1, "penda2", penda2)


    p2end1 = [unitcell[0][0]*a , unitcell[0][1]*b, unitcell[0][2]*c]
    p2end2 = [unitcell[0][0]*a + 2*a, p2end1[1], p2end1[2]]
    pla = (p2end2[1] - p2end1[1])/(p2end2[0] - p2end1[0])
    plb = p2end2[1] - pla*p2end2[0]
    pmid = [(p2end1[0] + p2end2[0])/2, (p2end1[1] + p2end2[1])/2, (p2end1[2] + p2end2[2])/2 - c]
    pendh1 = pmid[0]
    pendh2 = pmid[0]
    pendk1 = pmid[2]
    pendk2 = pmid[2]
    
    penda1 = (p2end1[2] - pendk1)/((p2end1[0] - pendh1)**2)
    penda2 = (p2end2[2] - pendk2)/((p2end2[0] - pendh2)**2)
    print("p2end1", p2end1, "p2end2", p2end2)
    print("penda1", penda1, "penda2", penda2)
    p2conx, p2cony, p2conz,p2sidex, p2sidey, p2sidez = linear2parabolaXZConnect(pla, plb, penda1, pendh1, pendk1, penda2, pendh2, pendk2, p2end1, p2end2, pmid, 1.6)
    print("p2sidex", p2sidex, "p2sidey", p2sidey, "p2sidez", p2sidez)
    print("p2conx", p2conx, "p2cony", p2cony, "p2conz", p2conz)
    pconnections1 = []
    cside = 0
    for i in range(len(pconx)):
        pconnections1.append([pconx[i], pcony[i], pconz[i]])
        if i%2 == 0:
            pconnections1.append([psidex[cside], psidey[cside], psidez[cside]])
            cside +=1

    pconnections2 = []
    cside = 0
    for i in range(len(p2conx)):
        pconnections2.append([p2conx[i], p2cony[i], p2conz[i]])
        if i%2 == 0:
            pconnections2.append([p2sidex[cside], p2sidey[cside], p2sidez[cside]])
            cside +=1

    count = 1
    for i in range(Rx):
        basex = i*a*2
        for j in range(Ry):
            basey = j*b*2
            for s in range(4):
                if s%2 == 0:
                    Rzlist = list(range(Rz))
                else:
                    Rzlist = list(range(Rz-1, -1,-1))
                for k in Rzlist:
                    basez = k*c
                    for ci in range(9):
                        mol = ci//9 + 1
                        #x = basex + (symmetries[s][0]*a + rcoords[ci][0]*a) + xoffset*(symmetries[s][2]*c + rcoords[ci][2]*c)*a
                        x = basex + unitcell[s*9 + ci,0]*a 
                        #print("x", x,"basex", basex, "symmetries[s][0]*a", symmetries[s][0]*a, "symmetries[s][0]", symmetries[s][0])
                        y = basey + unitcell[s*9 + ci,1]*b
                        #print("y", y, "basey", basey, "symmetries[s][1]*b", symmetries[s][1]*b, "symmetries[s][1]", symmetries[s][1])
                        #z = basez + (symmetries[s][2]*c + rcoords[ci][2]*c)*zoffset 
                        z = basez + unitcell[s*9 + ci,2]*c
                        tx = x + z*xoffset
                        ty = y
                        tz = z*zoffset
                        atomsinfo.append([count, s+1, 1, tx, ty, tz])
                        count = count + 1
                    if k == Rzlist[-1] and s%2 == 0:
                        for ii in range(len(vars()["connections"+str(s+1)])):
                            x = basex + vars()["connections" + str(s+1)][ii][0]
                            y = basey + vars()["connections" + str(s+1)][ii][1]
                            z = basez + vars()["connections" + str(s+1)][ii][2]
                            tx = x + z*xoffset
                            ty = y 
                            tz = z*zoffset
                            atomsinfo.append([count, 5+s, 1, tx, ty, tz])
                            count = count + 1
                    if k == Rzlist[-1] and s%2 ==1 and (s != 3 or j < Ry-1):
                        for ii in range(len(vars()["connections"+str(s+1)])):
                            x = basex + vars()["connections" + str(s+1)][ii][0]
                            y = basey + vars()["connections" + str(s+1)][ii][1]
                            z = vars()["connections" + str(s+1)][ii][2]
                            tx = x + z*xoffset
                            ty = y
                            tz = z*zoffset
                            atomsinfo.append([count, 5+s, 1, tx, ty, tz])
                            count = count + 1
                    if k == Rzlist[-1] and s==3 and j == Ry-1 and i%2 ==0:
                        for ii in range(len(pconnections1)):
                            x = basex + pconnections1[ii][0]
                            y = basey + pconnections1[ii][1]
                            z = pconnections1[ii][2]
                            tx = x + z*xoffset
                            ty = y
                            tz = z*zoffset

                            atomsinfo.append([count, 5+s, 1, tx, ty, tz])
                            count = count + 1

                    elif k == Rzlist[-1] and s==3 and j == Ry -1 and i%2 ==1:

                        for ii in range(len(pconnections2)):
                            x = basex  + pconnections2[ii][0]
                            y = pconnections2[ii][1]
                            z = pconnections2[ii][2]
                            tx = x + z*xoffset
                            ty = y
                            tz = z*zoffset
                            atomsinfo.append([count, 5+s, 1, tx, ty, tz])
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

def writeMoleculesinCrystal(xyzfile, moleculeClist, Rx, Ry, Rz):
    shownCs = []
    maxC = 0
    spaceC = 5
    #from moleculeClist
    for m in moleculeClist:
        clist = list(range(maxC, maxC + m))
        shownCs.extend(clist)
        maxC = maxC + m + spaceC

    print("shownCs", shownCs)
    atomsinfo = []
    a = 6.63
    b = 20.78
    c = 6.50
    atype = [1, 1, 1, 1, 1, 1, 1, 1, 1]
    beta = 99.5
    diffangle = 99.5 - 90
    # find unit vector of new c axis
    zoffset = np.cos(math.radians(diffangle))
    xoffset = -np.sin(math.radians(diffangle))
    rcoords = np.zeros((9, 3))
    rcoords[:, 0] = [float(i) for i in "-0.0727 -0.0765 -0.1021 -0.3087 -0.1146 -0.1044 0.2775 0.0872 0.1026".split()]
    rcoords[:, 1] = [float(i) for i in "0.2291 0.1592 0.1602 0.0589 0.0928 0.0854 0.0797 0.1156 0.1221".split()]
    rcoords[:, 2] = [float(i) for i in "0.2004 0.2788 0.5098 0.4941 0.6057 0.8428 0.9260 0.9730 1.2109".split()]
    widthy = max(rcoords[:, 1]) - min(rcoords[:, 1])
    widthx = max(rcoords[:, 0]) - min(rcoords[:, 0])
    spacing = widthy / 2

    xnegmin = min(-rcoords[:, 0])
    xposmin = min(rcoords[:, 0])
    diffxmin = xposmin - xnegmin
    symmetries = np.array([[1, 1, 1], [-1, -1, -1], [1, -1 + 0.5, 1 + 0.5], [-1, 1 + 0.5, -1 + 0.5]])
    symmetriesM = np.array([[1, 1, 1], [-1, -1, -1], [1, 1, 1], [-1, -1, -1]])
    symmetriesA = np.array(
        [[1, spacing, 0], [diffxmin, 1 - spacing, 1], [0, 1 + spacing, 0], [1 + diffxmin, 2 - spacing, 1]])
    unitcell = np.zeros((len(symmetries) * len(rcoords), 3))
    ucount = 0
    for s in range(len(symmetries)):
        print("symmetries", symmetries[s])
        for ci in range(len(rcoords)):
            unitcell[ucount, 1] = rcoords[ci, 1] * symmetriesM[s, 1] + symmetriesA[s, 1]
            unitcell[ucount, 2] = (rcoords[ci, 2] * symmetriesM[s, 2] + symmetriesA[s, 2])
            unitcell[ucount, 0] = rcoords[ci, 0] * symmetriesM[s, 0] + symmetriesA[s, 0]
            ucount += 1
    # calculate distances between ends:
    # pairs in unit cells
    ends1s = [[unitcell[8][0] * a, unitcell[8][1] * b, unitcell[8][2] * c],
              [unitcell[17][0] * a, unitcell[17][1] * b, unitcell[17][2] * c],
              [unitcell[26][0] * a, unitcell[26][1] * b, unitcell[26][2] * c],
              [unitcell[35][0] * a, unitcell[35][1] * b, unitcell[35][2] * c]]
    ends2s = [[unitcell[10][0] * a, unitcell[10][1] * b, unitcell[10][2] * c],
              [unitcell[19][0], unitcell[19][1] * b, unitcell[19][2] * c],
              [unitcell[28][0] * a, unitcell[28][1] * b, unitcell[28][2] * c],
              [unitcell[1, 0] * a, (unitcell[1, 1] + 2) * b, unitcell[1, 2] * c]]
    midpoints = []
    for i in range(len(ends1s)):
        # print("end1", ends1s[i], "end2", ends2s[i])
        xdist = (ends1s[i][0] - ends2s[i][0]) * a
        ydist = (ends1s[i][1] - ends2s[i][1]) * b
        zdist = (ends1s[i][2] - ends2s[i][2]) * c
        dist = (xdist ** 2 + ydist ** 2 + zdist ** 2) ** 0.5
    for i in range(len(ends1s)):
        if i % 2 == 0:
            midpoints.append([(ends1s[i][0] + ends2s[i][0]) / 2, (ends1s[i][1] + ends2s[i][1]) / 2,(ends1s[i][2] + ends2s[i][2]) / 2 + c])
        else:
            midpoints.append([(ends1s[i][0] + ends2s[i][0]) / 2, (ends1s[i][1] + ends2s[i][1]) / 2,  (ends1s[i][2] + ends2s[i][2]) / 2 - c])

    # connection 1
    # parabola in z y
    # linear in x
    Aparabola1 = []
    Hparabola1 = []
    Kparabola1 = []
    Aparabola2 = []
    Hparabola2 = []
    Kparabola2 = []
    Aline = []
    Bline = []

    for i in range(len(ends1s)):
        print(ends1s[i][1], midpoints[i][1], ends2s[i][1], ends1s[i][2], midpoints[i][2], ends2s[i][2])
        ph1 = midpoints[i][1]
        pk1 = midpoints[i][2]
        pa1 = (ends1s[i][2] - pk1) / ((ends1s[i][1] - ph1) ** 2)
        Aparabola1.append(pa1)
        Hparabola1.append(ph1)
        Kparabola1.append(pk1)
        ph2 = midpoints[i][1]
        pk2 = midpoints[i][2]
        pa2 = (ends2s[i][2] - pk2) / ((ends2s[i][1] - ph2) ** 2)
        Aparabola2.append(pa2)
        Hparabola2.append(ph2)
        Kparabola2.append(pk2)
        la = (ends2s[i][1] - ends1s[i][1]) / (ends2s[i][0] - ends1s[i][0])
        lb = ends2s[i][1] - la * ends2s[i][0]
        # [la,lb] = np.polyfit([ends1s[i][1], ends2s[i][1]],[ends1s[i][0], ends2s[i][0]],1)
        Aline.append(la)
        Bline.append(lb)

    connections1 = []
    connections2 = []
    connections3 = []
    connections4 = []
    npoints = [8, 6, 8, 6]
    spacing = [1.54, 1.54, 1.54, 1.54]
    for i in range(4):
        conx, cony, conz, sidex, sidey, sidez = linear2parabolaConnect(Aline[i], Bline[i], Aparabola1[i], Hparabola1[i],
                                                                       Kparabola1[i], Aparabola2[i], Hparabola2[i],
                                                                       Kparabola2[i], ends1s[i], ends2s[i],
                                                                       midpoints[i], spacing[i])
        sidec = 0
        for j in range(len(conx)):
            vars()["connections" + str(i + 1)].append([conx[j], cony[j], conz[j]])
            if j % 2 == 0:
                vars()["connections" + str(i + 1)].append([sidex[sidec], sidey[sidec], sidez[sidec]])
                sidec += 1

    pend1 = ends1s[3]
    pend2 = [ends1s[3][0] + 2 * a, ends1s[3][1], ends1s[3][2]]
    pla = (pend2[1] - pend1[1]) / (pend2[0] - pend1[0])
    plb = pend2[1] - pla * pend2[0]
    pmid = [(pend1[0] + pend2[0]) / 2, (pend1[1] + pend2[1]) / 2, (pend1[2] + pend2[2]) / 2 - c]
    pendh1 = pmid[0]
    pendh2 = pmid[0]
    pendk1 = pmid[2]
    pendk2 = pmid[2]

    penda1 = (pend1[2] - pendk1) / ((pend1[0] - pendh1) ** 2)
    penda2 = (pend2[2] - pendk2) / ((pend2[0] - pendh2) ** 2)
    pconx, pcony, pconz, psidex, psidey, psidez = linear2parabolaXZConnect(pla, plb, penda1, pendh1, pendk1, penda2,
                                                                           pendh2, pendk2, pend1, pend2, pmid, 1.6)

    p2end1 = [unitcell[0][0] * a, unitcell[0][1] * b, unitcell[0][2] * c]
    p2end2 = [unitcell[0][0] * a + 2 * a, p2end1[1], p2end1[2]]
    pla = (p2end2[1] - p2end1[1]) / (p2end2[0] - p2end1[0])
    plb = p2end2[1] - pla * p2end2[0]
    pmid = [(p2end1[0] + p2end2[0]) / 2, (p2end1[1] + p2end2[1]) / 2, (p2end1[2] + p2end2[2]) / 2 - c]
    pendh1 = pmid[0]
    pendh2 = pmid[0]
    pendk1 = pmid[2]
    pendk2 = pmid[2]

    penda1 = (p2end1[2] - pendk1) / ((p2end1[0] - pendh1) ** 2)
    penda2 = (p2end2[2] - pendk2) / ((p2end2[0] - pendh2) ** 2)
    p2conx, p2cony, p2conz, p2sidex, p2sidey, p2sidez = linear2parabolaXZConnect(pla, plb, penda1, pendh1, pendk1,
                                                                                 penda2, pendh2, pendk2, p2end1, p2end2,
                                                                                 pmid, 1.6)
    pconnections1 = []
    cside = 0
    for i in range(len(pconx)):
        pconnections1.append([pconx[i], pcony[i], pconz[i]])
        if i % 2 == 0:
            pconnections1.append([psidex[cside], psidey[cside], psidez[cside]])
            cside += 1

    pconnections2 = []
    cside = 0
    for i in range(len(p2conx)):
        pconnections2.append([p2conx[i], p2cony[i], p2conz[i]])
        if i % 2 == 0:
            pconnections2.append([p2sidex[cside], p2sidey[cside], p2sidez[cside]])
            cside += 1

    count = 1
    for i in range(Rx):
        basex = i * a * 2
        for j in range(Ry):
            basey = j * b * 2
            for s in range(4):
                if s % 2 == 0:
                    Rzlist = list(range(Rz))
                else:
                    Rzlist = list(range(Rz - 1, -1, -1))
                for k in Rzlist:
                    basez = k * c
                    for ci in range(9):
                        mol = ci // 9 + 1
                        x = basex + unitcell[s * 9 + ci, 0] * a
                        y = basey + unitcell[s * 9 + ci, 1] * b
                        z = basez + unitcell[s * 9 + ci, 2] * c
                        tx = x + z * xoffset
                        ty = y
                        tz = z * zoffset
                        if count in shownCs: 
                            atomsinfo.append([count, s + 1, 1, tx, ty, tz])
                        count = count + 1
                    if k == Rzlist[-1] and s % 2 == 0:
                        for ii in range(len(vars()["connections" + str(s + 1)])):
                            x = basex + vars()["connections" + str(s + 1)][ii][0]
                            y = basey + vars()["connections" + str(s + 1)][ii][1]
                            z = basez + vars()["connections" + str(s + 1)][ii][2]
                            tx = x + z * xoffset
                            ty = y
                            tz = z * zoffset
                            if count in shownCs:
                                atomsinfo.append([count, 5 + s, 1, tx, ty, tz])
                            count = count + 1
                    if k == Rzlist[-1] and s % 2 == 1 and (s != 3 or j < Ry - 1):
                        for ii in range(len(vars()["connections" + str(s + 1)])):
                            x = basex + vars()["connections" + str(s + 1)][ii][0]
                            y = basey + vars()["connections" + str(s + 1)][ii][1]
                            z = vars()["connections" + str(s + 1)][ii][2]
                            tx = x + z * xoffset
                            ty = y
                            tz = z * zoffset
                            if count in shownCs:
                                atomsinfo.append([count, 5 + s, 1, tx, ty, tz])
                            count = count + 1
                    if k == Rzlist[-1] and s == 3 and j == Ry - 1 and i % 2 == 0:
                        for ii in range(len(pconnections1)):
                            x = basex + pconnections1[ii][0]
                            y = basey + pconnections1[ii][1]
                            z = pconnections1[ii][2]
                            tx = x + z * xoffset
                            ty = y
                            tz = z * zoffset
                            if count in shownCs:
                                atomsinfo.append([count, 5 + s, 1, tx, ty, tz])
                            count = count + 1

                    elif k == Rzlist[-1] and s == 3 and j == Ry - 1 and i % 2 == 1:

                        for ii in range(len(pconnections2)):
                            x = basex + pconnections2[ii][0]
                            y = pconnections2[ii][1]
                            z = pconnections2[ii][2]
                            tx = x + z * xoffset
                            ty = y
                            tz = z * zoffset
                            if count in shownCs:
                                atomsinfo.append([count, 5 + s, 1, tx, ty, tz])
                            count = count + 1

    atomsinfo = np.array(atomsinfo)
    writexyz(xyzfile, atomsinfo)
    xlo = min(atomsinfo[:, 3]) - 5
    xhi = max(atomsinfo[:, 3]) + 5
    ylo = min(atomsinfo[:, 4]) - 5
    yhi = max(atomsinfo[:, 4]) + 5
    zlo = min(atomsinfo[:, 5]) - 5
    zhi = max(atomsinfo[:, 5]) + 5
    return atomsinfo, xlo, xhi, ylo, yhi, zlo, zhi


def main():
    writeMoleculesinCrystal("MolinCryst.xyz", [100, 100, 100], 4, 4, 4)
    #readlammpsbondsPPctypes("TrialInfa1iPPbonds.data", "TrialInfa1iPPC1type.data")
    #atomsinfo, xlo, xhi, ylo, yhi, zlo, zhi = setupInfiniteSystem("TrialInfa1iPP.xyz", 10, 10, 10)
    
    #msc.writelammpsdatajustatoms("TrialInfa1iPP.data",[xlo,xhi,ylo,yhi,zlo,zhi], [15], len(atomsinfo), atomsinfo)
    #atomsinfo, xlo, xhi, ylo, yhi, zlo, zhi = readxyz("custompp-iso.xyz")
    #msc.writelammpsdatajustatoms("custompp-iso.data", [xlo, xhi, ylo, yhi, zlo, zhi], [14, 1], len(atomsinfo), atomsinfo)
main()
