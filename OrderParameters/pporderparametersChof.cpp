#include <iostream>
#include <string>
#include <cstring>
#include <fstream>
#include <sstream>
#include <vector>
#include <stdlib.h>
#include <stdio.h>
#include <cmath>
#include <unistd.h>
#include <algorithm>
#include <map>
#include <numeric>
#include <assert.h>
#include <time.h>
#include <omp.h>
using namespace std;


struct Atominfo {
	int natoms;
	vector <int> atomnumber;
	vector <int> atomtype;
	vector <double> xcoords;
	vector <double> ycoords;
	vector <double>  zcoords;
	vector <int>  xi;
	vector <int> yi;
	vector <int> zi;
};

struct Bondinfo {
	int nbonds;
	vector <int> bondnumber;
	vector <int> bondtype;
	vector <int> bondfirst;
	vector <int> bondsecond;
};

struct Angleinfo {
	int nangles;
	vector <int> anglenumber;
	vector <int> angletype;
	vector <int> anglefirst;
	vector <int> anglesecond;
	vector <int> anglethird;

};

struct Lmpsabainfo {
	Atominfo ai;
	Bondinfo bi;
	Angleinfo ani;
	int atomtypes;
	int bondtypes;
	int angletypes;
	vector <double> masses;
	vector <int> molecule;
	double xboxlen;
	double yboxlen;
	double zboxlen;
	double xlo;
	double xhi;
	double ylo;
	double yhi;
	double zlo;
	double zhi;
};


struct LmpsTraj {
	vector <int> timesteps;
	vector <int> numbersatoms;
	vector<vector<vector<double> > > boxcoordsall;
	vector <Atominfo> atomsinfoall;


};

struct InfoLmpsTraj {
	vector <double> crystallinities;
	vector <double> orientationP2;
	vector <double> HermansOrientP2;
	vector<vector<vector<vector<double>>>> localP2all;

};

struct InfoLmpsData {
	double crystallinities;
	double orientationP2;
	vector<vector<vector<double>>> localP2;
	double HermansOrientP2;
};




Lmpsabainfo read_finalstructwithatombondangleinfo(string filename, int nlines)
/* filename of lammps data file, nlines doesn't need to be exact but needs to be bigger than actual number of lins*/
{
	ifstream infile;
	infile.open(filename.c_str());
	string line;
	char * pch;
	bool linecomment = true;
	bool linetypeinfo = false;
	bool linebox = false;
	bool linemasses = false;
	bool lineatom = false;
	bool linebond = false;
	bool lineangle = false;
	double boxcoords [3][2];
	int nindex;
	int bindex;
	int anindex;
	int natoms;
	int atomtypes;
	int nbonds;
	int bondtypes;
	int nbtype;
	int nangles;
	int angletypes;
	int nangtype;
	int countatoms = 0;
	int countbonds = 0;
	int countangles = 0;
	Atominfo ai;
	Bondinfo bi;
	Angleinfo ani;
	int nmass;
	Lmpsabainfo lmpi;
	for (int i = 0; i < nlines; i++) {
	//cout << "i " << i << endl;
	vector<string> vec_line;
	getline(infile,line);
	if (line.length() >0){
		//cout << line << endl;
	char lineArr[line.length()];
	strcpy(lineArr,line.c_str());
	pch = strtok(lineArr, " ");
	while(pch != NULL) { vec_line.push_back(pch); pch = strtok(NULL, " " );}
	//cout << vec_line[0] << endl;
	
	
	if (linecomment) {
		linecomment = false;
	}
	if (vec_line.size() >= 1) {
		if (vec_line[0].compare("Masses") == 0)
			linemasses = true;
		if (vec_line[0].compare("Atoms") == 0)
			lineatom = true;
		if (vec_line[0].compare("Bonds") == 0)
			linebond = true;
		if (vec_line[0].compare("Angles") == 0)
			lineangle = true;
	}
	if (vec_line.size() >= 4) {
		if (vec_line[2] == "xlo"){
			boxcoords[0][0] = atof (vec_line[0].c_str());
                        boxcoords[0][1] = atof (vec_line[1].c_str());
                        lmpi.xlo = boxcoords[0][0];
                        lmpi.xhi = boxcoords[0][1];
                        lmpi.xboxlen = boxcoords[0][1] - boxcoords[0][0];
                        //cout << "xboxlen " << lmpi.xboxlen << endl;
                }
		else if (vec_line[2] == "ylo"){
                        boxcoords[1][0] = atof (vec_line[0].c_str());
                        boxcoords[1][1] = atof (vec_line[1].c_str());
                        lmpi.ylo = boxcoords[1][0];
                        lmpi.yhi = boxcoords[1][1];
			lmpi.angletypes = angletypes;
                        lmpi.yboxlen = boxcoords[1][1] - boxcoords[1][0];
                        //cout << "yboxlen " << lmpi.yboxlen << endl;
                }
		else if (vec_line[2] == "zlo"){
                        boxcoords[2][0] = atof (vec_line[0].c_str());
                        boxcoords[2][1] = atof (vec_line[1].c_str());
                        lmpi.zlo = boxcoords[2][0];
                        lmpi.zhi = boxcoords[2][1];
                        linebox = i+1;
                        lmpi.zboxlen = boxcoords[2][1] - boxcoords[2][0];
                        //cout << "zboxlen " << lmpi.zboxlen << endl;
                        }
		}
	if (vec_line.size() >= 2 && vec_line[1] == "atoms"){
		//cout << "atoms" << endl;
                natoms = stoi (vec_line[0]);
                ai.natoms = natoms;
                //cout << "natoms: " << natoms << endl;
                ai.atomnumber.resize(natoms);
                        ai.atomtype.resize(natoms);
                        ai.xcoords.resize(natoms);
                        ai.ycoords.resize(natoms);
                        ai.zcoords.resize(natoms);
                        ai.xi.resize(natoms);
                        ai.yi.resize(natoms);
                        ai.zi.resize(natoms);
                        lmpi.molecule.assign(natoms,0);
        }
	else if (vec_line.size() >= 2 && vec_line[1] == "bonds"){
                nbonds = stoi(vec_line[0]);
                //cout << "nbonds: " << nbonds << endl;
                bi.nbonds = nbonds;
                bi.bondnumber.resize(nbonds);
                bi.bondtype.resize(nbonds);
                bi.bondfirst.resize(nbonds);
                bi.bondsecond.resize(nbonds);

        }
	else if (vec_line.size() >= 2 && vec_line[1] == "angles"){
		nangles = stoi(vec_line[0]);
		//cout << "nangles: " << nangles << endl;
		ani.nangles = nangles;
		ani.anglenumber.resize(nangles);
		ani.angletype.resize(nangles);
		ani.anglefirst.resize(nangles);
		ani.anglesecond.resize(nangles);
		ani.anglethird.resize(nangles); 
	}
	else if (vec_line.size() >= 2 && vec_line[1] == "atom"){
                atomtypes = stoi(vec_line[0]);
                        lmpi.masses.assign(atomtypes,0);
                lmpi.atomtypes = atomtypes;
        }
	else if (vec_line.size() >= 2 && vec_line[1] == "bond"){
		bondtypes = stoi(vec_line[0]);
	}
	else if (vec_line.size() >= 2 && vec_line[1] == "angle"){
		angletypes = stoi(vec_line[0]);
	}
	if (linemasses){
                if (vec_line.size() == 1){
                }
                else {
                        nmass = stoi(vec_line[0]);
                        //cout << "nmass " << nmass;
                        lmpi.masses[(nmass-1)] = atof(vec_line[1].c_str());
                        //cout << " mass " << lmpi.masses[(nmass-1)] << endl;
                        if (nmass == atomtypes){
                        linemasses = false;
                        }
                }
	}
	if (lineatom){
        	//cout << "lineatom true " << line << endl;
                        if (vec_line.size() >5 ){
				countatoms++;
                                nindex = stoi(vec_line[0]) - 1;
                                ai.atomnumber[nindex] = nindex+1;
				
                                //cout << "Atomnumber: " << ai.atomnumber[nindex] << endl;
				
                                ai.atomtype[nindex] = stoi(vec_line[2]);
				
                                //cout << "Atomtype: " << ai.atomtype[nindex] << endl;
				
                                lmpi.molecule[nindex] = stoi(vec_line[1]);
                                ai.xcoords[nindex] = atof (vec_line[3].c_str());
                                ai.ycoords[nindex] = atof (vec_line[4].c_str());
                                ai.zcoords[nindex] = atof (vec_line[5].c_str());
				//assert (ai.xcoords[nindex] >= lmpi.xlo && ai.xcoords[nindex] <= lmpi.xhi);
				//assert (ai.ycoords[nindex] >= lmpi.ylo && ai.ycoords[nindex] <= lmpi.yhi);
				//assert (ai.zcoords[nindex] >= lmpi.zlo && ai.zcoords[nindex] <= lmpi.zhi); 
				//cout << "x original " << vec_line[3] << " y original " << vec_line[4] << " z original " << vec_line[5] << endl;
			        //cout << "Atom coords info" <<  ai.atomnumber[nindex] << " " << ai.atomtype[nindex] << " x: " << ai.xcoords[nindex] << " y: " << ai.ycoords[nindex] << " z: " << ai.zcoords[nindex] << endl;

                                if (vec_line.size() == 9){
                                ai.xi[nindex] = stoi(vec_line[6]);
                                ai.yi[nindex] = stoi(vec_line[7]);
                                ai.zi[nindex] = stoi(vec_line[8]);      }
                                else if (vec_line.size() < 9) {
                                ai.xi[nindex] = 0;
                                ai.yi[nindex] = 0;
                                ai.zi[nindex] = 0;
				}
                                if (countatoms >= natoms){
					//cout << "looked at all atoms" << endl;
                                        lineatom = false;

                                }

                        
                    }
                }

	        
        else if (linebond){
                if (vec_line.size() >= 4) {
                        bindex = stoi(vec_line[0]) - 1;
                        bi.bondnumber[bindex] = bindex +1;
                        bi.bondtype[bindex] = stoi(vec_line[1]);
                        bi.bondfirst[bindex] = stoi(vec_line[2]);
                        bi.bondsecond[bindex] = stoi(vec_line[3]);
                        //cout << "bindex " << bindex << " bondnumber " << bi.bondnumber[bindex] << endl;
                        countbonds++;

                        if (countbonds >= nbonds){
				//cout << "looked at all bonds" << endl;
                                linebond = false;
                        }
                }

        }
	else if (lineangle){
		if (vec_line.size() >= 5){
			anindex = stoi(vec_line[0]) -1;
			ani.anglenumber[anindex] = anindex+1;
			ani.angletype[anindex] = stoi(vec_line[1]);
			ani.anglefirst[anindex] = stoi(vec_line[2]);
			ani.anglesecond[anindex] = stoi(vec_line[3]);
			ani.anglethird[anindex] = stoi(vec_line[4]);
			//cout << " angle index " << anindex << " anglenumber " << ani.anglenumber[anindex] << " first atom " << ani.anglefirst[anindex] << " second atom " << ani.anglesecond[anindex] << " third atom " << ani.anglethird[anindex] << endl;	
			
			countangles++;
			if (countangles >= nangles){
				//cout << "looked at all angles" << endl;
				lineangle = false;
				//nlines = i+1;	
				//cout << "nlines " << nlines; 
			}
		}
	}
	

	}

	}	
	//cout << "out of loop in lammps data file" << endl;
	lmpi.ai = ai;
	//cout << " set ai " <<endl;
	lmpi.bi = bi;
	//cout << " set bi " << endl;
	lmpi.ani = ani;
	//cout << " set ani " << endl;
	lmpi.atomtypes = atomtypes;
	//cout << " set atomtypes " << endl;
	lmpi.bondtypes = bondtypes;	
	//cout << " set bondtypes " << endl;
	lmpi.angletypes = angletypes;
	//cout << " set angletypes " << endl;
	return lmpi;

}

vector <double> pbcdist(double (& coords1) [3], double (& coords2) [3], double boxcoords [3][2]){
	//cout << "pbcdist" << endl;
        double xlo = boxcoords[0][0];
        double xhi = boxcoords[0][1];
        double ylo = boxcoords[1][0];
        double yhi = boxcoords[1][1];
        double zlo = boxcoords[2][0];
        double zhi = boxcoords[2][1];
        double boxlenx = xhi - xlo;
        double boxleny = yhi - ylo;
        double boxlenz = zhi - zlo;
        vector <double> dists;
        double distx = coords1[0] - coords2[0];
        double disty = coords1[1] - coords2[1];
        double distz = coords1[2] - coords2[2];
        while (distx > boxlenx*0.5){
                distx = distx - boxlenx;
        }
        while (distx < -boxlenx*0.5){
                distx = distx + boxlenx;
        }
        while (disty > boxleny*0.5){
                disty = disty - boxleny;
        }
        while (disty < -boxleny*0.5){
                disty = disty + boxleny;
        }
        while (distz > boxlenz*0.5){
                distz = distz - boxlenz;
        }
        while (distz < -boxlenz*0.5){
                distz = distz + boxlenz;
        }
	assert (distx <= boxlenx*0.5 && distx >= -boxlenx*0.5);
	assert (disty <= boxleny*0.5 && disty >= -boxleny*0.5);
	assert (distz <= boxlenz*0.5 && distz >= -boxlenz*0.5);
        dists.push_back(distx);
        dists.push_back(disty);
        dists.push_back(distz);
	//cout << "distx " << distx << " disty " << disty << " distz " << distz << endl;
        return dists;

}

vector <double> makeunitvector(double & i, double &  j, double &  k){
        double magnitude = pow(pow(i,2) + pow(j,2) + pow(k,2),0.5);
        vector <double> v;
        v.push_back(i/magnitude);
        v.push_back(j/magnitude);
        v.push_back(k/magnitude);
        return v;

}

vector <double> localchainorientationvector(double (&bead1coords) [3], double (&bead3coords) [3], double boxcoords [3][2]){
        

	//cout << " x diff " << bead3coords[0] - bead1coords[0] << " y diff " << bead3coords[1] - bead1coords[1] << " z diff " << bead3coords[2] - bead1coords[2] << endl;
        //find pbc distance vector between bead1 and bead3
        vector <double> dists = pbcdist(bead1coords, bead3coords, boxcoords);
	//cout << " boxlenx " << boxlenx << " boxleny " << boxleny << " boxlenz " << boxlenz << endl;
	//cout << " pbc x dist " << dists[0] << " pbc y dist " << dists[1] << " pbc z dist " << dists[2] << endl;
        vector <double> v = makeunitvector(dists[0], dists[1], dists[2]);
	//cout << "unit vector v " << v[0] << " " << v[1] << " " << v[2] << endl;
	//cout << "returning v" << endl;
        return v;
}


double calccosthetaij(vector <double> & vi, vector <double> & vj){
double dotproduct = 0;
double vimagnitude = 0;
double vjmagnitude = 0;


for (int i=0;i<3;i++){
	//cout << " i " << i << " vi[i] " << vi[i] << " vj[i] " << vj[i] << endl;
        dotproduct = dotproduct + vi[i]*vj[i];
        vimagnitude = vimagnitude + pow(vi[i],2);
        vjmagnitude = vjmagnitude + pow(vj[i],2);
}
//cout << "dotproduct" << dotproduct << endl;
vimagnitude = pow(vimagnitude,0.5);
vjmagnitude = pow(vjmagnitude,0.5);
//cout << " vimagnitude " << vimagnitude << endl;
//cout << " vjmagnitude " << vjmagnitude << endl;
double costhetaij = dotproduct/(vimagnitude*vjmagnitude);
//cout << "costhetaij " << costhetaij << endl;
return costhetaij;

}



auto lookatanglescosthetaHermansz(Angleinfo& anglesinfo,   double boxcoords [3][2], vector <double> & xcoords, vector <double> & ycoords, vector <double> & zcoords){
	
	int angles = anglesinfo.nangles;

	vector <double> allcosthetaizaxSqd;
	
	int atomi1;
	
	int atomi3;
	
	vector <double> vi (3);
	vector <double> zi = {0,0,1};


	double xlo = boxcoords[0][0];
        double xhi = boxcoords[0][1];
        double ylo = boxcoords[1][0];
        double yhi = boxcoords[1][1];
        double zlo = boxcoords[2][0];
        double zhi = boxcoords[2][1];
        double boxlenx = xhi - xlo;
        double boxleny = yhi - ylo;
        double boxlenz = zhi - zlo;
	double dists1x, dists1y, dists1z;
	double mag1;
	
	for (int i=0; i < angles; i++){
		atomi1 = anglesinfo.anglefirst[i] - 1;
		
		atomi3 = anglesinfo.anglethird[i] - 1;
		double atomi1coords[3] = {xcoords[atomi1],ycoords[atomi1],zcoords[atomi1]};
		double atomi3coords[3] = {xcoords[atomi3],ycoords[atomi3],zcoords[atomi3]};
		dists1x = atomi1coords[0] - atomi3coords[0];
		dists1y = atomi1coords[1] - atomi3coords[1];
		dists1z = atomi1coords[2] - atomi3coords[2];
		//cout << "atomi1coords[0] " << atomi1coords[0] << " atomi3coords[0] " << atomi3coords[0] << endl;
		
		//cout << "atomi1coords[1] " << atomi1coords[1] << " atomi3coords[1] " << atomi3coords[1] << endl;
		
		//cout << "atomi1coords[2] " << atomi1coords[2] << " atomi3coords[2] " << atomi3coords[2] << endl;
		while (dists1x > boxlenx*0.5){
			dists1x = dists1x - boxlenx;
		}
		while (dists1x < -boxlenx*0.5){
			dists1x = dists1x + boxlenx;
		}
		while (dists1y > boxleny*0.5){
			dists1y = dists1y - boxleny;
		}
		while (dists1y < -boxleny*0.5){
			dists1y = dists1y + boxleny;
		}
		while (dists1z > boxlenz*0.5){
			dists1z = dists1z - boxlenz;
		}
		while (dists1z < -boxlenz*0.5){
			dists1z = dists1z + boxlenz;
		}
		mag1 = pow(pow(dists1x,2) + pow(dists1y,2) + pow(dists1z,2),0.5);
		
		vi[0]= dists1x/mag1;
		vi[1] = dists1y/mag1;
		vi[2] = dists1z/mag1;
		
		double costhetaij = 0;
	

		for (int ii=0;ii<3;ii++){
			//cout << " ii " << i << " vi[ii] " << vi[ii] << " vj[ii] " << vj[ii] << endl;
			costhetaij = costhetaij + vi[ii]*zi[ii];
		}

		allcosthetaizaxSqd.push_back(pow(costhetaij,2));
	}
	return allcosthetaizaxSqd;
	
}

vector<vector<vector<double>>> localP2PreSorted(Lmpsabainfo lmpi, Angleinfo& vectorsinfo, vector<vector<vector<vector<int>>>> sortedVectorsBins){
	int nbinsx = sortedVectorsBins.size();
	int nbinsy = sortedVectorsBins[0].size();
	int nbinsz = sortedVectorsBins[0][0].size();
	
	int atomi1;
	int atomi3;
	int atomj1;
	int atomj3;
	int ni;
	int nj;
	double atomi1coords [3];
	double atomi3coords [3];
	double atomj1coords [3];
	double atomj3coords [3];
	vector <double> vi (3);
	vector <double> vj (3);
	double boxlenx = lmpi.xhi - lmpi.xlo;
	double boxleny = lmpi.yhi - lmpi.ylo;
	double boxlenz = lmpi.zhi - lmpi.zlo;
	vector <double> xcoords = lmpi.ai.xcoords;
	vector <double> ycoords = lmpi.ai.ycoords;
	vector <double> zcoords = lmpi.ai.zcoords;
	
	vector<vector<vector<vector<double>>>> binAllcosthetaijsqrd; 
	binAllcosthetaijsqrd.resize(nbinsx);
	for (int i = 0; i < nbinsx; i++){
		binAllcosthetaijsqrd[i].resize(nbinsy);
		for (int j = 0; j < nbinsy; j++){
			binAllcosthetaijsqrd[i][j].resize(nbinsz);
			for (int k = 0; k < nbinsz; k++){
				binAllcosthetaijsqrd[i][j][k].resize(0);
			}
		}
	}

	double dists1x;
	double dists1y;
	double dists1z;
	double dists2x;
	double dists2y;
	double dists2z;

	for (int i = 0; i < nbinsx; i++){
		for (int j =0; j < nbinsy; j++){
			for (int k = 0; k < nbinsz; k++){
				for (int ii= 0; ii < sortedVectorsBins[i][j][k].size(); ii++){
					ni = sortedVectorsBins[i][j][k][ii];
					//cout << "ni "<< ni << endl;
					for (int jj=ii+1; jj < sortedVectorsBins[i][j][k].size(); jj++){
						nj = sortedVectorsBins[i][j][k][jj];
						//cout << "nj " << nj << endl;
						atomi1 = vectorsinfo.anglefirst[ni];
						atomi3 = vectorsinfo.anglethird[ni];
						atomj1 = vectorsinfo.anglefirst[nj];
						atomj3 = vectorsinfo.anglethird[nj];

						atomi1coords[0] = xcoords[atomi1];
						atomi1coords[1] = ycoords[atomi1];
						atomi1coords[2] = zcoords[atomi1];
						atomi3coords[0] = xcoords[atomi3];
						atomi3coords[1] = ycoords[atomi3];
						atomi3coords[2] = zcoords[atomi3];
						atomj1coords[0] = xcoords[atomj1];
						atomj1coords[1] = ycoords[atomj1];
						atomj1coords[2] = zcoords[atomj1];
						atomj3coords[0] = xcoords[atomj3];
						atomj3coords[1] = ycoords[atomj3];
						atomj3coords[2] = zcoords[atomj3];
						

						dists1x = atomi1coords[0] - atomi3coords[0];
						dists1y = atomi1coords[1] - atomi3coords[1];
						dists1z = atomi1coords[2] - atomi3coords[2];
						//cout << "atomi1coords[0] " << atomi1coords[0] << " atomi3coords[0] " << atomi3coords[0] << endl;
						
						//cout << "atomi1coords[1] " << atomi1coords[1] << " atomi3coords[1] " << atomi3coords[1] << endl;
						
						//cout << "atomi1coords[2] " << atomi1coords[2] << " atomi3coords[2] " << atomi3coords[2] << endl;
						while (dists1x > boxlenx*0.5){
							dists1x = dists1x - boxlenx;
						}
						while (dists1x < -boxlenx*0.5){
							dists1x = dists1x + boxlenx;
						}
						while (dists1y > boxleny*0.5){
							dists1y = dists1y - boxleny;
						}
						while (dists1y < -boxleny*0.5){
							dists1y = dists1y + boxleny;
						}
						while (dists1z > boxlenz*0.5){
							dists1z = dists1z - boxlenz;
						}
						while (dists1z < -boxlenz*0.5){
							dists1z = dists1z + boxlenz;
						}
						
						dists2x = atomj1coords[0] - atomj3coords[0];
						dists2y = atomj1coords[1] - atomj3coords[1];
						dists2z = atomj1coords[2] - atomj3coords[2];

						//cout << "atomj1coords[0] " << atomj1coords[0] << " atomj3coords[0] " << atomj3coords[0] << endl;
						
						//cout << "atomj1coords[1] " << atomj1coords[1] << " atomj3coords[1] " << atomj3coords[1] << endl;
						
						//cout << "atomj1coords[2] " << atomj1coords[2] << " atomj3coords[2] " << atomj3coords[2] << endl;

						while (dists2x > boxlenx*0.5){
							dists2x = dists2x - boxlenx;
						}
						while (dists2x < -boxlenx*0.5){
							dists2x = dists2x + boxlenx;
						}
						while (dists2y > boxleny*0.5){
							dists2y = dists2y - boxleny;
						}
						while (dists2y < -boxleny*0.5){
							dists2y = dists2y + boxleny;
						}
						while (dists2z > boxlenz*0.5){
							dists2z = dists2z - boxlenz;
						}
						while (dists2z < -boxlenz*0.5){
							dists2z = dists2z + boxlenz;
						}

						//cout << "dists2x " << dists2x << " dists2y " << dists2y << " dists2z " << dists2z << endl;

						//cout << " boxlenx " << boxlenx << " boxleny " << boxleny << " boxlenz " << boxlenz << endl;
						//cout << " pbc x dist " << dists[0] << " pbc y dist " << dists[1] << " pbc z dist " << dists[2] << endl;
						//vi = makeunitvector(dists1[0], dists1[1], dists1[2]);
						//vj = makeunitvector(dists2[0], dists2[1], dists2[2]);
						double mag1 = pow(pow(dists1x,2) + pow(dists1y,2) + pow(dists1z,2),0.5);
						double mag2 = pow(pow(dists2x,2) + pow(dists2y,2) + pow(dists2z,2),0.5);

						vi[0]= dists1x/mag1;
						vi[1] = dists1y/mag1;
						vi[2] = dists1z/mag1;
						vj[0] = dists2x/mag2;
						vj[1] = dists2y/mag2;
						vj[2] = dists2z/mag2;
						
						double costhetaij = 0;

						for (int ii=0;ii<3;ii++){
							//cout << " ii " << i << " vi[ii] " << vi[ii] << " vj[ii] " << vj[ii] << endl;
							costhetaij = costhetaij + vi[ii]*vj[ii];
						}

						binAllcosthetaijsqrd[i][j][k].push_back(pow(costhetaij,2));
						//cout << " costhetaij "<< costhetaij << endl;
								

					}
				}
			}
		}
	}


	//find average P2 for each bin 

	vector<vector<vector<double>>> avgP2bin;
	avgP2bin.resize(nbinsx);
	for (int i = 0; i < nbinsx; i++){
		avgP2bin[i].resize(nbinsy);
		for (int j =0; j < nbinsy; j++){
			avgP2bin[i][j].resize(nbinsz);
		}
	}

	double avgCosThetaSqrd;

	for (int i = 0; i < nbinsx; i++){
		for (int j = 0; j < nbinsy; j++){
			for (int k = 0; k < nbinsz; k++){
				avgCosThetaSqrd = accumulate(binAllcosthetaijsqrd[i][j][k].begin(), binAllcosthetaijsqrd[i][j][k].end(), 0.0)/binAllcosthetaijsqrd[i][j][k].size();
				avgP2bin[i][j][k] = ((3.0/2.0)*avgCosThetaSqrd - 1.0/2.0);
				if (binAllcosthetaijsqrd[i][j][k].size() > 1){
				//cout << "avgP2bin " << avgP2bin[i][j][k] << " nvectors " << binAllcosthetaijsqrd[i][j][k].size() << endl; 
				}
			}
		}
	}


	return avgP2bin;





}

vector<vector<vector<double>>> localcos2HOFZPreSorted(Lmpsabainfo lmpi, Angleinfo& vectorsinfo, vector<vector<vector<vector<int>>>> sortedVectorsBins, vector<double> zVector){
	int nbinsx = sortedVectorsBins.size();
	int nbinsy = sortedVectorsBins[0].size();
	int nbinsz = sortedVectorsBins[0][0].size();
	
	int atomi1;
	int atomi3;
	int atomj1;
	int atomj3;
	int ni;
	int nj;
	double atomi1coords [3];
	double atomi3coords [3];
	double atomj1coords [3];
	double atomj3coords [3];
	vector <double> vi (3);
	vector <double> vj (3);
	double boxlenx = lmpi.xhi - lmpi.xlo;
	double boxleny = lmpi.yhi - lmpi.ylo;
	double boxlenz = lmpi.zhi - lmpi.zlo;
	vector <double> xcoords = lmpi.ai.xcoords;
	vector <double> ycoords = lmpi.ai.ycoords;
	vector <double> zcoords = lmpi.ai.zcoords;
	
	vector<vector<vector<vector<double>>>> binAllcosthetaijsqrd; 
	binAllcosthetaijsqrd.resize(nbinsx);
	for (int i = 0; i < nbinsx; i++){
		binAllcosthetaijsqrd[i].resize(nbinsy);
		for (int j = 0; j < nbinsy; j++){
			binAllcosthetaijsqrd[i][j].resize(nbinsz);
			for (int k = 0; k < nbinsz; k++){
				binAllcosthetaijsqrd[i][j][k].resize(0);
			}
		}
	}

	double dists1x;
	double dists1y;
	double dists1z;
	double dists2x;
	double dists2y;
	double dists2z;

	for (int i = 0; i < nbinsx; i++){
		for (int j =0; j < nbinsy; j++){
			for (int k = 0; k < nbinsz; k++){
				for (int ii= 0; ii < sortedVectorsBins[i][j][k].size(); ii++){
					ni = sortedVectorsBins[i][j][k][ii];
					atomi1 = vectorsinfo.anglefirst[ni];
					atomi3 = vectorsinfo.anglethird[ni];
					atomj1 = vectorsinfo.anglefirst[nj];
					atomj3 = vectorsinfo.anglethird[nj];

					atomi1coords[0] = xcoords[atomi1];
					atomi1coords[1] = ycoords[atomi1];
					atomi1coords[2] = zcoords[atomi1];
					atomi3coords[0] = xcoords[atomi3];
					atomi3coords[1] = ycoords[atomi3];
					atomi3coords[2] = zcoords[atomi3];
						

					dists1x = atomi1coords[0] - atomi3coords[0];
					dists1y = atomi1coords[1] - atomi3coords[1];
					dists1z = atomi1coords[2] - atomi3coords[2];
					//cout << "atomi1coords[0] " << atomi1coords[0] << " atomi3coords[0] " << atomi3coords[0] << endl;
						
					//cout << "atomi1coords[1] " << atomi1coords[1] << " atomi3coords[1] " << atomi3coords[1] << endl;
						
					//cout << "atomi1coords[2] " << atomi1coords[2] << " atomi3coords[2] " << atomi3coords[2] << endl;
					while (dists1x > boxlenx*0.5){
						dists1x = dists1x - boxlenx;
					}
					while (dists1x < -boxlenx*0.5){
						dists1x = dists1x + boxlenx;
					}
					while (dists1y > boxleny*0.5){
						dists1y = dists1y - boxleny;
					}
					while (dists1y < -boxleny*0.5){
						dists1y = dists1y + boxleny;
					}
					while (dists1z > boxlenz*0.5){
						dists1z = dists1z - boxlenz;
					}
					while (dists1z < -boxlenz*0.5){
						dists1z = dists1z + boxlenz;
					}
						
					
					//vi = makeunitvector(dists1[0], dists1[1], dists1[2]);
					//vj = makeunitvector(dists2[0], dists2[1], dists2[2]);
					double mag1 = pow(pow(dists1x,2) + pow(dists1y,2) + pow(dists1z,2),0.5);
					
					double mag2 = pow(pow(zVector[0],2) + pow(zVector[1],2) + pow(zVector[2],2), 0.5);

					vi[0]= dists1x/mag1;
					vi[1] = dists1y/mag1;
					vi[2] = dists1z/mag1;
					vj[0] = zVector[0]/mag2;
					vj[1] = zVector[1]/mag2;
					vj[2] = zVector[2]/mag2;
						
					double costhetaij = 0;

					for (int ii=0;ii<3;ii++){
						//cout << " ii " << i << " vi[ii] " << vi[ii] << " vj[ii] " << vj[ii] << endl;
						costhetaij = costhetaij + vi[ii]*vj[ii];
					}

					binAllcosthetaijsqrd[i][j][k].push_back(pow(costhetaij,2));
					//cout << " costhetaij "<< costhetaij << endl;
								

					
				}
			}
		}
	}


	//find average P2 with respect to z axis

	vector<vector<vector<double>>> avgcos2bin;
	avgcos2bin.resize(nbinsx);
	for (int i = 0; i < nbinsx; i++){
		avgcos2bin[i].resize(nbinsy);
		for (int j =0; j < nbinsy; j++){
			avgcos2bin[i][j].resize(nbinsz);
		}
	}

	double avgCosThetaSqrd;

	for (int i = 0; i < nbinsx; i++){
		for (int j = 0; j < nbinsy; j++){
			for (int k = 0; k < nbinsz; k++){
				avgCosThetaSqrd = accumulate(binAllcosthetaijsqrd[i][j][k].begin(), binAllcosthetaijsqrd[i][j][k].end(), 0.0)/binAllcosthetaijsqrd[i][j][k].size();
				avgcos2bin[i][j][k] = avgCosThetaSqrd;
				//avgP2bin[i][j][k] = ((3.0/2.0)*avgCosThetaSqrd - 1.0/2.0);
				if (binAllcosthetaijsqrd[i][j][k].size() > 1){
				//cout << "avgP2bin " << avgP2bin[i][j][k] << " nvectors " << binAllcosthetaijsqrd[i][j][k].size() << endl; 
				}
			}
		}
	}


	return avgcos2bin;





}

auto lookatanglescostheta(Angleinfo& anglesinfo,  double boxcoords [3][2], vector <double> & xcoords, vector <double> & ycoords, vector <double> & zcoords, int & maxVecSize, int & ii){

	int angles = anglesinfo.nangles;
	
	
	vector <double> allcosthetaijsqrd;

	int atomi1;
	int atomi2;
	int atomi3;
	int atomj1;
	int atomj2;
	int atomj3;
	vector <double> vi (3);
	vector <double> vj (3);
	double dists1x;
	double dists1y;
	double dists1z;
	double dists2x;
	double dists2y;
	double dists2z;
	double atomi1coords[3];
	double atomi3coords[3];
	double atomj1coords[3];
	double atomj3coords[3];
	double costhetaij;
	double xlo = boxcoords[0][0];
        double xhi = boxcoords[0][1];
        double ylo = boxcoords[1][0];
        double yhi = boxcoords[1][1];
        double zlo = boxcoords[2][0];
        double zhi = boxcoords[2][1];
        double boxlenx = xhi - xlo;
        double boxleny = yhi - ylo;
        double boxlenz = zhi - zlo;
	double mag1;
	double mag2;
	int count = 0;
	int next = 0;
	for (int i=ii; i < angles; ++i){
		//cout << "i " << i << " angles " << angles << " count " << count << endl;
		next = (angles - i+1);
		if (allcosthetaijsqrd.size() + next < maxVecSize){ 	
			for (int j=i+1; j < angles; ++j){
				//cout << "i" << i << " j " << j << " angles " << angles << endl;
				
				atomi1 = anglesinfo.anglefirst[i] - 1;
				atomi3 = anglesinfo.anglethird[i] -1;
				atomj1 = anglesinfo.anglefirst[j] - 1;
				atomj3 = anglesinfo.anglethird[j] - 1;
				//cout << "atomi1 " << atomi1 +1 << " atomi3 " << atomi3 + 1 << " atomj1 " << atomj1 << " atomj3 " << atomj3 << endl;
				
				atomi1coords[0] = xcoords[atomi1];
				atomi1coords[1] = ycoords[atomi1];
				atomi1coords[2] = zcoords[atomi1];
				atomi3coords[0] = xcoords[atomi3];
				atomi3coords[1] = ycoords[atomi3];
				atomi3coords[2] = zcoords[atomi3];
				atomj1coords[0] = xcoords[atomj1];
				atomj1coords[1] = ycoords[atomj1];
				atomj1coords[2] = zcoords[atomj1];
				atomj3coords[0] = xcoords[atomj3];
				atomj3coords[1] = ycoords[atomj3];
				atomj3coords[2] = zcoords[atomj3];
				//atomi1coords= {atomsinfo.xcoords[atomi1],atomsinfo.ycoords[atomi1],atomsinfo.zcoords[atomi1]};
				//atomi3coords = {atomsinfo.xcoords[atomi3],atomsinfo.ycoords[atomi3],atomsinfo.zcoords[atomi3]};
				//atomj1coords = {atomsinfo.xcoords[atomj1],atomsinfo.ycoords[atomj1], atomsinfo.zcoords[atomj1]};
				//atomj3coords = {atomsinfo.xcoords[atomj3],atomsinfo.ycoords[atomj3], atomsinfo.zcoords[atomj3]};
				
				//vi = localchainorientationvector(atomi1coords,atomi3coords,boxcoords);
				//vj = localchainorientationvector(atomj1coords,atomj3coords,boxcoords);
				//dists1 = pbcdist(atomi1coords, atomi3coords, boxcoords);
				//dists2 = pbcdist(atomj1coords, atomj3coords, boxcoords);

				dists1x = atomi1coords[0] - atomi3coords[0];
				dists1y = atomi1coords[1] - atomi3coords[1];
				dists1z = atomi1coords[2] - atomi3coords[2];
				//cout << "atomi1coords[0] " << atomi1coords[0] << " atomi3coords[0] " << atomi3coords[0] << endl;
				
				//cout << "atomi1coords[1] " << atomi1coords[1] << " atomi3coords[1] " << atomi3coords[1] << endl;
				
				//cout << "atomi1coords[2] " << atomi1coords[2] << " atomi3coords[2] " << atomi3coords[2] << endl;
				while (dists1x > boxlenx*0.5){
					dists1x = dists1x - boxlenx;
				}
				while (dists1x < -boxlenx*0.5){
					dists1x = dists1x + boxlenx;
				}
				while (dists1y > boxleny*0.5){
					dists1y = dists1y - boxleny;
				}
				while (dists1y < -boxleny*0.5){
					dists1y = dists1y + boxleny;
				}
				while (dists1z > boxlenz*0.5){
					dists1z = dists1z - boxlenz;
				}
				while (dists1z < -boxlenz*0.5){
					dists1z = dists1z + boxlenz;
				}
				
				dists2x = atomj1coords[0] - atomj3coords[0];
				dists2y = atomj1coords[1] - atomj3coords[1];
				dists2z = atomj1coords[2] - atomj3coords[2];

				//cout << "atomj1coords[0] " << atomj1coords[0] << " atomj3coords[0] " << atomj3coords[0] << endl;
				
				//cout << "atomj1coords[1] " << atomj1coords[1] << " atomj3coords[1] " << atomj3coords[1] << endl;
				
				//cout << "atomj1coords[2] " << atomj1coords[2] << " atomj3coords[2] " << atomj3coords[2] << endl;

				while (dists2x > boxlenx*0.5){
					dists2x = dists2x - boxlenx;
				}
				while (dists2x < -boxlenx*0.5){
					dists2x = dists2x + boxlenx;
				}
				while (dists2y > boxleny*0.5){
					dists2y = dists2y - boxleny;
				}
				while (dists2y < -boxleny*0.5){
					dists2y = dists2y + boxleny;
				}
				while (dists2z > boxlenz*0.5){
					dists2z = dists2z - boxlenz;
				}
				while (dists2z < -boxlenz*0.5){
					dists2z = dists2z + boxlenz;
				}

				//cout << "dists2x " << dists2x << " dists2y " << dists2y << " dists2z " << dists2z << endl;

				//cout << " boxlenx " << boxlenx << " boxleny " << boxleny << " boxlenz " << boxlenz << endl;
				//cout << " pbc x dist " << dists[0] << " pbc y dist " << dists[1] << " pbc z dist " << dists[2] << endl;
				//vi = makeunitvector(dists1[0], dists1[1], dists1[2]);
				//vj = makeunitvector(dists2[0], dists2[1], dists2[2]);
				mag1 = pow(pow(dists1x,2) + pow(dists1y,2) + pow(dists1z,2),0.5);
				mag2 = pow(pow(dists2x,2) + pow(dists2y,2) + pow(dists2z,2),0.5);

				vi[0]= dists1x/mag1;
				vi[1] = dists1y/mag1;
				vi[2] = dists1z/mag1;
				vj[0] = dists2x/mag2;
				vj[1] = dists2y/mag2;
				vj[2] = dists2z/mag2;
				
				costhetaij = 0;

				for (int ii=0;ii<3;ii++){
					//cout << " ii " << i << " vi[ii] " << vi[ii] << " vj[ii] " << vj[ii] << endl;
					costhetaij = costhetaij + vi[ii]*vj[ii];
				}
				
				//cout << "costhetaij " << costhetaij << endl;
				allcosthetaijsqrd.push_back((pow(costhetaij,2)));
				count++;
				
				
			} 
		}
		else {
			ii = i;
			return allcosthetaijsqrd;
		}
	}
	return allcosthetaijsqrd;
	

}


double calculateP2(double boxcoords [3][2], Atominfo atomsinfo, Angleinfo anglesinfo){
	double nTerms = anglesinfo.nangles -1;
	cout << "angles " << anglesinfo.nangles << " nTerms " << nTerms << endl;
	double totalTerms = (nTerms/2.0)*(anglesinfo.nangles);
	int tTerms = int(totalTerms);
	cout << " totalTerms " << totalTerms << endl;
	int maxVecSize = 2000000000;
	cout << "maxVecSize " << maxVecSize << endl;
	int nVectors = totalTerms/maxVecSize + 1;
	int ii = 0;
	cout << "nVectors " << nVectors << endl;
	
	vector <double> allcosthetaijsqrd;
	double average; 
        for (int i = 0; i < nVectors; i++){
		cout << "i " << i << endl;
		allcosthetaijsqrd = (lookatanglescostheta(anglesinfo, boxcoords, atomsinfo.xcoords, atomsinfo.ycoords, atomsinfo.zcoords, maxVecSize, ii));
		average += accumulate(allcosthetaijsqrd.begin(), allcosthetaijsqrd.end(),0.0);
	}
      
	average = average/totalTerms;
	//accumulate(allcosthetaijsqrd.begin(),allcosthetaijsqrd.end(),0.0)/allcosthetaijsqrd.size();
	//cout << "average: " << average << endl;
        double P2 = (3.0*average - 1.0)/2.0;
	cout << "orientation P2 " << P2 << endl;
        return P2;
}







double calculateHermansP2orientation(double boxcoords [3][2], Atominfo atomsinfo, Angleinfo anglesinfo){
	// Calculate angle between 1,3 and [0,0,1] vector not between pairs of angles
	vector <double> allcosthetaijsqrd = lookatanglescosthetaHermansz(anglesinfo, boxcoords, atomsinfo.xcoords, atomsinfo.ycoords, atomsinfo.zcoords);
	
	double average  = accumulate(allcosthetaijsqrd.begin(), allcosthetaijsqrd.end(),0.0)/allcosthetaijsqrd.size();
	//cout << "average: " << average << endl;
	double HermansP2 = (3.0*average -1.0)/2.0;
	//cout << "HermansOrientation P2 " << HermansP2 << endl;
	return HermansP2;
	
}

int lastindexlessthanval(vector <double> boundvals, double value){
	int index = 0;
	for (int i = 0; i < boundvals.size()-1; i++){
		if (boundvals[i] <= value){
			index = i;
		}

	} 
	return index;

}

vector <int> squarewellconvbin(double coords [3], vector <double> xbounds, vector <double> ybounds, vector <double> zbounds, double resolutions [3]){
        double x = coords[0];
        double y = coords[1];
        double z = coords[2];
	//cout << "squarewellconvbin " << " ybounds.size() " << ybounds.size() << endl;
        //int xbin = upper_bound(xbounds.begin(),xbounds.end(),x) -xbounds.begin() - 1;
	int xbin = lastindexlessthanval(xbounds, x);

        //int ybin = upper_bound(ybounds.begin(),ybounds.end(),y) -ybounds.begin() - 1;
	int ybin = lastindexlessthanval(ybounds,y);
        //int zbin = upper_bound(zbounds.begin(),zbounds.end(),z) -zbounds.begin() - 1;
	int zbin = lastindexlessthanval(zbounds,z);
	//cout << "xbin: " << xbin << " ybin: " << ybin << " zbin: " << zbin << endl;
	//cout << " xbounds[0] " << xbounds[0] << " xbounds[xbounds.size()-1] " << xbounds[xbounds.size()-1] << endl;
	//cout << " ybounds[0] " << ybounds[0] << " ybounds[xbounds.size()-1] " << ybounds[ybounds.size()-1]
//<< endl;
	//cout << " zbounds[0] " << zbounds[0] << " zbounds[zbounds.size()-1] " << zbounds[zbounds.size()-1]
//<< endl;
	//cout << " x " << x << " xbin " << xbin << " xbounds[xbin] " << xbounds[xbin] << endl;
	//cout << " y " << y << " ybin " << ybin << " ybounds[ybin] " << ybounds[ybin] << endl;
	//cout << " z " << z << " zbin " << zbin << " zbounds[zbin] " << zbounds[zbin] << endl;
        vector <int> binlocs;
	
	//cout << "xbin " << xbin << " ybin " << ybin << " zbin " << zbin << endl;
        binlocs.push_back(xbin);
        binlocs.push_back(ybin);
        binlocs.push_back(zbin);
        return binlocs;


}



vector<vector<vector<vector<double>>>> calcP2localanglescostheta(vector <double> xbounds, vector <double> ybounds, vector <double> zbounds, Angleinfo anglesinfo, Atominfo atomsinfo, double boxcoords [3][2], double resolutions [3]){
        vector <double> vi;
        vector <double> vj;
        int nbinsx = xbounds.size() -1;
        int nbinsy = ybounds.size() -1;
        int nbinsz = zbounds.size() -1;
        //cout << "nbinsx " << nbinsx << " nbinsy " << nbinsy << " nbinsz " << nbinsz << endl;
        vector<vector<vector<vector<double>>>> allcosthetaij;
        //set up sizes (nbinsx*nbinsy*nbinsz)
        allcosthetaij.resize(nbinsx);
	//cout << "resized allcosthetaij x " << allcosthetaij.size() << endl;
        for (int i=0; i <  nbinsx; i++){
                allcosthetaij[i].resize(nbinsy);
		//cout << "resized allcosthetaij y " << allcosthetaij[i].size() << endl;
                for (int j = 0; j < nbinsy; j++){
                        allcosthetaij[i][j].resize(nbinsz);
			//cout << "resized allcosthetaij z " << allcosthetaij[i][j].size() << endl;

			for (int k =0; k < nbinsz; k++){
			allcosthetaij[i][j][k].resize(0);
			//cout << "resized element in i,j,k allcosthetaij" << endl;
			}
			
                }
        }
	//cout << "resized allcosthetaij" << endl;
        int angles = anglesinfo.nangles;
	cout << "angles " << angles << endl;
        for (int i=0; i < angles; i++){
                for (int j=i+1; j < angles; j++){
			//cout << " i " << i << " j " << j << endl;
                        //choose the bins based on the coordinates

                        int atomi1 = anglesinfo.anglefirst[i] - 1;
                        int atomi2 = anglesinfo.anglesecond[i] -1;
                        int atomi3 = anglesinfo.anglethird[i] -1;
                        int atomj1 = anglesinfo.anglefirst[j] - 1;
                        int atomj2 = anglesinfo.anglesecond[j] -1;
                        int atomj3 = anglesinfo.anglethird[j] - 1;
			//cout << "atomi1 " << atomi1 << " atomj1  " << atomj1 << endl;
                        double atomi1coords[3] = {atomsinfo.xcoords[atomi1],atomsinfo.ycoords[atomi1],atomsinfo.zcoords[atomi1]};
                        double atomi2coords[3] = {atomsinfo.xcoords[atomi2],atomsinfo.ycoords[atomi2],atomsinfo.zcoords[atomi2]};
                        double atomi3coords[3] = {atomsinfo.xcoords[atomi3],atomsinfo.ycoords[atomi3],atomsinfo.zcoords[atomi3]};
                        double atomj1coords[3] = {atomsinfo.xcoords[atomj1],atomsinfo.ycoords[atomj1], atomsinfo.zcoords[atomj1]};
                        double atomj2coords[3] = {atomsinfo.xcoords[atomj2],atomsinfo.ycoords[atomj2],atomsinfo.zcoords[atomj2]};
                        double atomj3coords[3] = {atomsinfo.xcoords[atomj3],atomsinfo.ycoords[atomj3], atomsinfo.zcoords[atomj3]};
			//cout << "atomi2y " << atomi2coords[1] << " atomj2y " << atomj2coords[1] << endl;
                        vector <int> ibins = squarewellconvbin(atomi2coords,xbounds, ybounds, zbounds, resolutions);
                        vector <int> jbins = squarewellconvbin(atomj2coords,xbounds,ybounds,zbounds, resolutions);
                        //only do a calculation if in same bin
                        if (ibins == jbins){
			//cout << "ibins (" << ibins[0] << "," << ibins[1] << "," << ibins[2] << ") jbins (" << jbins[0] << "," << jbins[1] << "," << jbins[2] << ")" << endl;
			//cout << "going to localchainorientationvector function" << endl;
                        vi = localchainorientationvector(atomi1coords,atomi3coords,boxcoords);	
			//cout << "got vi" << endl;
                        vj = localchainorientationvector(atomj1coords,atomj3coords,boxcoords);
			//cout << "got vj" << endl;
			//cout << "before pushing" << endl;
			//cout << "ibins[0] " << ibins[0] << " ibins[1] " << ibins[1] << " ibins[2] " << ibins[2] << endl;
			//cout << " nbinsx " << nbinsx << " nbinsy " << nbinsy << " nbinsz " << nbinsz << endl;
			//cout << "before pushing: allcosthetaij[ibins[0]][ibins[1]][ibins[2]].size()" << allcosthetaij[ibins[0]][ibins[1]][ibins[2]].size() << endl;
			//allcosthetaij[ibins[0]][ibins[1]][ibins[2]].resize(allcosthetaij[ibins[0]][ibins[1]][ibins[2]].size()+1);
                        //allcosthetaij[ibins[0]][ibins[1]][ibins[2]][allcosthetaij[ibins[0]][ibins[1]][ibins[2]].size()-1] = calccosthetaij(vi,vj);
			allcosthetaij[ibins[0]][ibins[1]][ibins[2]].push_back(calccosthetaij(vi,vj));
			//cout << "size allcosthetaij x dimension " << allcosthetaij.size() << endl;
			//cout << "size allcosthetaij y dimension " << allcosthetaij[0].size() << endl;
			//cout << "size allcosthetaij z dimension " << allcosthetaij[0][0].size() << endl;
			//cout << "size allcosthetaij fourth dimension here "  << allcosthetaij[ibins[0]][ibins[1]][ibins[2]].size() << endl;
			//cout << "costhetavalue" << calccosthetaij(vi,vj) << endl;
                        }

                }
        }

        return allcosthetaij;
}


vector<vector<vector<vector<double>>>> calcmassesinbins(vector <double> xbounds, vector <double> ybounds, vector <double> zbounds, Lmpsabainfo lmpi, Atominfo atomsinfo, double resolutions [3]){
        int nbinsx = xbounds.size() -1;
        int nbinsy = ybounds.size() -1;
        int nbinsz = zbounds.size() -1;
        //cout << "nbinsx " << nbinsx << " nbinsy " << nbinsy << " nbinsz " << nbinsz << endl;
        vector<vector<vector<vector<double>>>> allmassesinbins;
        //set up sizes (nbinsx*nbinsy*nbinsz)
        allmassesinbins.resize(nbinsx);
        for (int i=0; i <  nbinsx; i++){
                allmassesinbins[i].resize(nbinsy);
                for (int j = 0; j < nbinsy; j++){
                        allmassesinbins[i][j].resize(nbinsz);

                        for (int k =0; k < nbinsz; k++){
                        allmassesinbins[i][j][k].resize(0);
                        }

                }
        }
        int natoms = atomsinfo.natoms;
	//cout << "natoms: " << natoms << endl;
        for (int i=0; i < natoms; i++){
                        //choose the bins based on the coordinates
			//cout << "atomindex: " << i << endl;
                        double atomcoords[3] = {atomsinfo.xcoords[i],atomsinfo.ycoords[i],atomsinfo.zcoords[i]};
			//cout << " x: " << atomcoords[0] << " y: " << atomcoords[1] << " z: " << atomcoords[2] << endl;
                        vector <int> bin = squarewellconvbin(atomcoords,xbounds,ybounds,zbounds,resolutions);
			//cout << "bin i: " << bin[0] << " bin j: " << bin[1] << " bin k: " << bin[2] << endl;
                        int atomtype = atomsinfo.atomtype[i];
                        double mass = lmpi.masses[atomtype-1];
                        allmassesinbins[bin[0]][bin[1]][bin[2]].push_back(mass);


        }

        return allmassesinbins;


}

vector<vector<vector<double>>> localDensity(double boxcoords[3][2], Atominfo atomsinfo,Lmpsabainfo lmpi, vector <int> coordinateindices, double resolutions [3]){
	double xlo = boxcoords[0][0];
	double xhi = boxcoords[0][1];
	double ylo = boxcoords[1][0];
	double yhi = boxcoords[1][1];
	double zlo = boxcoords[2][0];
	double zhi = boxcoords[2][1];
	double lx = xhi -xlo;
	double ly = yhi - ylo;
	double lz = zhi - zlo;
	//set up for local density array
	int ndim = coordinateindices.size();
	assert (ndim==1 || ndim==2 || ndim ==3);
	vector<int> nbins(3,1);
	vector<double> xbounds = {xlo,xhi};
	vector<double> ybounds = {ylo,yhi};
	vector <double> zbounds = {zlo,zhi};
	for (int i=0; i < ndim; ++i){
		int c = coordinateindices[i];
		if (c==0){
                        xbounds[0] = xlo;
                        //cout <<  "xlo " << xlo << endl;
                        xbounds[1] = xlo+resolutions[0];
                        //cout << "xlo+resolutions[0] " << xlo+resolutions[0] << endl;
                        double xval = xlo+resolutions[0];
                        //cout << "xval " << xval << endl;
                        while (xval+resolutions[0] < xhi+resolutions[0]){
                                //cout << "xval " << xval << endl;
                                xval = xval + resolutions[0];

                                xbounds.push_back(xval);
                        }
                        //cout << "len(xbounds) " << xbounds.size() << endl;
                        for (int ii=0; ii < xbounds.size(); ii++){
                        //cout << "xbounds ii " << ii << " value " << xbounds[ii] << endl;
                        }
                }
                else if (c==1){
                        ybounds[0] = ylo;
                        //cout << "ylo " << ylo << endl;
                        //cout << "yhi " << yhi << endl;
                        ybounds[1] = ylo+resolutions[1];
                        //cout << "ylo+resolutions[1] " << ylo+resolutions[1] << endl;
                        double yval = ylo+resolutions[1];
                        //cout << "yhi+resolutions[1] " << yhi+resolutions[1] << endl;
                        while (yval+resolutions[1] < yhi+resolutions[1]){
                                //cout << "yval " << yval << endl;
                                yval = yval + resolutions[1];
                                ybounds.push_back(yval);
                        }
                        //cout << "yval " << yval << endl;
                        //cout << "len(ybounds) " << ybounds.size() << endl;

                }
                else if (c==2){
                        zbounds[0] = zlo;
                        //cout << "zlo " << zlo << endl;
                        zbounds[1] = zlo+resolutions[2];
                        //cout << "zlo+resolutions[0] " << zlo+resolutions[0] << endl;
                        double zval = zlo+resolutions[2];
                        while (zval+resolutions[2] < zhi+resolutions[2]){
                                //cout << "zval " << zval << endl;
                                zval = zval + resolutions[2];

                                zbounds.push_back(zval);
                        }
                        //cout << "len(zbounds) " << zbounds.size() << endl;
                }
                else {
                cout << "Error in coordinateindicies input" << endl;
                }
        }
	//find masses in each bin, then divide by volume of bin to get lacal density
        vector<vector<vector<vector<double>>>> allmassesinbins =  calcmassesinbins(xbounds, ybounds, zbounds,lmpi, atomsinfo,resolutions);
	//calculate P2 for each bin
        int nbinsx = xbounds.size() - 1;
        int nbinsy = ybounds.size() - 1;
        int nbinsz = zbounds.size() - 1;

        assert (nbinsx == allmassesinbins.size());
        assert (nbinsy == allmassesinbins[0].size());
        assert (nbinsz == allmassesinbins[0][0].size());

        vector<vector<vector<double>>> localDensity;
        localDensity.resize(nbinsx);
        for (int i=0;   i < nbinsx; i++){
               localDensity[i].resize(nbinsy);
                for (int j=0; j < nbinsy; j++){
                        localDensity[i][j].resize(nbinsz);
                }
        }
        vector <double> massesinbin;
	double xwidth;
	double ywidth;
	double zwidth;
        for (int i=0; i < nbinsx; i++){
                for (int j=0; j < nbinsy; j++){
                        for (int k=0; k < nbinsz; k++){
                                massesinbin = allmassesinbins[i][j][k];
                                for (int iii = 0; iii < massesinbin.size(); iii++){
                                }
                                if (massesinbin.size() >=1){
                                        double totalmass = accumulate(massesinbin.begin(),massesinbin.end(),0.0);
					if (xbounds[i+1] <= xhi){
					xwidth = xbounds[i+1] - xbounds[i];
					}
					else{
					xwidth = xhi - xbounds[i];
					}
					if (ybounds[j+1] <= yhi){
					ywidth = ybounds[j+1] - ybounds[j];
					}
					else{
					ywidth = yhi - ybounds[j];
					}
					if (zbounds[k+1] <= zhi){
					zwidth = zbounds[k+1] - zbounds[k];
					}
					else{
					zwidth = zhi - zbounds[k];
					}
					double volA = xwidth*ywidth*zwidth;
					double convertgmolAtogmL = (pow(10,30)/(6.022140857*pow(10,23)))/(1000*1000);
                                        localDensity[i][j][k] = (totalmass/volA)*convertgmolAtogmL;
                                        //cout << "P2local[i][j][k] " << P2local[i][j][k] << " i " << i << " j " << j << " k " << k << endl;

                                }
                                else{
                                        localDensity[i][j][k] = std::numeric_limits<double>::quiet_NaN();
                                        //cout << "no need to calculate P2local here" << endl;
                                }

                        }
                }
        }

        return localDensity;
	

}


vector<vector<vector<double>>> localP2(double boxcoords[3][2], Atominfo atomsinfo, Angleinfo anglesinfo, vector <int> coordinateindices, double resolutions [3]){
	double xlo = boxcoords[0][0];
        double xhi = boxcoords[0][1];
        double ylo = boxcoords[1][0];
        double yhi = boxcoords[1][1];
        double zlo = boxcoords[2][0];
        double zhi = boxcoords[2][1];
	//cout << "xlo " << xlo << " xhi " << xhi << " ylo " << ylo << " yhi " << yhi << " zlo " << zlo << " zhi " << zhi << endl;
        double lx = xhi - xlo;
        double ly = yhi - ylo;
        double lz = zhi - zlo;
	//set up for P2 array 
	int ndim = coordinateindices.size();
	assert (ndim==1 || ndim==2 || ndim ==3);
	vector<int> nbins(3,1);
	vector <double> xbounds = {xlo,xhi};
	vector <double> ybounds = {ylo,yhi};
	vector <double> zbounds = {zlo,zhi};
	for (int i=0; i <ndim; ++i){
		int c = coordinateindices[i];
		if (c==0){
			xbounds[0] = xlo;
			//cout <<  "xlo " << xlo << endl;
			xbounds[1] = xlo+resolutions[0];
			//cout << "xlo+resolutions[0] " << xlo+resolutions[0] << endl;
			double xval = xlo+resolutions[0];
			//cout << "xval " << xval << endl;
			while (xval+resolutions[0] < xhi+resolutions[0]){
				//cout << "xval " << xval << endl;
				xval = xval + resolutions[0];
			
				xbounds.push_back(xval);
			}
			//cout << "len(xbounds) " << xbounds.size() << endl;
			for (int ii=0; ii < xbounds.size(); ii++){
			//cout << "xbounds ii " << ii << " value " << xbounds[ii] << endl;
			}
		}
		else if (c==1){
			ybounds[0] = ylo;	
			//cout << "ylo " << ylo << endl;
			//cout << "yhi " << yhi << endl;
			ybounds[1] = ylo+resolutions[1];
			//cout << "ylo+resolutions[1] " << ylo+resolutions[1] << endl;
			double yval = ylo+resolutions[1];
			//cout << "yhi+resolutions[1] " << yhi+resolutions[1] << endl;
			while (yval+resolutions[1] < yhi+resolutions[1]){
				//cout << "yval " << yval << endl;
				yval = yval + resolutions[1];
				ybounds.push_back(yval);
			}
			//cout << "yval " << yval << endl;
			//cout << "len(ybounds) " << ybounds.size() << endl;
		
		}	
		else if (c==2){
			zbounds[0] = zlo;
			//cout << "zlo " << zlo << endl;
			zbounds[1] = zlo+resolutions[2];	
			//cout << "zlo+resolutions[0] " << zlo+resolutions[0] << endl;
			double zval = zlo+resolutions[2];
			while (zval+resolutions[2] < zhi+resolutions[2]){
				//cout << "zval " << zval << endl;
				zval = zval + resolutions[2];
				
				zbounds.push_back(zval);
			}
			//cout << "len(zbounds) " << zbounds.size() << endl;
		}
		else {
		cout << "Error in coordinateindicies input" << endl;
		}
	}



	//cout << "going into calcP2localanglescostheta function " << endl;
	vector<vector<vector<vector<double>>>> allcosthetaij =  calcP2localanglescostheta(xbounds, ybounds, zbounds,anglesinfo, atomsinfo, boxcoords, resolutions);
	//calculate P2 for each bin
	int nbinsx = xbounds.size() - 1;
	int nbinsy = ybounds.size() - 1;
	int nbinsz = zbounds.size() - 1;
	
	assert (nbinsx == allcosthetaij.size());
	assert (nbinsy == allcosthetaij[0].size());
	assert (nbinsz == allcosthetaij[0][0].size());
	//cout << "nbinsx: " << nbinsx << endl;
	//cout << "nbinsy: " << nbinsy << endl;
	//cout << "nbinsz: " << nbinsz << endl;
	
	vector<vector<vector<double>>> P2local;
	P2local.resize(nbinsx);
	for (int i=0; 	i < nbinsx; i++){
		P2local[i].resize(nbinsy);
		for (int j=0; j < nbinsy; j++){
			P2local[i][j].resize(nbinsz);
		}
	}
	vector <double> costhetas;
	for (int i=0; i < nbinsx; i++){
		for (int j=0; j < nbinsy; j++){
			for (int k=0; k < nbinsz; k++){
				costhetas = allcosthetaij[i][j][k];
				//cout << "len(costhetas) " << costhetas.size() << endl;		
				//for (int iii = 0; iii < costhetas.size(); iii++){
					//cout << costhetas[iii] << endl;
				//}
				if (costhetas.size() >=1){
					//cout << "in if costhetas.size() " << costhetas.size() << endl;
					vector <double> costhetasqrd;
					for (int ii=0; ii < costhetas.size(); ii++){
						
						costhetasqrd.push_back(pow(costhetas[ii],2));
					}
					double average = accumulate(costhetasqrd.begin(),costhetasqrd.end(),0.0)/costhetasqrd.size();
					//cout << "average " << average << endl;
					P2local[i][j][k] = (3.0*average -1.0)/2.0;
					//cout << "P2local[i][j][k] " << P2local[i][j][k] << " i " << i << " j " << j << " k " << k << endl;
					
				}
				else{
					P2local[i][j][k] = std::numeric_limits<double>::quiet_NaN();
					//cout << "no need to calculate P2local here" << endl;
				}
				
			}
		}
	}
	
	return P2local;
	
}



double calculatecrystallinity(vector<vector<vector<double>>> P2localxyz,double criteria = 0.6){
	int nbinsx = P2localxyz.size();
	int nbinsy = P2localxyz[0].size();
	int nbinsz = P2localxyz[0][0].size();
	vector<vector<vector<double>>> CrystallinityLocal;
	CrystallinityLocal.resize(nbinsx);
	for (int i = 0; i < nbinsx; i++){
		CrystallinityLocal[i].resize(nbinsy);
		for (int j = 0; j < nbinsy; j++){
			CrystallinityLocal[i][j].resize(nbinsz);
		}
	}
	int validpoints=0;
	double CrystallinityPoints=0.0;
	for (int i = 0; i < nbinsx; i++){
		for (int j= 0; j < nbinsy; j++){
			for (int k=0; k < nbinsz; k++){
				if (P2localxyz[i][j][k] > criteria){
					CrystallinityLocal[i][j][k] = 1.0;
					CrystallinityPoints = CrystallinityPoints + 1.0;
					validpoints++;
				}
				else if (::isnan(P2localxyz[i][j][k])==false){
					CrystallinityLocal[i][j][k] = 0.0;
					validpoints++;
				}
		
			}
		}
	}
	double crystallinity = CrystallinityPoints/validpoints;
	//cout << "crystallinity " << crystallinity << endl;
	return crystallinity;

}





LmpsTraj readlammpstraj(string lammpstraj,int nlines){
	
	bool looktimestep = false;
	bool looknumber = false;
	bool atomcoords = false;
	bool lookbox = false;
	int atomscount = 0;
	int nboxcoords = 0;
	int natoms;
	int xind  = 0;
	int yind = 0;
	int zind = 0;
	Atominfo ai;
	vector <vector <double>> boxcoords;
	boxcoords.resize(3);
	for (int i =0; i < 3; i++){
		boxcoords[i].resize(2);
	}
	vector <int> timesteps;
	vector <int> numbersatoms;
	vector <Atominfo> atomsinfoall;
	vector <vector < vector <double>>> boxcoordsall;
	
	ifstream infile;
	infile.open(lammpstraj.c_str());
	string line;
	char * pch;
	for (int i =0; i < nlines; i++){
	//cout << " i " << i << endl; 
	vector<string> vec_line;
	getline(infile,line);
	//cout << "line " << line << endl;
	if (line.length() > 0){

		char lineArr[line.length()];
		strcpy(lineArr, line.c_str());
		pch = strtok(lineArr, " " );
		while (pch != NULL) { vec_line.push_back(pch); pch = strtok(NULL, " " );}
		//cout << "vec_line" << vec_line[0] << endl;
		if (vec_line[0] == "ITEM:"){
			if (vec_line[1].compare("TIMESTEP") == 0){
				//cout << "TIMESTEP read" << endl;
				looktimestep = true;
				atomcoords = false;
			}
			else if (vec_line[1].compare("NUMBER")==0){
				//cout << "NUMBER read" << endl;
				looknumber = true;
			}
			else if (vec_line[1].compare("BOX")==0){
				lookbox = true;
				vector <vector <double>> boxcoords;
				nboxcoords = 0;
       				boxcoords.resize(3);
				for (int i =0; i < 3; i++){
                                boxcoords[i].resize(2);
                                }
				//cout << "Resized boxcoords going to resize boxcoordsall " << endl;
				boxcoordsall.resize(boxcoordsall.size()+1);
                		boxcoordsall[boxcoordsall.size()-1].resize(3);
				for (int i = 0; i < 3; i++){
              				  boxcoordsall[boxcoordsall.size()-1][i].resize(2);
        			}
        		}	
			else if (vec_line[1].compare("ATOMS")==0){
				xind = find(vec_line.begin(),vec_line.end(),"x") - vec_line.begin() -2;
				yind = find(vec_line.begin(),vec_line.end(),"y") - vec_line.begin() -2;
				zind = find(vec_line.begin(),vec_line.end(),"z") - vec_line.begin() -2;
				atomscount = 0;
				atomcoords = true;
				natoms = numbersatoms.back();
				ai.natoms = natoms;
				ai.atomnumber.resize(natoms);
				ai.atomtype .resize(natoms);
                        	ai.xcoords.resize(natoms);
                        	ai.ycoords.resize(natoms);
                        	ai.zcoords.resize(natoms);
			}
		}
		else if (looktimestep){
			timesteps.push_back(stoi(vec_line[0].c_str()));
			//cout << " timestep " << vec_line[0] << endl;
			looktimestep = false;
		}
		else if (looknumber){
			numbersatoms.push_back(stoi(vec_line[0].c_str()));
			looknumber = false;
		}
		else if (lookbox and nboxcoords < 6){
			boxcoords[floor(nboxcoords/2)][0] = atof(vec_line[0].c_str());
			//cout << " boxcoords element " << atof(vec_line[0].c_str()) << endl;
			nboxcoords++;
			boxcoords[floor(nboxcoords/2)][1] = atof(vec_line[1].c_str());
			//cout << " boxcoords element " << atof(vec_line[1].c_str()) << endl;
			nboxcoords++;
			if (nboxcoords >=6){
				lookbox = false;
				boxcoordsall[boxcoordsall.size()-1][0][0] = boxcoords[0][0];
				boxcoordsall[boxcoordsall.size()-1][0][1] = boxcoords[0][1];
				boxcoordsall[boxcoordsall.size()-1][1][0] = boxcoords[1][0];
				boxcoordsall[boxcoordsall.size()-1][1][1] = boxcoords[1][1];
				boxcoordsall[boxcoordsall.size()-1][2][0] = boxcoords[2][0];
				boxcoordsall[boxcoordsall.size()-1][2][1] = boxcoords[2][1];
			}
		}
		else if (atomcoords && atomscount < natoms){
			int atomind = stoi(vec_line[0].c_str()) -1;
			//cout << "atomind: " << atomind << endl;
			for (int jj=0; jj < vec_line.size(); jj++){
				if (jj == xind){
					ai.xcoords[atomind] = atof(vec_line[jj].c_str());
					//cout << "xcoord " << ai.xcoords[atomind] << endl;
				}
				else if (jj == yind){
					ai.ycoords[atomind] = atof(vec_line[jj].c_str());
					//cout << "ycoord " << ai.ycoords[atomind] << endl;
				}
				else if (jj == zind){
					ai.zcoords[atomind] = atof(vec_line[jj].c_str());
					//cout << "zcoord " <<  ai.zcoords[atomind] << endl;
				}
				else if (jj == 0){
					ai.atomnumber[atomind] = stoi(vec_line[jj].c_str());
					//cout << "atomnumber " << ai.atomnumber[atomind] << endl;
					
				}
				else if (jj == 1){
					ai.atomtype[atomind] = stoi(vec_line[jj].c_str());
					//cout << "atomtype " << ai.atomtype[atomind] << endl;
				} 
			}
			atomscount++;
			if (atomscount == natoms){
				atomsinfoall.push_back(ai);
			}
		}

	}
	}
	
	LmpsTraj lmpitrj;
	lmpitrj.timesteps = timesteps;
	lmpitrj.numbersatoms =	numbersatoms;
	lmpitrj.boxcoordsall = boxcoordsall;
	lmpitrj.atomsinfoall = atomsinfoall;	
		
	return lmpitrj;	
	
}
	 

double calcCrystallinityLammpsData(string lammpsdatafile, double resolutions [3], int nlinesdata){
	Lmpsabainfo lmpi = read_finalstructwithatombondangleinfo(lammpsdatafile,nlinesdata);
        cout << "read lammps data file " << lammpsdatafile << endl;
	double boxcoords[3][2];
	boxcoords[0][0] = lmpi.xlo;
	boxcoords[0][1] = lmpi.xhi;
	boxcoords[1][0] = lmpi.ylo;
	boxcoords[1][1] = lmpi.yhi;
	boxcoords[2][0] = lmpi.zlo;
	boxcoords[2][1] = lmpi.zhi;
	vector <int> coordinateindices = {0,1,2};
	vector<vector<vector<double>>> P2localxyzhere = localP2(boxcoords, lmpi.ai, lmpi.ani, coordinateindices, resolutions);
	cout << "calculated P2localxyzhere" << endl;

	double crystallinity = calculatecrystallinity(P2localxyzhere);
	return crystallinity;
}

vector <double> calcCrystallinityLammpstrajdata(string lammpsdatafile, string lammpstrajfile, double resolutions [3], int nlinesdata, int nlinestraj, int itframes=1){

	Lmpsabainfo lmpi = read_finalstructwithatombondangleinfo(lammpsdatafile,nlinesdata);	
	cout << "read lammps data file " << lammpsdatafile <<  endl;
	
	LmpsTraj lmpitrj = readlammpstraj(lammpstrajfile,nlinestraj);
	cout << "read lammps trajectory file " <<lammpstrajfile <<  endl;
	cout << "size timesteps " << lmpitrj.timesteps.size() << " size numbersatoms " << lmpitrj.numbersatoms.size() << endl;
	cout << "size boxcoordsall " << lmpitrj.boxcoordsall.size() << endl;
	assert (lmpitrj.timesteps.size() == lmpitrj.numbersatoms.size());
	int nframes = lmpitrj.timesteps.size();
	cout << "nframes " << nframes << endl;
	vector <double> crystallinities;
	vector <int> coordinateindicies = {0,1,2};
	for (int i = 0; i < nframes; i+=itframes){
		cout << "frame " << i << " timestep " << lmpitrj.timesteps[i] << endl;
		double boxcoords[3][2];
		boxcoords[0][0] = lmpitrj.boxcoordsall[i][0][0];
		boxcoords[0][1] = lmpitrj.boxcoordsall[i][0][1];
		boxcoords[1][0] = lmpitrj.boxcoordsall[i][1][0];
		boxcoords[1][1] = lmpitrj.boxcoordsall[i][1][1];
		boxcoords[2][0] = lmpitrj.boxcoordsall[i][2][0];
		boxcoords[2][1] = lmpitrj.boxcoordsall[i][2][1];
		vector<vector<vector<double>>> P2localxyzhere = localP2(boxcoords, lmpitrj.atomsinfoall[i], lmpi.ani, coordinateindicies, resolutions);
		cout << "calculated P2localxyzhere" << endl;
		double crystallinityhere =  calculatecrystallinity(P2localxyzhere);
		crystallinities.push_back(crystallinityhere);
	}
	return crystallinities;

}

vector <vector <int>> recursiveBackbone(vector <vector <int>> &allpaths, vector <int> &path, int &backboneBonds, vector <vector <int>> & graph, vector <bool> &visited){
	int u = path.back();
	
	for (int v: graph[u]){
		
		if (!visited[v] && v > u){
			

			path.push_back(v);
			visited[v] = true;

			if (path.size() < backboneBonds +1){
				recursiveBackbone(allpaths, path, backboneBonds, graph, visited);
			}
			else if (path.size() == backboneBonds +1){
				allpaths.push_back(path);
				
			}
		}
		
	}
	//cout << "" << endl;
	return allpaths;
}

vector<vector<vector<vector <int>>>>  sortVectorsBins(Lmpsabainfo lmpi, Angleinfo vectorsinfo, double resolutions [3], vector <int> coordinateindices){
	// sort the vectors defined by anglesinfo into bins based on coordinates
	int ndim = coordinateindices.size();
	vector <double> xbounds (1);
	vector <double> ybounds (1);
	vector <double> zbounds (1);
	double xlo = lmpi.xlo;
	double ylo = lmpi.ylo;
	double zlo = lmpi.zlo;
	double xhi = lmpi.xhi;
	double yhi = lmpi.yhi;
	double zhi = lmpi.zhi;
	xbounds[0] = lmpi.xlo;
	ybounds[0] = lmpi.ylo;
	zbounds[0] = lmpi.zlo;
	for (int i=0; i <ndim; ++i){
		int c = coordinateindices[i];
		if (c==0){
		
			//cout << "xlo+resolutions[0] " << xlo+resolutions[0] << endl;
			double xval = xlo+resolutions[0];
			//cout << "xval " << xval << endl;
			while (xval < xhi - resolutions[0]){
				//cout << "xval " << xval << endl;
				xbounds.push_back(xval);
				xval = xval + resolutions[0];
				
			}
			
		}
		else if (c==1){
			
			//cout << "ylo " << ylo << endl;
			//cout << "yhi " << yhi << endl;
			
			double yval = ylo+resolutions[1];
			//cout << "yhi+resolutions[1] " << yhi+resolutions[1] << endl;
			while (yval < yhi - resolutions[1] ){
				//cout << "yval " << yval << endl;
				ybounds.push_back(yval);
				yval = yval + resolutions[1];
				
			}
			//cout << "yval " << yval << endl;
			//cout << "len(ybounds) " << ybounds.size() << endl;
		
		}	
		else if (c==2){
		
			//cout << "zlo " << zlo << endl;
			double zval = zlo+resolutions[2];
			while (zval < zhi - resolutions[2]){
				//cout << "zval " << zval << endl;
				zbounds.push_back(zval);
				zval = zval + resolutions[2];
				
				
			}
			//cout << "len(zbounds) " << zbounds.size() << endl;
		}
		else {
		cout << "Error in coordinateindicies input" << endl;
		} 
	}

		int nbinsx = xbounds.size();
		int nbinsy = ybounds.size();
		int nbinsz = zbounds.size();

		vector<vector<vector<vector<int>>>> sameBinVectors;

		sameBinVectors.resize(nbinsx);
		for (int i=0; i <  nbinsx; i++){
			sameBinVectors[i].resize(nbinsy);
			for (int j = 0; j < nbinsy; j++){
				sameBinVectors[i][j].resize(nbinsz);

				for (int k =0; k < nbinsz; k++){
				sameBinVectors[i][j][k].resize(0);
				}
				
			}
		}


		//go through all vectors saved in vectorsinfo

		int vectors = vectorsinfo.nangles;
		int atomM;		
		
		int xbin;
		int ybin;
		int zbin;
	
		for (int i = 0; i < vectors; i++){
			atomM = vectorsinfo.anglesecond[i]-1;

			double atomcoords[3] = {lmpi.ai.xcoords[atomM],lmpi.ai.ycoords[atomM], lmpi.ai.zcoords[atomM]};
			//cout << " x: " << atomcoords[0] << " y: " << atomcoords[1] << " z: " << atomcoords[2] << endl;
                        vector <int> bin = squarewellconvbin(atomcoords,xbounds,ybounds,zbounds,resolutions);
			xbin = bin[0];
			ybin = bin[1];
			zbin = bin[2];
			sameBinVectors[xbin][ybin][zbin].push_back(i);
			//cout << "xbin " << xbin << " ybin " << ybin << " zbin " << zbin << " i " << i << endl;
		}

	return sameBinVectors; 

}

double calculateCrystHermansP2orientation(vector<vector<vector<double>>> P2localxyzhere, vector<vector<vector<double>>> cos2HOFZlocalhere,  double criteria){
//Calculates Hermans orientation factor only with respect to the crystalline domains
	vector <double> cos2HOFz; 
	int nbinsx = P2localxyzhere.size();
	int nbinsy = P2localxyzhere[0].size();
	int nbinsz = P2localxyzhere[0][0].size();

	for (int i = 0; i < nbinsx; i++){
		for (int j= 0; j < nbinsy; j++){
			for (int k=0; k < nbinsz; k++){
				if (P2localxyzhere[i][j][k] > criteria){
					cos2HOFz.push_back(cos2HOFZlocalhere[i][j][k]);	
					//cout << cos2HOFZlocalhere[i][j][k] << endl;
				}
		
			}
		}
	}

	double avgCosP2Z = accumulate(cos2HOFz.begin(), cos2HOFz.end(), 0.0)/cos2HOFz.size();

	double CrystP2Z = (3.0*avgCosP2Z -1.0)/2.0;

	return CrystP2Z;

}


InfoLmpsData calcP2CrystandlocalLammpsData(string lammpsdatafile, double resolutions [3], int nlinesdata, vector <int> coordinateindicies, int backboneBonds, vector <int> backboneAtypes, vector <double> zVector){
	Lmpsabainfo lmpi = read_finalstructwithatombondangleinfo(lammpsdatafile, nlinesdata);
	double resolutions10 [3] = {4,4,4};
	//cout << "read lammps data file " << lammpsdatafile << endl;
	Angleinfo anglesinfo;
	vector <vector <int>> graph (lmpi.ai.natoms);
	int atom1;
	int atom2;
	int atomMid;
	int atype1;
	int atype2;
	for (int i = 0; i < lmpi.bi.nbonds; i++){
		
		atom1 = lmpi.bi.bondfirst[i]-1;
		atype1 = lmpi.ai.atomtype[atom1];
		
		atom2 = lmpi.bi.bondsecond[i]-1;
		atype2 = lmpi.ai.atomtype[atom2];
		if (find(backboneAtypes.begin(), backboneAtypes.end(), atype1) != backboneAtypes.end() && find(backboneAtypes.begin(), backboneAtypes.end(), atype2) != backboneAtypes.end()){
			graph[atom1].push_back(atom2);
			graph[atom2].push_back(atom1);
		
		}
	}
	//cout << "made graph " << endl;
	anglesinfo.nangles = 0;
	for (int i = 0; i < lmpi.ai.natoms; i++){
		vector <vector <int>> allpaths;
		vector <int> start = {i};
		vector <bool> visited (lmpi.ai.natoms, false);
		//if (i == 0){
		//cout << "start " << i << endl;
		//}
		allpaths = recursiveBackbone(allpaths, start, backboneBonds, graph, visited);
		for (int j = 0; j < allpaths.size(); j++){
			anglesinfo.nangles++;
			atom1 = allpaths[j][0]+1;
			atom2 = allpaths[j][backboneBonds]+1;
			atomMid = allpaths[j][backboneBonds/2]+1;
			if (find(anglesinfo.anglefirst.begin(), anglesinfo.anglefirst.end(), atom1) == anglesinfo.anglefirst.end() && find(anglesinfo.anglethird.begin(), anglesinfo.anglethird.end(), atom2) == anglesinfo.anglethird.end()){
				anglesinfo.anglefirst.push_back(atom1);
				anglesinfo.anglesecond.push_back(atomMid);
				anglesinfo.anglethird.push_back(atom2);
				//if (i == 0){
				//	cout << " atom1 " << anglesinfo.anglefirst.back() << " atom2 " << anglesinfo.anglethird.back() ;
				//}
			}
			

		}
	}
	//cout << "atom1[0] " << anglesinfo.anglefirst[0] << " atom2[0] " << anglesinfo.anglethird[0]  << endl;



	double boxcoords[3][2];
	boxcoords[0][0] = lmpi.xlo;
	boxcoords[0][1] = lmpi.xhi;
	boxcoords[1][0] = lmpi.ylo;
	boxcoords[1][1] = lmpi.yhi;
	boxcoords[2][0] = lmpi.zlo;
	boxcoords[2][1] = lmpi.zhi;
	//vector<vector<vector<double>>> P2localhere = localP2(boxcoords, lmpi.ai, lmpi.ani, coordinateindicies, resolutions);
	//cout << "localP2 specified: ";
        //        for (int ii=0; ii <P2localhere.size(); ii++){
        //                for (int jj=0; jj < P2localhere[ii].size(); jj++){
        //                        for (int kk=0; kk < P2localhere[ii][jj].size(); kk++){
        //                                cout << P2localhere[ii][jj][kk] << " ii " << ii << " jj " << jj << " kk " << kk << endl;
        //                       }
        //                }
        //        }
		//vector <int> xyzcoordinates = {0,1,2};
	vector <int> xyzcoordinates = {0,1,2};
	vector <int> profileZcoordinates = {2};
	//cout << "going to sortVectorsBins " << endl;
	vector<vector<vector<vector<int>>>> sortedVectorsBins =  sortVectorsBins(lmpi, anglesinfo, resolutions, xyzcoordinates);
	vector<vector<vector<vector<int>>>> sortedVectorsBins10 = sortVectorsBins(lmpi, anglesinfo, resolutions10, xyzcoordinates);
	
	//cout << "Done sorting" << endl;
	vector<vector<vector<double>>> P2localxyzhere = localP2PreSorted(lmpi,  anglesinfo, sortedVectorsBins);
	double crystallinityhere = calculatecrystallinity(P2localxyzhere);

	vector<vector<vector<double>>> P2localxyzhere10 = localP2PreSorted(lmpi,  anglesinfo, sortedVectorsBins10);

	vector<vector<vector<double>>> cos2HOFZlocalhere10 =  localcos2HOFZPreSorted(lmpi, anglesinfo, sortedVectorsBins10, zVector);

	//cout << "crystallinity: " << crystallinityhere << endl; 

	//vector<vector<vector<vector<int>>>> sortedVProfile = sortVectorsBins(lmpi, anglesinfo, resolutions, profileZcoordinates);
	//vector<vector<vector<double>>> P2localProfile = localP2PreSorted(lmpi, anglesinfo, sortedVProfile);
	//for (int ii=0; ii < P2localProfile.size(); ii++){
	//	for (int jj=0; jj < P2localProfile[ii].size(); jj++){
	//		for (int kk=0; kk < P2localProfile[ii][jj].size(); kk++){
	//			cout << P2localProfile[ii][jj][kk] << " ii " << ii << " jj " << jj << " kk " << kk << endl;
	//		}
	//	}
	//}




	//cout << "going to calculateP2" << endl;
	//double P2here = calculateP2(boxcoords,lmpi.ai, anglesinfo);
	//cout << "P2 orientation: " << P2here << endl;
	//double hermansOrientP2 = calculateHermansP2orientation(boxcoords, lmpi.ai, anglesinfo);
	double hermansCOrientP2 = calculateCrystHermansP2orientation(P2localxyzhere10, cos2HOFZlocalhere10, 0.6);
	//cout << "Herman's Orientation factor: " << hermansOrientP2 << endl;
	cout << lammpsdatafile << ", " << crystallinityhere << ", " << hermansCOrientP2 << endl;



	InfoLmpsData info;
	
	info.crystallinities = crystallinityhere;
	//info.orientationP2 = P2here;
	//info.localP2 = P2localProfile;
	info.HermansOrientP2 = hermansCOrientP2;
	return info;
		
	
	

}

vector<vector<vector<vector<double>>>> calclocalDensityLammpstrajdata(string lammpsdatafile, string lammpstrajfile, double resolutions [3], int nlinesdata, int nlinestraj, vector <int> coordinateindices, string outfile, int itframes=1){
	Lmpsabainfo lmpi = read_finalstructwithatombondangleinfo(lammpsdatafile,nlinesdata);
	cout << "read lammps data file " << lammpsdatafile << endl;
	LmpsTraj lmpitrj = readlammpstraj(lammpstrajfile, nlinestraj);
	cout << "read lammps trajectory file " << lammpstrajfile << endl;
	assert (lmpitrj.timesteps.size() == lmpitrj.numbersatoms.size());
	int nframes = lmpitrj.timesteps.size();
	cout << "localDensity nframes: " << nframes << endl;
	vector<vector<vector<vector<double>>>> localDensityall;
	#pragma omp parallel for
	for (int i=0; i < nframes; i+=itframes){
		ofstream fout;
		fout.open("frame"+to_string(i)+outfile, ofstream::app);
		
		int omp_get_thread_num();
                int omp_get_num_threads();
		int omp_get_num_procs();
		int omp_get_max_threads();
		printf("Max Threads: %d\n", omp_get_max_threads());
		printf("Num of CPU: %d\n", omp_get_num_procs());
                printf("NumberofThreads:%d\n", omp_get_num_threads());
                printf("Thread rank: %d\n", omp_get_thread_num());
		fout << "frame " << i << " timestep " << lmpitrj.timesteps[i] << endl;
		double boxcoords[3][2];
                boxcoords[0][0] = lmpitrj.boxcoordsall[i][0][0];
                boxcoords[0][1] = lmpitrj.boxcoordsall[i][0][1];
                boxcoords[1][0] = lmpitrj.boxcoordsall[i][1][0];
                boxcoords[1][1] = lmpitrj.boxcoordsall[i][1][1];
                boxcoords[2][0] = lmpitrj.boxcoordsall[i][2][0];
                boxcoords[2][1] = lmpitrj.boxcoordsall[i][2][1];
                cout << "xlo " << boxcoords[0][0] << " xhi " << boxcoords[0][1] << " ylo " << boxcoords[1][0] << " yhi " << boxcoords[1][1] << " zlo " << boxcoords[2][0] << " zhi " << boxcoords[2][1] << endl;
		int ndim = coordinateindices.size();
                double xlo = boxcoords[0][0];
                double xhi = boxcoords[0][1];
                double ylo = boxcoords[1][0];
                double yhi = boxcoords[1][1];
                double zlo = boxcoords[2][0];
                double zhi = boxcoords[2][1];
        assert (ndim==1 || ndim==2 || ndim ==3);
        vector<int> nbins(3,1);
        vector <double> xbounds = {xlo,xhi};
        vector <double> ybounds = {ylo,yhi};
        vector <double> zbounds = {zlo,zhi};
        for (int i=0; i <ndim; ++i){
                int c = coordinateindices[i];
                if (c==0){
                        xbounds[0] = xlo;
                        //cout <<  "xlo " << xlo << endl;
                        xbounds[1] = xlo+resolutions[0];
                        //cout << "xlo+resolutions[0] " << xlo+resolutions[0] << endl;
                        double xval = xlo+resolutions[0];
                        //cout << "xval " << xval << endl;
                        while (xval+resolutions[0] < xhi+resolutions[0]){
                                //cout << "xval " << xval << endl;
                                xval = xval + resolutions[0];

                                xbounds.push_back(xval);
                        }
                        //cout << "len(xbounds) " << xbounds.size() << endl;
                        for (int ii=0; ii < xbounds.size(); ii++){
                        //cout << "xbounds ii " << ii << " value " << xbounds[ii] << endl;
                        }
                }
		else if (c==1){
                        ybounds[0] = ylo;
                        //cout << "ylo " << ylo << endl;
                        //cout << "yhi " << yhi << endl;
                        ybounds[1] = ylo+resolutions[1];
                        //cout << "ylo+resolutions[1] " << ylo+resolutions[1] << endl;
                        double yval = ylo+resolutions[1];
                        //cout << "yhi+resolutions[1] " << yhi+resolutions[1] << endl;
                        while (yval+resolutions[1] < yhi+resolutions[1]){
                                //cout << "yval " << yval << endl;
                                yval = yval + resolutions[1];
                                ybounds.push_back(yval);
                        }
                        //cout << "yval " << yval << endl;
                        //cout << "len(ybounds) " << ybounds.size() << endl;

                }
                else if (c==2){
                        zbounds[0] = zlo;
                        //cout << "zlo " << zlo << endl;
                        zbounds[1] = zlo+resolutions[2];
                        //cout << "zlo+resolutions[0] " << zlo+resolutions[0] << endl;
                        double zval = zlo+resolutions[2];
                        while (zval+resolutions[2] < zhi+resolutions[2]){
                                //cout << "zval " << zval << endl;
                                zval = zval + resolutions[2];

                                zbounds.push_back(zval);
                        }
                }


	}
	int nbinsx = xbounds.size() -1;
                int nbinsy = ybounds.size() -1;
                int nbinsz = zbounds.size() -1;
                fout << "nbinsx: " << nbinsx << endl;
                fout << "nbinsy: " << nbinsy << endl;
                fout << "nbinsz: " << nbinsz << endl;
		fout << "xlo " << xlo << " xhi " << xhi << " ylo " << ylo << " yhi " << yhi << " zlo " << zlo << " zhi " << zhi << endl;	

	vector<vector<vector<double>>> localDensityhere = localDensity(boxcoords, lmpitrj.atomsinfoall[i], lmpi, coordinateindices, resolutions);
                localDensityall.push_back(localDensityhere);
                fout << "localDensity specified: ";
                for (int ii=0; ii <localDensityhere.size(); ii++){
                        for (int jj=0; jj < localDensityhere[ii].size(); jj++){
                                for (int kk=0; kk < localDensityhere[ii][jj].size(); kk++){
                                        fout << localDensityhere[ii][jj][kk] << " ii " << ii << " jj " << jj << " kk " << kk << endl;
                                }
                        }
                }

	}

	return localDensityall;

}

InfoLmpsTraj calcP2CrystandlocalLammpstrajdata(string lammpsdatafile, string lammpstrajfile, double resolutions [3], int nlinesdata, int nlinestraj, vector <int> coordinateindices,string outfile, int itframes=1){
	Lmpsabainfo lmpi = read_finalstructwithatombondangleinfo(lammpsdatafile,nlinesdata);
	cout << "read lammps data file " << lammpsdatafile << endl;
	LmpsTraj lmpitrj = readlammpstraj(lammpstrajfile, nlinestraj);
	cout << "read lammps trajectory file " << lammpstrajfile << endl;
	cout << "size timesteps " << lmpitrj.timesteps.size() << " size numbersatoms " << lmpitrj.numbersatoms.size() << endl;
	cout << "size boxcoordsall " << lmpitrj.boxcoordsall.size() << endl;
	assert (lmpitrj.timesteps.size() == lmpitrj.numbersatoms.size());
	int nframes = lmpitrj.timesteps.size();
	cout << "OrderParam nframes: " << nframes << endl;
	vector <double> crystallinities;
	vector <double> orientationP2;
	vector <double> HermansorientationP2;
	vector<vector<vector<vector<double>>>> P2localall;
	vector <int> xyzcoordinates = {0,1,2};
	InfoLmpsTraj info;
	#pragma omp parallel for
	for (int i =0; i < nframes; i+=itframes){
		ofstream fout;
                fout.open("frame"+to_string(i)+outfile, ofstream::app);
		int omp_get_thread_num();
                int omp_get_num_threads();
		int omp_get_num_procs();
		int omp_get_max_threads();
                printf("Max Threads: %d\n", omp_get_max_threads());
                printf("Num of CPU: %d\n", omp_get_num_procs());
                printf("NumberofThreads:%d\n", omp_get_num_threads());
                printf("Thread rank: %d\n", omp_get_thread_num());
		fout << "frame " << i << " timestep " << lmpitrj.timesteps[i] << endl;
		double boxcoords[3][2];
                boxcoords[0][0] = lmpitrj.boxcoordsall[i][0][0];
                boxcoords[0][1] = lmpitrj.boxcoordsall[i][0][1];
                boxcoords[1][0] = lmpitrj.boxcoordsall[i][1][0];
                boxcoords[1][1] = lmpitrj.boxcoordsall[i][1][1];
                boxcoords[2][0] = lmpitrj.boxcoordsall[i][2][0];
                boxcoords[2][1] = lmpitrj.boxcoordsall[i][2][1];
		fout << "xlo " << boxcoords[0][0] << " xhi " << boxcoords[0][1] << " ylo " << boxcoords[1][0] << " yhi " << boxcoords[1][1] << " zlo " << boxcoords[2][0] << " zhi " << boxcoords[2][1] << endl;
		int ndim = coordinateindices.size();
		double xlo = boxcoords[0][0];
		double xhi = boxcoords[0][1];
		double ylo = boxcoords[1][0];
		double yhi = boxcoords[1][1];
		double zlo = boxcoords[2][0];
		double zhi = boxcoords[2][1];
        assert (ndim==1 || ndim==2 || ndim ==3);
        vector<int> nbins(3,1);
        vector <double> xbounds = {xlo,xhi};
        vector <double> ybounds = {ylo,yhi};
        vector <double> zbounds = {zlo,zhi};
        for (int i=0; i <ndim; ++i){
                int c = coordinateindices[i];
                if (c==0){
                        xbounds[0] = xlo;
                        //cout <<  "xlo " << xlo << endl;
                        xbounds[1] = xlo+resolutions[0];
                        //cout << "xlo+resolutions[0] " << xlo+resolutions[0] << endl;
                        double xval = xlo+resolutions[0];
                        //cout << "xval " << xval << endl;
                        while (xval+resolutions[0] < xhi+resolutions[0]){
                                //cout << "xval " << xval << endl;
                                xval = xval + resolutions[0];

                                xbounds.push_back(xval);
                        }
                        //cout << "len(xbounds) " << xbounds.size() << endl;
                        for (int ii=0; ii < xbounds.size(); ii++){
                        //cout << "xbounds ii " << ii << " value " << xbounds[ii] << endl;
                        }
                }
                else if (c==1){
                        ybounds[0] = ylo;
                        //cout << "ylo " << ylo << endl;
                        //cout << "yhi " << yhi << endl;
                        ybounds[1] = ylo+resolutions[1];
                        //cout << "ylo+resolutions[1] " << ylo+resolutions[1] << endl;
                        double yval = ylo+resolutions[1];
                        //cout << "yhi+resolutions[1] " << yhi+resolutions[1] << endl;
                        while (yval+resolutions[1] < yhi+resolutions[1]){
                                //cout << "yval " << yval << endl;
                                yval = yval + resolutions[1];
                                ybounds.push_back(yval);
                        }
                        //cout << "yval " << yval << endl;
                        //cout << "len(ybounds) " << ybounds.size() << endl;

                }
                else if (c==2){
                        zbounds[0] = zlo;
                        //cout << "zlo " << zlo << endl;
                        zbounds[1] = zlo+resolutions[2];
                        //cout << "zlo+resolutions[0] " << zlo+resolutions[0] << endl;
                        double zval = zlo+resolutions[2];
                        while (zval+resolutions[2] < zhi+resolutions[2]){
                                //cout << "zval " << zval << endl;
                                zval = zval + resolutions[2];

                                zbounds.push_back(zval);
                        }
		}
	}
                        //cout << "len(zbounds) " << zbounds.size() 
		int nbinsx = xbounds.size() -1;
        	int nbinsy = ybounds.size() -1;
        	int nbinsz = zbounds.size() -1;	
		fout << "nbinsx: " << nbinsx << endl;
		fout << "nbinsy: " << nbinsy << endl;
		fout << "nbinsz: " << nbinsz << endl;
		fout << "xlo " << xlo << " xhi " << xhi << " ylo " << ylo << " yhi " << yhi << " zlo " << zlo << " zhi " << zhi << endl;
		vector<vector<vector<double>>> P2localhere = localP2(boxcoords, lmpitrj.atomsinfoall[i], lmpi.ani, coordinateindices, resolutions);
		P2localall.push_back(P2localhere);
		fout << "localP2 specified: ";
		for (int ii=0; ii <P2localhere.size(); ii++){
			for (int jj=0; jj < P2localhere[ii].size(); jj++){
				for (int kk=0; kk < P2localhere[ii][jj].size(); kk++){
					fout << P2localhere[ii][jj][kk] << " ii " << ii << " jj " << jj << " kk " << kk << endl;
				}
			}
		}
		
		vector<vector<vector<double>>> P2localxyzhere = localP2(boxcoords, lmpitrj.atomsinfoall[i], lmpi.ani, xyzcoordinates, resolutions);
		double crystallinityhere = calculatecrystallinity(P2localxyzhere);
		fout << "crystallinity " << crystallinityhere << endl;
		crystallinities.push_back(crystallinityhere);
		double P2here = calculateP2( boxcoords,lmpitrj.atomsinfoall[i], lmpi.ani);
		fout << "orientation P2 " << P2here << endl;
		orientationP2.push_back(P2here);
		double HermansP2here = calculateHermansP2orientation(boxcoords, lmpitrj.atomsinfoall[i], lmpi.ani);
		fout << "Hermansorientation P2 " << HermansP2here << endl;
		HermansorientationP2.push_back(HermansP2here);
	}
	info.crystallinities = crystallinities;
	info.orientationP2 = orientationP2;
	info.HermansOrientP2 = HermansorientationP2;
	info.localP2all = P2localall;
	
	return info;

}

int main(){
	time_t start,end;
	time (&start);
	double resolutions[3] = {4,4,4};
	vector <int> coordinateindices  = {2};
	string outfiledens = "PPCrystStartdens.txt";
	string outfileorder = "PPCrystStartorder.txt";
	
	//vector<vector<vector<vector<double>>>> localDensityall = calclocalDensityLammpstrajdata("PEmeltPolyHDPE10_66_smallerinboxcen.data", "meltcrystalizing_10_66_Poly_smaller.lammpstrj", resolutions, 91670, 688410, coordinateindices,outfiledens,4);
	//InfoLmpsData info = calcP2CrystandlocalLammpsData("../PEsemicrystalline/longermatch/PEcryst201EMClongrmrinfdabready.data", resolutions,  11078, coordinateindicies);
	//InfoLmpsTraj info = calcP2CrystandlocalLammpstrajdata("PEmeltPolyHDPE10_66_smallerinboxcen.data","meltcrystalizing_10_66_Poly_smaller.lammpstrj", resolutions, 91670, 688410, coordinateindices,outfileorder,4);
	//vector <double> crystallinities = calcCrystallinityLammpstrajdata("PEcryst201odbadready.data", "../PEsemicrystalline/fromhen2lessdensemiddle/chainsrem/PEsetrmrsim.lammpstraj", resolutions, 11078, 834372,5);
	 //double crystallinity = calcCrystallinityLammpsData("p613CrystStart.data", resolutions, 600676);
	int backboneBonds = 6;
	vector <int> backboneAtypes {2,3};
	cout << "Datafile, Crystallinity, Hermans Crystalline Orientation Factor" << endl;
	vector<double> zVector = {-0.16504760677,0,0.986285601537};
	InfoLmpsData info = calcP2CrystandlocalLammpsData("p613CrystStart.data", resolutions,641792, coordinateindices, backboneBonds, backboneAtypes, zVector);
	
	//cout << " P2 " << info.orientationP2 << endl;
	//cout << "crystallinity " << info.crystallinities << endl;
	
	//double crystallinity = calcCrystallinityLammpsData("PEcryst201odbadready.data", resolutions, 11078);
	
	
	
	
	time (&end);
	double dif = difftime(end,start);
	printf("Elapsed time is %.1f seconds.\n",dif);
	//cout << "going to calculate crystallinites from three frames" << endl;
	//vector <double> crystallinities = calcCrystallinityLammpstrajdata("PEcryst201odbadready.data", "../PEsemicrystalline/fromhen2lessdensemiddle/chainsrem/PEsetrmrsim.lammpstraj", resolutions, 11078, 834372,100);

	return 0;
}
