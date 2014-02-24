/*
 * 
 * 
    Copyright (C) 2012  Emad Noorizadeh   emad.noorizadeh@gmail.com

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "md.h"

/**=====================================================
 *            initializing variables   
 */

void mdpack::md::initVars(void){
  
  atoms= NULL;
  bonds= NULL;
  angles=NULL;
  dihedrals=NULL;
  impDihedrals= NULL;
  nonBondedPairs = NULL;
//   r= NULL;
//   v= NULL;
  f= NULL;
  strcpy(exclude, "1-4");
  cutoff = 12.0;
  switchdist = 11.0;
  
  // step size for integration
  h=TIMEFACTOR*0.001;
  // time step of simulation
  timeStep = h/TIMEFACTOR;
  // the nth step
  stepn = 0;
  initialTime =0.0;
  simTime= initialTime;
	
  //energies
  Ebonds=0.0;
  Eangles=0.0;
  Edihedrals =0.0;
  EimpDihedrals =0.0;
  EnonBonded =0.0;
  Ekinetic =0.0;
  Epotential = 0.0;
  Etotal =0.0;
  EsphBC =0.0;
  Evdw =0.0;
  Eelec =0.0;
  Emisc =0.0;
	
  //spherical boundary conditions
  sphBC = true;
  sphBCr1=8.0;
  sphBCk1 =0.01;
  sphBCexp1 =2.0;
    
  // variables for electorstatic forces
  e0 =1.0;
  inv_e0 =1.0/e0;
  e14 =1.0;
  
  // constant variables for non-bonded force computations
  inv_Rc2 = 1.0/(cutoff*cutoff);
  fct = 1.0/pow3(cutoff*cutoff-switchdist*switchdist);
  
  // temperatures
  tempInstant =0.0;
  tempAvg =0.0;
  targetTemp =0.0;
  
  // user defined names for simulation and output files
  strcpy(outputName, "run");
  strcpy(simName, "no simulation name");
  
  // conditional variables for writing output files
  isWritePosOn = false;
  isWriteVelOn = false;
  isWritePhiOn = false;
  isBinaryOn   = false;
  isNonBondedPairs = false;
  isCells = false;
  
  // conditional variable for using cell lists (neighbour lists)
  useCellLists = false;
  
  // setting up random number generator
  seed = (int)time(0); 
  W.setSeed(seed);
  
  // Langevin parameters
  gamma = 0.01;  //10 ps
  langFactor = NULL;
  langFactorh = NULL;
  // Nose-Hoover-Langevin parameters
  nhl_gamma = 0.01; // 10 ps
  nhl_mu= 1.0;
  nhl_imu =1.0/nhl_mu;
  nhl_xi=0.0;
  
}
/**============================================================*/

/**=====================================================
 *            setting up arrays and variables   
 */
void mdpack::md::setup(){
  nDegsFreedom = nAtoms*3;
  
  // constant factors for Langevin 
  expGamma_h = exp(-gamma*h);
  expGamma_hh = exp(-0.5*gamma*h);
  langFactor = new mdpack::dreal [nAtoms];
  langFactorh = new mdpack::dreal [nAtoms];
  for(int i=0; i<nAtoms; i++){
    langFactor[i] = sqrt(targetTemp*BOLTZMANN*atoms[i].im*(1.0-exp(-2.0*gamma*h)));
    langFactorh[i] = sqrt(targetTemp*BOLTZMANN*atoms[i].im*(1.0-exp(-gamma*h)));
  }
  // Constants factors for Nose-Hoover-Langevin
  expNhl= exp(-nhl_gamma*h);
  nhlFactor= sqrt((targetTemp*BOLTZMANN*(1.0-exp(-2.0*nhl_gamma*h)))/nhl_mu);
  // constant factor for over-damped dynamics
  overDampedFactor = sqrt(2.0*h*targetTemp*BOLTZMANN);
  
//   r = new mdpack::vector [nAtoms];
//   v = new mdpack::vector [nAtoms];
  f = new mdpack::vector3 [nAtoms];
  this->setLJtable();
	
  centerOfMass_m =0.0;
  max_x = -9999999999.0;
  max_y = -9999999999.0;
  max_z = -9999999999.0;
  min_x =  9999999999.0;
  min_y =  9999999999.0;
  min_z =  9999999999.0;
  for(int i=0; i<nAtoms; i++){
    centerOfMass_m += atoms[i].m;
    centerOfMassCoords += atoms[i].r*atoms[i].m;
    
    if(atoms[i].r.x > max_x)
      max_x = atoms[i].r.x;
    else if(atoms[i].r.x < min_x)
      min_x = atoms[i].r.x;
	  
    if( (atoms[i].r.y) > max_y)
      max_y = atoms[i].r.y;
    else if( (atoms[i].r.y) < min_y)
      min_y = atoms[i].r.y;
	  
    if( (atoms[i].r.z) > max_z)
      max_z = (atoms[i].r.z);
    else if(atoms[i].r.z < min_z)
      min_z = (atoms[i].r.z);
  }
  
  centerOfMassCoords *=(1.0/centerOfMass_m);
  sphBCCenter = centerOfMassCoords;
  if(sphBC)
    boxSide =2* sphBCr1;
  else{
    mdpack::dreal mx, my, mz;
    mx= max(max_x - centerOfMassCoords.x, centerOfMassCoords.x - min_x);
    my= max(max_y - centerOfMassCoords.y, centerOfMassCoords.y - min_y);
    mz= max(max_z - centerOfMassCoords.z, centerOfMassCoords.z - min_z);
    // std::cout<<" mx my mz "<<mx <<" "<<my<<" "<<mz<<std::endl;
    boxSide = 2.0* (max(max(mx, my),mz) + HALFCUTOFF);
    //boxSide = 2.0* max(max(max_x- min_x, max_y- min_y), max_z- min_z);
  }
	
}
/**==========================================================*/

/**=====================================================
 *            setting up cell lists   
 * 
 * this functions is not finished and to be completed and tested later
 */
void mdpack::md::setupCellLists(void){
  double divd = boxSide/cutoff;
  std::cout<<" divd "<< divd<<" ";
  double up = ceil(divd);
  std::cout<<" up "<< up<<" ";
  boxSide += (up-divd)*cutoff;
  std::cout<<" boxSide "<< boxSide<<" ";
  int nEachSide = (int) ((boxSide/cutoff) +0.1);
  std::cout<<" nEachSide "<< nEachSide<<" ";
  nCells = nEachSide*nEachSide*nEachSide;
	 
  cells = new mdpack::cell[nCells];
  isCells = true;
  double leftx  = centerOfMassCoords.x - 0.5*boxSide;
  // double rightx = centerOfMassCoords.x + 0.5*boxSide;
  double lefty  = centerOfMassCoords.y - 0.5*boxSide;
  // double righty = centerOfMassCoords.y + 0.5*boxSide;
  double leftz  = centerOfMassCoords.z - 0.5*boxSide;
  // double rightz = centerOfMassCoords.z + 0.5*boxSide;
	
	
  double x,y,z;
  for(int k=0; k<nEachSide; k++)
    for(int j=0; j<nEachSide; j++)
      for(int i=0; i<nEachSide; i++){
	x= leftx+(i+1.0)*cutoff -0.5*cutoff;
	y= lefty+(j+1.0)*cutoff -0.5*cutoff;
	z= leftz+(k+1.0)*cutoff -0.5*cutoff;
	cells[i + j*nEachSide +k*(nEachSide*nEachSide)]= mdpack::cell(x,y,z, cutoff);
      }
       
  
   
  double rc2 =(3.0*cutoff*cutoff);
  std::list<int> neighbour;
  for(int i=0; i<nCells; i++){
    neighbour.clear(); 
    for(int j=0; j<nCells; j++){
      if(j!=i && ( pow2(cells[i].x-cells[j].x) + pow2(cells[i].y-cells[j].y) 
	    + pow2(cells[i].z-cells[j].z)) <= rc2)
	neighbour.push_back(j);
    }
    cells[i].setNeighbours(neighbour);
  }
	
  for (int i=0; i<nAtoms; i++)
    for(int j=0; j<nCells; j++){
      if(cells[j].isIn(atoms[i].r))
	cells[j].addAtom(i+1);
    }
    
}
/**=======================================================*/


void mdpack::md::readinpsf(const char * psf){
  
  bool hitAtomSection= false;
  char part1 [50];
  char part2 [50];
  char sname [10];
  char rname [10];
  char an    [10];
  char at    [10];
  int aid, rid, idx1, idx2, idx3, idx4, unused;
  mdpack::dreal mass;
  mdpack::dreal charge;
  
  std::ifstream input;
  input.open(psf, std::ios::in);
  if (!input)
    terminate(" cant open psf file ");
  
  input >> part1;
  if (strcmp(part1,"PSF")==0)
      input.ignore(250,'\n');
  while(!input.eof() && !hitAtomSection){
    input >> part1 >> part2;
    if(strcmp(part2,"!NTITLE")==0)
      input.ignore(250,'\n');
    else if(strcmp(part1,"REMARKS")==0)
      input.ignore(250,'\n');
    else if(strcmp(part2,"!NATOM")==0){
      nAtoms = atoi(part1);
      hitAtomSection =true;
    }
    else
      terminate(" problem with reading psf file, cannot recognize the format  ");
  }
  // setting up array of atoms
 
  atoms = new mdpack::atom [nAtoms];
  for(int i=0; i<nAtoms; i++){
    input >> aid >> sname >> rid >> rname >> an >> at >> charge >> mass >> unused;
    atoms[i].set(aid,sname,rid,rname,an,at,charge, mass);
  }
  
  input >> part1 >> part2;
  input.ignore(250,'\n');
  if(strcmp(part2,"!NBOND:")==0)
    nBonds = atoi(part1);
  else
    terminate(" can't read !NBOND: in the psf file ");
 // std::cout<<"nBonds "<<nBonds<<std::endl;
  bonds = new mdpack::bond [nBonds];
  for(int i=0; i< nBonds; i++){
    input >> idx1 >>idx2;
    bonds[i].set(idx1, idx2, atoms[idx1-1].type, atoms[idx2-1].type);
    // adding atoms to their bond list, which is used later 
    // for excluding nonbonded interactions 
    atoms[idx1-1].list1_2.push_back(idx2);
    atoms[idx2-1].list1_2.push_back(idx1);
  }
  
  input >> part1 >> part2;
  input.ignore(250,'\n');
  if(strcmp(part2,"!NTHETA:")==0)
    nAngles = atoi(part1);
  else
    terminate(" can't read !NTHETA: in the psf file ");
   //std::cout<<"nAngles "<<nAngles<<std::endl;
  angles = new mdpack::angle [nAngles];
  for(int i=0; i< nAngles; i++){
    input >> idx1 >>idx2 >> idx3;
    angles[i].set(idx1, idx2,idx3, atoms[idx1-1].type, atoms[idx2-1].type, atoms[idx3-1].type);
    // adding atoms to their 1-3 list (angle list), which is used later 
    // for excluding nonbonded interactions 
    atoms[idx1-1].list1_3.push_back(idx3);
    atoms[idx3-1].list1_3.push_back(idx1);
  }
  
  input >> part1 >> part2;
  input.ignore(250,'\n');
  if(strcmp(part2,"!NPHI:")==0)
    nDihedrals = atoi(part1);
  else
    terminate(" can't read !NPHI: in the psf file ");
  
  //std::cout<<"nDihedrals "<<nDihedrals<<std::endl;
  dihedrals = new mdpack::dihedral [nDihedrals];
  for (int i=0; i< nDihedrals ; i++){
    input >> idx1 >>idx2 >> idx3 >> idx4;
    dihedrals[i].set(idx1, idx2,idx3, idx4, atoms[idx1-1].type, 
		      atoms[idx2-1].type, atoms[idx3-1].type, atoms[idx4-1].type);
    // adding atoms to their 1-4 list (dihedral list), which is used later 
    // for excluding nonbonded interactions 
    atoms[idx1-1].list1_4.push_back(idx4);
    atoms[idx4-1].list1_4.push_back(idx1);
    
  }
  
  input >> part1 >> part2;
  input.ignore(250,'\n');
  if(strcmp(part2,"!NIMPHI:")==0)
    nImpDihedrals = atoi(part1);
  else
    terminate(" can't read !NPHI: in the psf file ");
  
  //std::cout<<"nImpDihedrals "<<nImpDihedrals<<std::endl;
  impDihedrals = new mdpack::impDihedral [nImpDihedrals];
  for (int i=0; i< nImpDihedrals ; i++){
    input >> idx1 >>idx2 >> idx3 >> idx4;
    impDihedrals[i].set(idx1, idx2,idx3, idx4, atoms[idx1-1].type, 
				   atoms[idx2-1].type, atoms[idx3-1].type,atoms[idx4-1].type);
    
    atoms[idx1-1].list1_4.push_back(idx4);
    atoms[idx4-1].list1_4.push_back(idx1);
  }
  
  
  // removing repeated entries from 1_2list 1_3, 1_4 list of atoms and sorting them
   //std::cout<<" here "<<nImpDihedrals<<std::endl;
  for(int i=0; i<nAtoms; i++){
    atoms[i].list1_2.sort();
    atoms[i].list1_2.unique();
    
    atoms[i].list1_3.sort();
    atoms[i].list1_3.unique();
    
    atoms[i].list1_4.sort();
    atoms[i].list1_4.unique();
  }
  //  std::cout<<" here "<<nImpDihedrals<<std::endl; 
  
  /**
   * End of readin psf file    
   * Note that the cross terms (CMAP) are not implemented in this version
   * 
   */
}
//========================================================================



//========================================================================
//  readinprm
/**
   * reads in parameter file and assigen parameter for bond, angle, diheral, impdihedral
   * and nonbonded potential
   * 
   */

void mdpack::md::readinprm( const char * prm){
  
  bool hitBondSection = false;
  bool hitAngleSection = false;
  bool hitDihedralSection = false;
  bool hitImproperSection = false;
  bool hitCMAPsection = false;
  bool hitNonBondedSection= false;
  bool endOfNonBonded = false;
  
  /*
   * at the end of reading prm file, if checkBonds == nBonds; 
   * then we have found parameters for all bonds
   * otherwise terminate the program
   */
  int checkBonds=0;  
  int checkAngles=0; 
  int checkDihedrals =0;
  int checkImpropers =0;
  int checkNonBonded =0;
  std::string line;
  char part1 [50];
  //char part2 [50];
  
  char str1 [10];
  char str2 [10];
  char str3 [10];
  char str4 [10];
  const double zero =0.0;
  double unused;
  int nParts;
  // bond parameters
  mdpack::dreal kb;
  mdpack::dreal b0;
  
  // angle parameters
  mdpack::dreal k_theta;
  mdpack::dreal theta_0;
  mdpack::dreal k_ub;
  mdpack::dreal s_0;
  
  // dihedral parameters
  mdpack::dreal k_chi;
  mdpack::dreal multip ;
  mdpack::dreal delta ;
  
  // improper dihedral parameters
  mdpack::dreal k_psi;
  mdpack::dreal psi_0 ;
  
  // nonbonded parameters
  mdpack::dreal epsilon;
  mdpack::dreal rmin ;
   
  
  std::ifstream input;
  input.open(prm, std::ios::in);
  if (!input)
    terminate(" cant open prm file ");
  
  while(input >> part1 && !hitBondSection ){
    if(part1[0]=='!' || part1[0]=='*' || strcmp(part1, "MASS")==0 
       || strcmp(part1, "ATOMS")==0)
      input.ignore(250,'\n');
    else if( strcmp(part1, "BONDS")==0)
      hitBondSection = true;
    else
      terminate("  cannot find BONDS in prm file  ");
  } // end of while
  
  while(!input.eof() && !hitAngleSection){
    getline(input, line);
    nParts = countParts (line);
    if(nParts ==4){
      readBond(line, str1, str2, kb, b0);
      for(int i=0; i<nBonds; i++){
	if(!bonds[i].isSet){
	  if(bonds[i].compare(str1, str2)){
	    bonds[i].isSet=true;
	    bonds[i].setPrm(kb, b0);
	    checkBonds +=1;
	  }
	}
      }
    }
    else if(nParts ==1){
      if(line.compare("ANGLES")==0)
	hitAngleSection = true;
      else
	terminate("  cannot find ANGLES in prm file  ");
    }
  } // end of while
  
  if(checkBonds != nBonds)
    terminate("  cannot find parameters for bonds  ");
  
  while(!input.eof() && !hitDihedralSection){
    getline(input,line);
    nParts = countParts (line);
    if(nParts ==5){
      readAngle(line, str1, str2, str3, k_theta, theta_0);
      for(int i=0; i< nAngles ; i++){
	if(!angles[i].isSet){
	  if(angles[i].compare(str1, str2, str3)){
	    angles[i].isSet = true;
	    angles[i].setPrm(k_theta, theta_0, zero, zero);
	    checkAngles +=1;
	  }
	}
      }
    }
    
    else if( nParts == 7){
      readAngle(line, str1, str2, str3, k_theta, theta_0, k_ub, s_0);
      for(int i=0; i< nAngles ; i++){
	if(!angles[i].isSet){
	  if(angles[i].compare(str1, str2, str3)){
	    angles[i].isSet = true;
	    angles[i].setPrm(k_theta, theta_0, k_ub, s_0);
	    checkAngles +=1;
	  }
	}
      }
    }
    
    else if(nParts ==1){
      if(line.compare("DIHEDRALS")==0)
	hitDihedralSection = true;
      else
	terminate("  cannot find DIHEDRALS in prm file  ");
    }
  } // end of while
  
  if(checkAngles != nAngles)
    terminate("  cannot find parameters for Angles  ");
  
  //------------ reading dihedral parameres---------------------
  
  while(!input.eof() && !hitImproperSection){
    getline(input,line);
    nParts = countParts (line);
    if(nParts ==7){
      readDihedral(line, str1, str2, str3, str4,  k_chi, multip, delta);
      for(int i=0; i< nDihedrals ; i++){
	if(!dihedrals[i].isSet){
	  if(dihedrals[i].compare(str1, str2, str3, str4)){
	    dihedrals[i].isSet = true;
	    dihedrals[i].setPrm(k_chi,multip, delta);
	    checkDihedrals +=1;
	  }
	}
      }
    }
    
    else if(nParts ==1){
      if(line.compare("IMPROPER")==0)
	hitImproperSection = true;
      else
	terminate("  cannot find IMPROPER in prm file  ");
    }
  } // end of while
  
  
  if(checkDihedrals != nDihedrals)
    terminate("  cannot find parameters for Dihedrals  ");
  
  // -------------- reading improper dihedral parameters ---------------
  
  while(!input.eof() && !hitCMAPsection){
    getline(input,line);
    nParts = countParts (line);
    if(nParts ==7){
      readImproper(line, str1, str2, str3, str4,  k_psi, unused, psi_0);
      for(int i=0; i< nImpDihedrals ; i++){
	if(!impDihedrals[i].isSet){
	  if(impDihedrals[i].compare(str1, str2, str3, str4)){
	    impDihedrals[i].isSet = true;
	    impDihedrals[i].setPrm(k_psi, psi_0);
	    checkImpropers +=1;
	  }
	}
      }
    }
    
    else if(nParts ==1){
      if(line.compare("CMAP")==0)
	hitCMAPsection = true;
      else
	terminate("  cannot find CMAP in prm file  ");
    }
  } // end of while
  
  if(checkImpropers != nImpDihedrals)
    terminate("  cannot find parameters for ImpDihedrals  ");
  
  // skipping untill we reach the nonbonded section
  while(input >> part1 && !hitNonBondedSection){
    if(strcmp(part1,"cutnb")==0){
      input.ignore(250,'\n');
      hitNonBondedSection = true;
    }
  }
  
  // -----------------  reading nonbonded parameres ---------------------
  
  while( !input.eof() && !endOfNonBonded){
    getline(input,line);
    nParts = countParts (line);
    if(nParts > 3){
      readNonBonded(line, str1, unused, epsilon, rmin);
      for(int i=0; i<nAtoms ; i++){
	if(!atoms[i].isSet){
	  if(strcmp(atoms[i].type, str1)==0){
	    atoms[i].isSet = true;
	    atoms[i].setPrm(epsilon, rmin);
	    checkNonBonded +=1;
	  }
	}
      }
    }
    
    else if(nParts == 3){
      if( getWord(line,part1, 1) && strcmp(part1, "HBOND")==0)
	endOfNonBonded = true;
      else
	terminate("  cannot find HBOND in prm file  ");
    }
  } // end of while
  
  if(checkNonBonded != nAtoms)
    terminate("  cannot find parameters for NonBonded  ");
	
}// end of readinprm
//=====================================================================	



//========================================================================
//  readinpdb
/**
   * reads in pdb file and assigen coordinates for atoms
   * 
   */
void mdpack::md::readinpdb( const char * pdb){
 
  int checkAtoms=0;  
  
  std::string line;
  char part1 [50];
  char part2 [50];
  int nParts, id;
  mdpack::dreal x,y,z;
  mdpack::dreal tFactor, occup;
  
  std::ifstream input;
  input.open(pdb, std::ios::in);
  if (!input)
    terminate(" cannot open pdb file ");
  
  while(!input.eof()){
    getline(input,line);
    if(getWord(line,part1, 1) && strcmp(part1,"ATOM")==0){
      nParts = countWords(line);
      if(nParts == 13){
	getWord(line,part2, 2);
	id= atoi(part2);
	getWord(line,part2, 7);
	x= atof(part2);
	getWord(line,part2, 8);
	y= atof(part2);
	getWord(line,part2, 9);
	z= atof(part2);
	getWord(line,part2, 10);
	occup= atof(part2);
	getWord(line,part2, 11);
	tFactor= atof(part2);
      }
      else if (nParts == 12){
	getWord(line,part2, 2);
	id= atoi(part2);
	getWord(line,part2, 6);
	x= atof(part2);
	getWord(line,part2, 7);
	y= atof(part2);
	getWord(line,part2, 8);
	z= atof(part2);
	getWord(line,part2, 9);
	occup= atof(part2);
	getWord(line,part2, 10);
	tFactor= atof(part2);
      }
      else
	terminate(" cannot read pdb file ");
      atoms[id-1].setCoords(x,y,z);
      atoms[id-1].setTempFactor(tFactor);
      atoms[id-1].setOccupancy(occup);
      checkAtoms +=1;
    }
  }
  
  if(checkAtoms != nAtoms)
    terminate("  cannot find coordinates for atoms from PDB file ");
}

//========================================================================
// isIn1_2 
/**
   * Returns: 
   * \true, if a pairs of atoms are connected by a bond
   * \flase, otherwise
   * \c input:
   * int id1, id2
   * where id1 id the atomId for atom1 and id2 is the atomId for 
   * atom2
   */
bool mdpack::md::isIn1_2(const int& id1, const int& id2){
  return atoms[id1-1].isIn1_2(id2);
  
//   bool temp= false;
//   for(int i=0; i< nBonds; i++){
//     if(bonds[i].atom1 ==id1 && bonds[i].atom2 == id2)
//       temp= true;
//     else if(bonds[i].atom2 ==id1 && bonds[i].atom1 == id2)
//       temp= true;
//   }
//   return temp;
}
//========================================================================
// isIn1_3 
/**
   * Returns 
   * \true, if a pairs of atoms are connected by a bond
   *  or by two consecutive bonds (that is both are connceted by a bond  to
   *  a middle atom)
   * \false, otherwise
   * \c input:
   * int id1, id2
   * where id1 id the atomId for atom1 and id2 is the atomId for 
   * atom2 
   */
bool mdpack::md::isIn1_3(const int& id1, const int& id2){
  return ( atoms[id1-1].isIn1_2(id2) || atoms[id1-1].isIn1_3(id2) );
//   if(isIn1_2(id1,id2))
//     return true;
//   else{
//     bool temp =false;
//     for(int i=0; i< nAngles; i++){
//       if(angles[i].atom1 ==id1 && angles[i].atom3 == id2)
// 	temp = true;
//       else if(angles[i].atom3 ==id1 && angles[i].atom1 == id2)
// 	temp= true;
//     }
//     return temp;
//   }
}

//========================================================================
// isIn1_4
/**
   * Returns 
   * \true, if a pairs of atoms are connected by a bond
   *  or by two consecutive bonds (that is both are connceted by a bond  to
   *  a middle atom) or three consecutive bonds
   * \false, otherwise
   * \c input:
   * int id1, id2
   * where id1 id the atomId for atom1 and id2 is the atomId for 
   * atom2 
   */
bool mdpack::md::isIn1_4(const int& id1, const int& id2){
  return ( atoms[id1-1].isIn1_2(id2) || atoms[id1-1].isIn1_3(id2) 
          || atoms[id1-1].isIn1_4(id2));

//   if(isIn1_2(id1, id2))
//     return true;
//   else if(isIn1_3(id1,id2))
//     return true;
//   else{
//     bool temp=false;
//     for(int i=0; i< nDihedrals; i++){
//       if(dihedrals[i].atom1 ==id1 && dihedrals[i].atom4 == id2)
// 	temp =true;
//       else if(dihedrals[i].atom4 ==id1 && dihedrals[i].atom1 == id2)
// 	temp =true;
//     }
//     for(int i=0; i< nImpDihedrals; i++){
//       if(impDihedrals[i].atom1 ==id1 && impDihedrals[i].atom4 == id2)
// 	temp= true;
//       else if(impDihedrals[i].atom4 ==id1 && impDihedrals[i].atom1 == id2)
// 	temp= true;
//     }
//     return temp;
//   }
}


//========================================================================
//  setLJtable
/**
   * setting the parameres A, B for Lennard-Jones potential
   * also excludes the pairs using the bonded-exclusion option by
   * setting paramere inc into 0 (if the pair is excluded)
   * 
   */

void mdpack::md::setLJtable(){
  mdpack::dreal eps, sigma;
  nonBondedPairs = new  mdpack::ljTable * [nAtoms];
  for (int i = 0; i < nAtoms; ++i)
     nonBondedPairs[i] = new mdpack::ljTable [nAtoms];
  isNonBondedPairs = true;
  
  
  for (int i=0; i< nAtoms-1 ; i++)
    for (int j=i+1; j<nAtoms ; j++){
      
      if(strcmp(exclude,"1-2")==0){
	if(!isIn1_2(atoms[i].id, atoms[j].id)){
	  nonBondedPairs[i][j].inc=1;
	  sigma = atoms[i].halfSigma + atoms[j].halfSigma;
	  eps= sqrt(atoms[i].eps*atoms[j].eps);
	  (nonBondedPairs[i][j]).A= 4.0*pow12(sigma)*eps;
	  (nonBondedPairs[i][j]).B= 4.0*pow6(sigma)*eps;
	  (nonBondedPairs[i][j]).id1= atoms[i].id;
	  (nonBondedPairs[i][j]).id2= atoms[j].id;
	//  nonBondedPairsList.push_back(nonBondedPairs[i][j]);
	  nonBondedPairs[j][i] = nonBondedPairs[i][j];
	}
      }
      
      else if(strcmp(exclude,"1-3")==0){
	if(!isIn1_3(atoms[i].id, atoms[j].id)){
	  nonBondedPairs[i][j].inc=1;
	  sigma = atoms[i].halfSigma + atoms[j].halfSigma;
	  eps= sqrt(atoms[i].eps*atoms[j].eps);
	  (nonBondedPairs[i][j]).A= 4.0*pow12(sigma)*eps;
	  (nonBondedPairs[i][j]).B= 4.0*pow6(sigma)*eps;
	  (nonBondedPairs[i][j]).id1= atoms[i].id;
	  (nonBondedPairs[i][j]).id2= atoms[j].id;
	//  nonBondedPairsList.push_back(nonBondedPairs[i][j]);
	  nonBondedPairs[j][i] = nonBondedPairs[i][j] ;
	}
      }
      
      else if(strcmp(exclude,"1-4")==0){
	if(!isIn1_4(atoms[i].id, atoms[j].id)){
// 	  if(i> nAtoms-1 || j> nAtoms -1){
// 	    std::cout<<" i "<< i <<" j "<<j<<std::endl;
// 	    exit(0);
// 	  }
	  nonBondedPairs[i][j].inc=1;
	  sigma = atoms[i].halfSigma + atoms[j].halfSigma;
	  eps= sqrt(atoms[i].eps*atoms[j].eps);
	  (nonBondedPairs[i][j]).A= pow12(sigma)*eps;
 	  (nonBondedPairs[i][j]).B= 2.0*pow6(sigma)*eps;
	  (nonBondedPairs[i][j]).id1= atoms[i].id;
 	  (nonBondedPairs[i][j]).id2= atoms[j].id;
	  //mdpack::ljTable tt = nonBondedPairs[j][i];
	  // std::cout<<nonBondedPairs[i][j].id1<<std::endl;
 	 // nonBondedPairsList.push_back(nonBondedPairs[j][i] );
	  nonBondedPairs[j][i] = nonBondedPairs[i][j];
	}
      }
      
      else{
	nonBondedPairs[i][j].inc=1;
	sigma = atoms[i].halfSigma + atoms[j].halfSigma;
	eps= sqrt(atoms[i].eps*atoms[j].eps);
	(nonBondedPairs[i][j]).A= 4.0*pow12(sigma)*eps;
	(nonBondedPairs[i][j]).B= 4.0*pow6(sigma)*eps;
	(nonBondedPairs[i][j]).id1= atoms[i].id;
	(nonBondedPairs[i][j]).id2= atoms[j].id;
	//nonBondedPairsList.push_back(nonBondedPairs[i][j]);
	nonBondedPairs[j][i] = nonBondedPairs[i][j];
      }
    }
}

/**============================================================
 * readinConfig: reading and setting parameter from config file
 */ 

void mdpack::md::readinConfig(const char * conf){
  
  std::ifstream input;
  input.open(conf, std::ios::in);
  if (!input) {
  std::cerr <<  "Can't open config file : "<< conf << std::endl;
  exit(1);
  }
  char part[50];

  while(input >> part){
    if(part[0]=='#')
      input.ignore(250,'\n');
    else{
      if(strcmp(part, "structure")==0)
	input >> psfFile;
      else if(strcmp(part, "coordinates")==0)
	input >> pdbFile;
      else if(strcmp(part, "parameters")==0)
	input >> prmFile;
      else if(strcmp(part, "numSteps")==0){
	input >> numSteps;
      }
      // setting timeStep from config file, must be given in ps
      else if(strcmp(part, "timeStep")==0){
	input >> h;
	timeStep = h;
	h= TIMEFACTOR*h; // mdpack time unit, so the kinetic energy will be Kcal/mol
	
      }
      else if(strcmp(part, "initialTime")==0){
	input >> initialTime;
	simTime = initialTime;
      }
      else if(strcmp(part, "switching")==0){
	input >> part;
	if(strcmp(part,"on")==0)
	  switching = true;
	else if(strcmp(part,"off")==0)
	  switching = false;
	else
	  std::cerr <<"cannot read parameter for switching"<<
	  " from config file "<< std::endl;
      }
      else if(strcmp(part, "cutoff")==0)
	input >> cutoff;
      else if(strcmp(part, "switchdist")==0)
	input >> switchdist;
      else if(strcmp(part, "exclude")==0){
	input >> part;
	if(strcmp(part, "none")==0)
	  strcpy(exclude,"none");
	else if(strcmp(part, "1-2")==0)
	  strcpy(exclude,"1-2");
	else if(strcmp(part, "1-3")==0)
	  strcpy(exclude,"1-3");
	else if(strcmp(part, "1-4")==0)
	  strcpy(exclude,"1-4");
	else
	  std::cerr <<"cannot read parameter for exclude"<<
	  " from config file "<< std::endl;
      }
      else if(strcmp(part, "dielectric")==0)
	input >> e0;
      else if(strcmp(part, "sphericalBC")==0){
	input >> part;
	if(strcmp(part, "on")==0)
	  sphBC = true;
	else if(strcmp(part, "off")==0)
	  sphBC = false;
	else
	  std::cerr <<"cannot read parameter for sphericalBC "<<
	  " from config file, acceptable values are on or off"<< std::endl;
      }
      else if(strcmp(part, "sphericalBCcenter")==0){
	input >> sphBCCenter.x >> sphBCCenter.y >> sphBCCenter.z;
      }
      else if(strcmp(part, "sphericalBCr1")==0)
	input >> sphBCr1;
      else if(strcmp(part, "sphericalBCk1")==0)
	input >> sphBCk1;
      else if(strcmp(part, "sphericalBCexp1")==0)
	input >> sphBCexp1;  
      
      // output parameters  
      else if(strcmp(part, "outputName")==0)
	input >> outputName;
      else if(strcmp(part, "writePosition")==0){
	input >> part;
	if(strcmp(part, "on")==0)
	  isWritePosOn = true;
	else if(strcmp(part, "off")==0)
	  isWritePosOn = false;
	else
	  std::cerr <<"cannot read parameter writePosition "<<
	  " from config file, acceptable values are on or off"<< std::endl;
      }
      else if(strcmp(part, "writePosition")==0){
	input >> part;
	if(strcmp(part, "on")==0)
	  isWritePosOn = true;
	else if(strcmp(part, "off")==0)
	  isWritePosOn = false;
	else
	  std::cerr <<"cannot read parameter writePosition "<<
	  " from config file, acceptable values are on or off"<< std::endl;
      }
      else if(strcmp(part, "writeVelocities")==0){
	input >> part;
	if(strcmp(part, "on")==0)
	  isWriteVelOn = true;
	else if(strcmp(part, "off")==0)
	  isWriteVelOn = false;
	else
	  std::cerr <<"cannot read parameter writeVelocities "<<
	  " from config file, acceptable values are on or off"<< std::endl;
      }
      else if(strcmp(part, "writePhi")==0){
	input >> part;
	if(strcmp(part, "on")==0)
	  isWritePhiOn = true;
	else if(strcmp(part, "off")==0)
	  isWritePhiOn = false;
	else
	  std::cerr <<"cannot read parameter writePhi "<<
	  " from config file, acceptable values are on or off"<< std::endl;
      }
      else if(strcmp(part, "binary")==0){
	input >> part;
	if(strcmp(part, "on")==0)
	  isBinaryOn = true;
	else if(strcmp(part, "off")==0)
	  isBinaryOn = false;
	else
	  std::cerr <<"cannot read parameter binary "<<
	  " from config file, acceptable values are on or off"<< std::endl;
      }
      else if(strcmp(part, "cellLists")==0){
	input >> part;
	if(strcmp(part, "on")==0)
	  useCellLists = true;
	else if(strcmp(part, "off")==0)
	  useCellLists = false;
	else
	  std::cerr <<"cannot read parameter cellLists "<<
	  " from config file, acceptable values are on or off"<< std::endl;
      }
      
      else if(strcmp(part, "seed")==0){
	input >> seed;
	W.setSeed(seed);
      }
      else if(strcmp(part, "gamma")==0){
	input >> gamma;
      }
      else if(strcmp(part, "temperature")==0){
	input >> targetTemp;
      }
      
      else if(strcmp(part, "NHLgamma")==0){
	input >> nhl_gamma;
      }
      else if(strcmp(part, "NHLmu")==0){
	input >> nhl_mu;
	nhl_imu= 1.0/nhl_mu;
      }
	
      
    }// end of else
  } // end of while loop 
}// end of readinConfig
/**=========================================================*/  


/**============================================================
 * computeNonBondeds: computes forces and energies for Lennard-Jones
 * and electorstatic interactions
 */ 
void mdpack::md::computeNonBondeds(void){
 // std::list<mdpack::ljTable>::iterator it;
  mdpack::dreal A, B;
  int inc;
  Evdw =0.0;
  Eelec=0.0;
  //std::cout<<Eelec <<std::endl;
 // std::cout<< nAtoms <<std::endl;
 // std::cout<< MINATOMS <<std::endl;
 // if(!useCellLists || nAtoms < MINATOMS)
   // std::cout<< " true " <<std::endl;
    
  if(!useCellLists || nAtoms < MINATOMS){
    for(int i =0; i<nAtoms-1; i++)
      for(int j=i+1; j<nAtoms; j++){
    //for(it=nonBondedPairsList.begin(); it!= nonBondedPairsList.end(); ++it){
    //  std::cout<< " true " <<std::endl;
      
	idx1 = i;
	idx2 = j;
	//std::cout<<idx1 << " "<<idx2<<std::endl;
	
	A = (nonBondedPairs[idx1][idx2]).A;
        B = (nonBondedPairs[idx1][idx2]).B;
        inc = (nonBondedPairs[idx1][idx2]).inc;
	//std::cout<<inc <<std::endl;
	
	if(inc !=0){
	//std::cout<<" idx1 "<<idx1<< " idx2 "<<idx2<<std::endl;
	// to save computing resources we write everything explicitly
	rij.x= atoms[idx2].r.x - atoms[idx1].r.x;
	rij.y= atoms[idx2].r.y - atoms[idx1].r.y;
	rij.z= atoms[idx2].r.z - atoms[idx1].r.z;
    
	qi = atoms[idx1].q;
	qj = atoms[idx2].q;
	dist  = sqrt(rij.x*rij.x + rij.y*rij.y + rij.z*rij.z);
	if(!switching){
	  inv_dist= 1.0/dist;
	  inv_dist6= (inv_dist*inv_dist*inv_dist*inv_dist*inv_dist*inv_dist);
	  inv_dist12 = inv_dist6*inv_dist6;
	  Evdw += ( A*inv_dist12 - B*inv_dist6 );
	  factor2 = e14*COULOMB*qi*qj*inv_e0*inv_dist;
	  Eelec += factor2;
 
	  dV = ( 12.0*A*inv_dist12*inv_dist - 6.0*B*inv_dist6*inv_dist );
	
	  f[idx1].x -=  inv_dist*dV*rij.x;
	  f[idx1].y -=  inv_dist*dV*rij.y;
	  f[idx1].z -=  inv_dist*dV*rij.z;
	  f[idx2].x +=  inv_dist*dV*rij.x; 
	  f[idx2].y +=  inv_dist*dV*rij.y; 
	  f[idx2].z +=  inv_dist*dV*rij.z; 
     
	  dV= factor2*inv_dist;
    
	  f[idx1].x -=  inv_dist*dV*rij.x;
	  f[idx1].y -=  inv_dist*dV*rij.y;
	  f[idx1].z -=  inv_dist*dV*rij.z;
	  f[idx2].x +=  inv_dist*dV*rij.x; 
	  f[idx2].y +=  inv_dist*dV*rij.y; 
	  f[idx2].z +=  inv_dist*dV*rij.z; 
	
	} //end of if(!switching)
      
	else{
    
	  if(dist <= switchdist){
	    inv_dist= 1.0/dist;
	    inv_dist6= (inv_dist*inv_dist*inv_dist*inv_dist*inv_dist*inv_dist);
	    inv_dist12 = inv_dist6*inv_dist6;
	    Evdw += ( A*inv_dist12 - B*inv_dist6 );
	    factor1 = 1.0-dist*dist*inv_Rc2;
	    factor2 = e14*COULOMB*qi*qj*inv_e0*inv_dist;
	    Eelec += (factor2)*(factor1*factor1);
 
	    dV = ( 12.0*A*inv_dist12*inv_dist - 6.0*B*inv_dist6*inv_dist );
	
	    f[idx1].x -=  inv_dist*dV*rij.x;
	    f[idx1].y -=  inv_dist*dV*rij.y;
	    f[idx1].z -=  inv_dist*dV*rij.z;
	    f[idx2].x +=  inv_dist*dV*rij.x; 
	    f[idx2].y +=  inv_dist*dV*rij.y; 
	    f[idx2].z +=  inv_dist*dV*rij.z; 
     
	    dV= ( (factor2*inv_dist)*(factor1*factor1) ) +
	      ( 4.0*(dist*inv_Rc2)*(factor1)*(factor2) ); 
    
	    f[idx1].x -=  inv_dist*dV*rij.x;
	    f[idx1].y -=  inv_dist*dV*rij.y;
	    f[idx1].z -=  inv_dist*dV*rij.z;
	    f[idx2].x +=  inv_dist*dV*rij.x; 
	    f[idx2].y +=  inv_dist*dV*rij.y; 
	    f[idx2].z +=  inv_dist*dV*rij.z; 
	
	  }
    
	  else if( dist < cutoff){
	    inv_dist= 1.0/dist;
	    inv_dist6= (inv_dist*inv_dist*inv_dist*inv_dist*inv_dist*inv_dist);
	    inv_dist12 = inv_dist6*inv_dist6;
      
	
	    factor1 = cutoff*cutoff-dist*dist;
	    factor2 = cutoff*cutoff+2.0*dist*dist -3.0*switchdist*switchdist;
	    factor3 = 1.0-dist*dist*inv_Rc2;
	    factor4 = e14*COULOMB*qi*qj*inv_e0*inv_dist;
	    Evdw += ( A*inv_dist12 - B*inv_dist6 )*
		  ( fct*(factor1*factor1)*(factor2) );
	    Eelec += (factor4)*(factor3*factor3); 
    
      
	    dV = ( 12.0*A*inv_dist12*inv_dist - 6.0*B*inv_dist6*inv_dist )*
		( fct*(factor1*factor1)*(factor2) ) +
		( A*inv_dist12 - B*inv_dist6  )*
		( 4.0*fct*( factor1*(factor2) - dist*factor1*factor1  ) );   

	    f[idx1].x -=  inv_dist*dV*rij.x;
	    f[idx1].y -=  inv_dist*dV*rij.y;
	    f[idx1].z -=  inv_dist*dV*rij.z;
	    f[idx2].x +=  inv_dist*dV*rij.x; 
	    f[idx2].y +=  inv_dist*dV*rij.y; 
	    f[idx2].z +=  inv_dist*dV*rij.z; 
      
	    dV= ( (factor4*inv_dist)*(factor3*factor3) ) +
	      ( 4.0*(dist*inv_Rc2)*(factor3)*(factor4) );
     
	    f[idx1].x -=  inv_dist*dV*rij.x;
	    f[idx1].y -=  inv_dist*dV*rij.y;
	    f[idx1].z -=  inv_dist*dV*rij.z;
	    f[idx2].x +=  inv_dist*dV*rij.x; 
	    f[idx2].y +=  inv_dist*dV*rij.y; 
	    f[idx2].z +=  inv_dist*dV*rij.z; 
	  }
	} // end else
      } // if(inc !=0)
    } // end of loop over nonBondedPairs
  }// end of nonBondedPairs without cellList
  
    else{
   // std::cout<< " else "<<std::endl;
//     for(int i=0; i<nCells ; i++)
//       cells[i].update(atoms);
//     for(int i=0; i<nCells; i++)
//       cells[i].update(cells, atoms);
      for(int i=0; i<nCells; i++)
	cells[i].atomIds.clear();

      for (int i=0; i<nAtoms; i++){
	//atoms[i].r.print();
	
	  for(int j=0; j<nCells; j++){
// 	     std::cout<<cells[j].leftx<<" "<<cells[j].rightx<<" "<< 
// 	               cells[j].lefty<<" "<<cells[j].righty<<" "<< 
// 	               cells[j].leftz<<" "<<cells[j].rightz<<" "<<std::endl;
	    if(cells[j].isIn(atoms[i].r))
	      cells[j].addAtom(i+1);
	  }
      }
	  int allAtoms=0.0;
	  for(int i=0; i<nCells; i++)
            allAtoms += cells[i].getSize();
	  //  std::cout<<" all atoms "<< allAtoms << "  "<<nCells<<std::endl;
	//  std::cout<<" all atoms "<< allAtoms << "  "<<nCells<<std::endl;
      std::list<int>::iterator it1;
      std::list<int>::iterator id1;
      std::list<int>::iterator id2;
    //std::cout<<" here "<<std::endl;
      for(int i=0; i<nCells; i++)
	//if(cells[i].getSize()>0)
      //std::cout<<" here "<<std::endl;
	  for(it1=cells[i].neighbours.begin(); it1!=cells[i].neighbours.end(); ++it1){
	   // if(cells[*it1].getSize()>0)
	  // std::cout<<" here "<<std::endl;
	      for(id1=cells[i].atomIds.begin(); id1!=cells[i].atomIds.end(); ++id1)
		for(id2=cells[*it1].atomIds.begin(); id2!=cells[*it1].atomIds.end(); ++id2){
		//std::cout<<" here "<<std::endl;
		  idx1 =  *id1 -1;
		  idx1 =  *id2 -1;
		  A = (nonBondedPairs[idx1][idx2]).A;
		  B = (nonBondedPairs[idx1][idx2]).B;
		  inc = (nonBondedPairs[idx1][idx2]).inc;
		  if(inc !=0){
		    rij.x= atoms[idx2].r.x - atoms[idx1].r.x;
		    rij.y= atoms[idx2].r.y - atoms[idx1].r.y;
		    rij.z= atoms[idx2].r.z - atoms[idx1].r.z;
    
		    qi = atoms[idx1].q;
		    qj = atoms[idx2].q;
		    dist  = sqrt(rij.x*rij.x + rij.y*rij.y 
			    + rij.z*rij.z);
    
		    if(dist <= switchdist){
		      inv_dist= 1.0/dist;
		      inv_dist6= (inv_dist*inv_dist*inv_dist*inv_dist*inv_dist*inv_dist);
		      inv_dist12 = inv_dist6*inv_dist6;
		      Evdw += ( A*inv_dist12 - B*inv_dist6 );
		      factor1 = 1.0-dist*dist*inv_Rc2;
		      factor2 = e14*COULOMB*qi*qj*inv_e0*inv_dist;
		      Eelec += (factor2)*(factor1*factor1);
 
		      dV = ( 12.0*A*inv_dist12*inv_dist - 6.0*B*inv_dist6*inv_dist );
	
		      f[idx1].x -=  inv_dist*dV*rij.x;
		      f[idx1].y -=  inv_dist*dV*rij.y;
		      f[idx1].z -=  inv_dist*dV*rij.z;
		      f[idx2].x +=  inv_dist*dV*rij.x; 
		      f[idx2].y +=  inv_dist*dV*rij.y; 
		      f[idx2].z +=  inv_dist*dV*rij.z; 
     
		      dV= (( (factor2*inv_dist)*(factor1*factor1) ) +
			  ( 4.0*(dist*inv_Rc2)*(factor1)*(factor2) ) ); 
    
		      f[idx1].x -=  inv_dist*dV*rij.x;
		      f[idx1].y -=  inv_dist*dV*rij.y;
		      f[idx1].z -=  inv_dist*dV*rij.z;
		      f[idx2].x +=  inv_dist*dV*rij.x; 
		      f[idx2].y +=  inv_dist*dV*rij.y; 
		      f[idx2].z +=  inv_dist*dV*rij.z; 
	
		    }
    
		    else if( dist < cutoff){
		      inv_dist= 1.0/dist;
		      inv_dist6= (inv_dist*inv_dist*inv_dist*inv_dist*inv_dist*inv_dist);
		      inv_dist12 = inv_dist6*inv_dist6;
      
	  
		      factor1 = cutoff*cutoff-dist*dist;
		      factor2 = cutoff*cutoff+2.0*dist*dist -3.0*switchdist*switchdist;
		      factor3 = 1.0-dist*dist*inv_Rc2;
		      factor4 = e14*COULOMB*qi*qj*inv_e0*inv_dist;
		      Evdw += ( A*inv_dist12 - B*inv_dist6 )*
			  ( fct*(factor1*factor1)*(factor2) );
		      Eelec += (factor4)*(factor3*factor3); 
    
      
		      dV =( 12.0*A*inv_dist12*inv_dist - 6.0*B*inv_dist6*inv_dist )*
			  ( fct*(factor1*factor1)*(factor2) ) +
			  ( A*inv_dist12 - B*inv_dist6  )*
			  ( 4.0*fct*( factor1*(factor2) - dist*factor1*factor1  ) );   

		      f[idx1].x -=  inv_dist*dV*rij.x;
		      f[idx1].y -=  inv_dist*dV*rij.y;
		      f[idx1].z -=  inv_dist*dV*rij.z;
		      f[idx2].x +=  inv_dist*dV*rij.x; 
		      f[idx2].y +=  inv_dist*dV*rij.y; 
		      f[idx2].z +=  inv_dist*dV*rij.z; 
      
		      dV= ( (factor4*inv_dist)*(factor3*factor3) ) +
			  ( 4.0*(dist*inv_Rc2)*(factor3)*(factor4) );
     
		      f[idx1].x -=  inv_dist*dV*rij.x;
		      f[idx1].y -=  inv_dist*dV*rij.y;
		      f[idx1].z -=  inv_dist*dV*rij.z;
		      f[idx2].x +=  inv_dist*dV*rij.x; 
		      f[idx2].y +=  inv_dist*dV*rij.y; 
		      f[idx2].z +=  inv_dist*dV*rij.z; 
		    }
		  } // if(inc!=0)
		}// end of for
	  }
      }// end of else for cell list
		   
}
/**==================================================================*/
   
/**============================================================
 * computeBonds: computes bond forces and bond energies for all
 * bonds
 */
void mdpack::md::computeBonds(void){
  Ebonds =0.0;
  for(int i=0; i< nBonds; i++){
    idx1 = bonds[i].atom1-1;
    idx2 = bonds[i].atom2-1;
    rij.x = atoms[idx2].r.x - atoms[idx1].r.x;
    rij.y = atoms[idx2].r.y - atoms[idx1].r.y;
    rij.z = atoms[idx2].r.z - atoms[idx1].r.z;
    dist  = sqrt(rij.x*rij.x + rij.y*rij.y + rij.z*rij.z);
    
    factor1 = dist - bonds[i].r0;
    factor2 = 2.0*bonds[i].k*(factor1)/dist;
    f[idx1].x += factor2*rij.x;
    f[idx1].y += factor2*rij.y;
    f[idx1].z += factor2*rij.z;
    f[idx2].x -= factor2*rij.x;
    f[idx2].y -= factor2*rij.y;
    f[idx2].z -= factor2*rij.z;
   
    Ebonds += bonds[i].k*(factor1*factor1);
  }
}
/**===========================================================*/
  
/**============================================================
 * computeAngles: computes angle forces and angle energies for all
 * angles
 * ( note there is a singularity in 1/sin, this should be fixed later,
 *   will do for now)
 */
void mdpack::md::computeAngles2(void){
  Eangles =0.0;
  for(int i=0; i< nAngles; i++){
    idx1 = angles[i].atom1-1;
    idx2 = angles[i].atom2-1;
    idx3 = angles[i].atom3-1;
    rij.x = atoms[idx2].r.x - atoms[idx1].r.x;
    rij.y = atoms[idx2].r.y - atoms[idx1].r.y;
    rij.z = atoms[idx2].r.z - atoms[idx1].r.z;
    rkj.x = atoms[idx2].r.x - atoms[idx3].r.x;
    rkj.y = atoms[idx2].r.y - atoms[idx3].r.y;
    rkj.z = atoms[idx2].r.z - atoms[idx3].r.z;
   
    //drij = inorm(rij);
    //drkj = inorm(rkj);
    
    drij= 1.0/sqrt(rij.x*rij.x + 
                        rij.y*rij.y + rij.z*rij.z);
   
    drkj = 1.0/sqrt(rkj.x*rkj.x + 
                        rkj.y*rkj.y + rkj.z*rkj.z);
    dot_ij_kj= (rij.x*rkj.x + rij.y*rkj.y + rij.z*rkj.z);
    cos_theta = (drij*drkj)*dot_ij_kj;
    if(cos_theta >1.0) cos_theta =1.0;
    if(cos_theta <-1.0) cos_theta = -1.0;
    atheta= acos(cos_theta);
    sin_theta = sin(atheta);
    factor2 = atheta- degTorad*angles[i].theta0;
    factor1 = (2.0*angles[i].ktheta*(factor2))/sin_theta;
    if(abs(sin_theta)<0.0001){
      std::cerr<< " sin(theta) is zero in angle potetial "<<i<<std::endl;
      exit(1);
    }
    
    //fi = (factor1/(sin_theta*drij))*( (rij*(cos_theta/drij)) -(rkj/drkj));
    fi = (drij*factor1)*( (rij*(cos_theta*drij)) -(rkj*drkj));
    f[idx1] += fi;
    //fk = (factor1/(sin_theta*drkj))* ( (rkj*(cos_theta/drkj)) -(rij/drij));
    fk = (drkj*factor1)* ( (rkj*(cos_theta*drkj)) -(rij*drij));
    f[idx3] += fk;
    f[idx2] -=(fi +fk);
    
    
//     fi.x = (drij*factor1)*( (rij.x*(cos_theta*drij)) -(rkj.x*drkj));
//     fi.y = (drij*factor1)*( (rij.y*(cos_theta*drij)) -(rkj.y*drkj));
//     fi.z = (drij*factor1)*( (rij.z*(cos_theta*drij)) -(rkj.z*drkj));
//     f[idx1].x += fi.x;
//     f[idx1].y += fi.y;
//     f[idx1].z += fi.z;
//     
//     fk.x = (drkj*factor1)* ( (rkj.x*(cos_theta*drkj)) -(rij.x*drij));
//     fk.y = (drkj*factor1)* ( (rkj.y*(cos_theta*drkj)) -(rij.y*drij));
//     fk.z = (drkj*factor1)* ( (rkj.z*(cos_theta*drkj)) -(rij.z*drij));
//     f[idx3].x += fk.x;
//     f[idx3].y += fk.y;
//     f[idx3].z += fk.z;
//     f[idx2].x -=(fi.x +fk.x);
//     f[idx2].y -=(fi.y +fk.y);
//     f[idx2].z -=(fi.z +fk.z);
    
    Eangles += angles[i].ktheta*pow2(atheta- degTorad*angles[i].theta0);
    if(angles[i].kub >0.0){
      //rik = atoms[idx3].r - atoms[idx1].r;
      rik.x = atoms[idx3].r.x - atoms[idx1].r.x;
      rik.y = atoms[idx3].r.y - atoms[idx1].r.y;
      rik.z = atoms[idx3].r.z - atoms[idx1].r.z;
      //drik = norm(rik);
      drik = sqrt( rik.x*rik.x + rik.y*rik.y+rik.z*rik.z);
      factor3= (1.0/drik)*(2.0*angles[i].kub*(drik-angles[i].rub));
      //f[idx1] += 2.0*angles[i].kub*(drik-angles[i].rub)*(rik/drik);
      f[idx1].x +=  factor3*rik.x;
      f[idx1].y +=  factor3*rik.y;
      f[idx1].z +=  factor3*rik.z; 
      //f[idx3] -= 2.0*angles[i].kub*(drik-angles[i].rub)*(rik/drik);
      f[idx3].x -=  factor3*rik.x;
      f[idx3].y -=  factor3*rik.y;
      f[idx3].z -=  factor3*rik.z;
      Eangles += angles[i].kub*pow2(drik-angles[i].rub);
    }
    
//     Eangles += angles[i].ktheta*(factor2*factor2);
//     if(angles[i].kub >0.0){
//       rik.x = atoms[idx3].r.x - atoms[idx1].r.x;
//       rik.y = atoms[idx3].r.y - atoms[idx1].r.y;
//       rik.z = atoms[idx3].r.z - atoms[idx1].r.z;
//       //drik = norm(rik);
//       drik = sqrt( rik.x*rik.x + rik.y*rik.y+rik.z*rik.z);
//       factor3= (2.0*angles[i].kub*(drik-angles[i].rub))/drik;
//       f[idx1].x +=  factor3*rik.x; 
//       f[idx1].y +=  factor3*rik.y; 
//       f[idx1].z +=  factor3*rik.z; 
//       f[idx3].x -=  factor3*rik.x; 
//       f[idx3].y -=  factor3*rik.y;
//       f[idx3].z -=  factor3*rik.z;
//       Eangles += angles[i].kub*(drik-angles[i].rub)*(drik-angles[i].rub);
//     }
  }
}
    
    
    
//     factor3 = (atheta- degTorad*angles[i].theta0);
//     factor1 = ( 2.0*angles[i].ktheta*(factor3) )/sin_theta;
//     
//     factor2 = cos_theta*inv_dist_ij;
//     fix = (inv_dist_ij*factor1)*( (rij.x*(factor2)) -(rkj.x*inv_dist_kj));
//     fiy = (inv_dist_ij*factor1)*( (rij.y*(factor2)) -(rkj.y*inv_dist_kj));
//     fiz = (inv_dist_ij*factor1)*( (rij.z*(factor2)) -(rkj.z*inv_dist_kj));
//     f[idx1].x += fix;
//     f[idx1].y += fiy;
//     f[idx1].z += fiz;
//     
//     factor2= cos_theta*inv_dist_kj;
//     fkx = (inv_dist_kj*factor1)* ( (rkj.x*factor2) -(rij.x*inv_dist_ij));
//     fky = (inv_dist_kj*factor1)* ( (rkj.y*factor2) -(rij.y*inv_dist_ij));
//     fkz = (inv_dist_kj*factor1)* ( (rkj.z*factor2) -(rij.z*inv_dist_ij));
//     f[idx3].x += fkx;
//     f[idx3].y += fky;
//     f[idx3].z += fkz;
//     f[idx2].x -=(fix +fkx);
//     f[idx2].y -=(fiy +fky);
//     f[idx2].z -=(fiz +fkz);
//     
//     Eangles += angles[i].ktheta*(factor3*factor3);
//     if(angles[i].kub >0.0){
//       rik.x = atoms[idx3].r.x - atoms[idx1].r.x;
//       rik.y = atoms[idx3].r.y - atoms[idx1].r.y;
//       rik.z = atoms[idx3].r.z - atoms[idx1].r.z;
//       dist_ik = sqrt(rik.x*rik.x + 
//                         rik.y*rik.y + rik.z*rik.z); 
//       factor5 = dist_ik-angles[i].rub;
//       factor4 = 2.0*angles[i].kub*(factor5)/dist_ik;
//       f[idx1].x += factor4*(rik.x); 
//       f[idx1].y += factor4*(rik.y);
//       f[idx1].z += factor4*(rik.z); 
//       f[idx3].x -= factor4*(rik.x);
//       f[idx3].y -= factor4*(rik.y);
//       f[idx3].z -= factor4*(rik.z);
//       
//       Eangles += angles[i].kub*(factor5)*(factor5);
//     }
//   }
// }
/**======================================================*/	


/**============================================================
 * computeAngles2: computes angle forces and angle energies for all
 * angles
 * ( note this version is singularity free)
 */
void mdpack::md::computeAngles(void){
  Eangles =0.0;
  for(int i=0; i< nAngles; i++){
    idx1 = angles[i].atom1-1;
    idx2 = angles[i].atom2-1;
    idx3 = angles[i].atom3-1;
    rji = atoms[idx1].r - atoms[idx2].r;
    rjk = atoms[idx3].r - atoms[idx2].r;
    inv_dist_rji = 1.0/norm(rji);
    inv_dist_rjk = 1.0/norm(rjk);
    //MM = dotProd(rji, rjk);
    cos_theta =  inv_dist_rji*inv_dist_rjk*dotProd(rji, rjk);;
    if(cos_theta>1.0)
      cos_theta =1.0;
    else if(cos_theta <-1.0)
      cos_theta =-1.0;
    atheta = acos(cos_theta);
    sin_theta = sin(atheta);
    factor1 = atheta-degTorad*angles[i].theta0 ;
    dV_dtheta = 2.0*angles[i].ktheta*factor1;
    Eangles += angles[i].ktheta*factor1*factor1;
    
    if(angles[i].kub >0.0){
      //rik = atoms[idx3].r - atoms[idx1].r;
      rik.x = atoms[idx3].r.x - atoms[idx1].r.x;
      rik.y = atoms[idx3].r.y - atoms[idx1].r.y;
      rik.z = atoms[idx3].r.z - atoms[idx1].r.z;
      //drik = norm(rik);
      drik = sqrt( rik.x*rik.x + rik.y*rik.y+rik.z*rik.z);
      factor3= (1.0/drik)*(2.0*angles[i].kub*(drik-angles[i].rub));
      //f[idx1] += 2.0*angles[i].kub*(drik-angles[i].rub)*(rik/drik);
      f[idx1].x +=  factor3*rik.x;
      f[idx1].y +=  factor3*rik.y;
      f[idx1].z +=  factor3*rik.z; 
      //f[idx3] -= 2.0*angles[i].kub*(drik-angles[i].rub)*(rik/drik);
      f[idx3].x -=  factor3*rik.x;
      f[idx3].y -=  factor3*rik.y;
      f[idx3].z -=  factor3*rik.z;
      Eangles += angles[i].kub*pow2(drik-angles[i].rub);
    }
   
   if( sin_theta >0.1 || sin_theta <-0.1){
     factor2 = dV_dtheta*(1.0/sin_theta);
     factor3 = inv_dist_rji* inv_dist_rjk;
     
     fi = factor2*( factor3*rjk - (cos_theta*inv_dist_rji*inv_dist_rji)*rji );
     fk = factor2*( factor3*rji - (cos_theta*inv_dist_rjk*inv_dist_rjk)*rjk );
     
     f[idx1] += fi;
     f[idx3] += fk;
     f[idx2] -= (fi+fk);
   }
   else if(abs(cos_theta) >0.1){
     R= crossProd(rji, rjk);
     normR = norm(R);
     inv_normR = 1.0/ normR;
     factor2 = -dV_dtheta*SIG(sin_theta)*(1.0/cos_theta);
     factor3 = inv_dist_rji* inv_dist_rjk;
     
     fi.x = factor2*( factor3*inv_normR*(-rjk.z*R.y + rjk.y*R.z)
                  -(factor3*(inv_dist_rji*inv_dist_rji)*normR)*rji.x );
     fi.y = factor2*( factor3*inv_normR*(rjk.z*R.x - rjk.x*R.z)
                  -(factor3*(inv_dist_rji*inv_dist_rji)*normR)*rji.y );
     fi.z = factor2*( factor3*inv_normR*(-rjk.y*R.x + rjk.x*R.y)
                  -(factor3*(inv_dist_rji*inv_dist_rji)*normR)*rji.z );
     
     fk.x = factor2*( factor3*inv_normR*(rji.z*R.y - rji.y*R.z)
                  -(factor3*(inv_dist_rjk*inv_dist_rjk)*normR)*rjk.x );
     fk.y = factor2*( factor3*inv_normR*(-rji.z*R.x + rji.x*R.z)
                  -(factor3*(inv_dist_rjk*inv_dist_rjk)*normR)*rjk.y );
     fk.z = factor2*( factor3*inv_normR*(rji.y*R.x - rji.x*R.y)
                  -(factor3*(inv_dist_rjk*inv_dist_rjk)*normR)*rjk.z );
     
     f[idx1] += fi;
     f[idx3] += fk;
     f[idx2] -= (fi+fk);
   }
   else{
     std::cerr<<" both sin and cos of theta are less than 0.1 "<<std::endl;
     exit(0);
   }
  }
}


    
/**============================================================
 * computeDihedrals: computes dihedral angle forces and dihedral angle energies for all
 * dihedral angles
 * (note this version is singularity free)
 */
void mdpack::md::computeDihedrals2(void){
  Edihedrals =0.0;

  for(int i=0; i< nDihedrals; i++){
    idx1 = dihedrals[i].atom1-1;
    idx2 = dihedrals[i].atom2-1;
    idx3 = dihedrals[i].atom3-1;
    idx4 = dihedrals[i].atom4-1;
    rij.x = atoms[idx1].r.x - atoms[idx2].r.x;
    rij.y = atoms[idx1].r.y - atoms[idx2].r.y;
    rij.z = atoms[idx1].r.z - atoms[idx2].r.z;
    rkj.x = atoms[idx3].r.x - atoms[idx2].r.x;
    rkj.y = atoms[idx3].r.y - atoms[idx2].r.y;
    rkj.z = atoms[idx3].r.z - atoms[idx2].r.z;
    rkl.x = atoms[idx3].r.x - atoms[idx4].r.x;
    rkl.y = atoms[idx3].r.y - atoms[idx4].r.y;
    rkl.z = atoms[idx3].r.z - atoms[idx4].r.z;
//     temp(0)=v1(1)*v2(2)-v1(2)*v2(1);
//   temp(1)=  v1(2)*v2(0)-v1(0)*v2(2);
//   temp(2)=v1(0)*v2(1)-v1(1)*v2(0);
    mm.x = rij.y*rkj.z - rij.z*rkj.y;
    mm.y = rij.z*rkj.x - rij.x*rkj.z;
    mm.z = rij.x*rkj.y - rij.y*rkj.x;
   // mm = crossProd( rij, rkj);
   // nn= crossProd(rkj, rkl);
    nn.x = rkj.y*rkl.z - rkj.z*rkl.y;
    nn.y = rkj.z*rkl.x - rkj.x*rkl.z;
    nn.z = rkj.x*rkl.y - rkj.y*rkl.x;
    
   // norm_rkj= norm(rkj);
    norm_rkj= sqrt(rkj.x*rkj.x + rkj.y*rkj.y +rkj.z*rkj.z);
    inv_norm_rkj =  1.0/norm_rkj;
    //inv_n = 1.0/dotProd(nn,nn);
    inv_n = 1.0/(nn.x*nn.x + nn.y*nn.y +nn.z*nn.z);
    //inv_m = 1.0/dotProd(mm,mm);
    inv_m = 1.0/(mm.x*mm.x + mm.y*mm.y +mm.z*mm.z);
    factor1= inv_norm_rkj*inv_norm_rkj;
    dot_ij_kj = rij.x*rkj.x + rij.y*rkj.y + rij.z*rkj.z;
    factor2 = factor1* dot_ij_kj;
    R.x= rij.x -( factor2*rkj.x );
    R.y= rij.y -( factor2*rkj.y );
    R.z= rij.z -( factor2*rkj.z );
    //S= rlk -( (pow2(inv_norm_rkj)*dotProd(rlk,rkj))*rkj );
    dot_kl_kj = rkl.x*rkj.x + rkl.y*rkj.y + rkl.z*rkj.z;
    factor3 = factor1*dot_kl_kj; 
    S.x= ( factor3*rkj.x ) -rkl.x;
    S.y= ( factor3*rkj.y ) -rkl.y;
    S.z= ( factor3*rkj.z ) -rkl.z;
    //inv_normR= inorm(R); // 1.0/norm(R);
    inv_normR = 1.0/sqrt(R.x*R.x + R.y*R.y+ R.z*R.z);
   // inv_normS= inorm(S); //1.0/norm(S);
    inv_normS = 1.0/sqrt(S.x*S.x + S.y*S.y+ S.z*S.z);
    factor4 = rij.x*nn.x + rij.y*nn.y +rij.z*nn.z;
    if( factor4>= 0.0)
      phi = acos( (inv_normR*inv_normS)*dotProd(R,S) );
    else if( factor4==0.0)
      phi=0.0;
    else
      phi = -1.0*acos( (inv_normR*inv_normS)*dotProd(R,S) );
  
    dihedrals[i].phi = phi;
    
    Edihedrals += dihedrals[i].k*( 1.0+ cos(dihedrals[i].n*phi 
                    +degTorad*dihedrals[i].delta) );
    KK = -dihedrals[i].n*dihedrals[i].k*
    sin(dihedrals[i].n*phi +degTorad*dihedrals[i].delta);
//      K = -md1.dihedrals[i].n*md1.dihedrals[i].k*
//     sin(md1.dihedrals[i].n*phi +degTorad*md1.dihedrals[i].delta);
    
    fi = -(KK*norm_rkj*inv_m)*mm;
    fl =  (KK*norm_rkj*inv_n)*nn;
    f[idx1] +=fi;
    f[idx4] +=fl;
    
   
    temp1 = inv_norm_rkj*inv_norm_rkj*dotProd(rij, rkj);
    temp2 = inv_norm_rkj*inv_norm_rkj*dotProd(rkl, rkj);
    
    //f[md1.dihedrals[i].atom2-1] += -fi + (temp1*fi) -(temp2*fl);
    f[idx2] -= fi;
    f[idx2] += (temp1*fi);
    f[idx2] -= (temp2*fl);
   
    //f[md1.dihedrals[i].atom3-1] += -fl - (temp1*fi) + (temp2*fl);
    f[idx3] -=fl;
    f[idx3] -= (temp1*fi);
    f[idx3] += (temp2*fl);
    
//     factor5 = KK*norm_rkj*inv_m;
//     
//     fi.x = -(factor5)*mm.x;
//     fi.y = -(factor5)*mm.y;
//     fi.z = -(factor5)*mm.z;
//     factor6 = KK*norm_rkj*inv_n;
//     fl.x =  (factor6)*nn.x;
//     fl.y =  (factor6)*nn.y;
//     fl.z =  (factor6)*nn.z;
//     f[idx1].x +=fi.x;
//     f[idx1].y +=fi.y;
//     f[idx1].z +=fi.z;
//     f[idx4].x +=fl.x;
//     f[idx4].y +=fl.y;
//     f[idx4].z +=fl.z;
//     
//     temp1 = inv_norm_rkj*inv_norm_rkj*dot_ij_kj;
//     temp2 = inv_norm_rkj*inv_norm_rkj*dot_kl_kj;
//     
//    
//     f[idx2].x += -fi.x + (temp1*fi.x) -(temp2*fl.x);
//     f[idx2].y += -fi.y + (temp1*fi.y) -(temp2*fl.y);
//     f[idx2].z += -fi.z + (temp1*fi.z) -(temp2*fl.z);
// //     
// //   
//     f[idx3].x += -fl.x - (temp1*fi.x) + (temp2*fl.x);
//     f[idx3].y += -fl.y - (temp1*fi.y) + (temp2*fl.y);
//     f[idx3].z += -fl.z - (temp1*fi.z) + (temp2*fl.z);
//    
//     f[dihedrals[i].atom2-1] -= fi;
//     f[dihedrals[i].atom2-1] += (temp1*fi);
//     f[dihedrals[i].atom2-1] -= (temp2*fl);
   
    //f[md1.dihedrals[i].atom3-1] += -fl - (temp1*fi) + (temp2*fl);
//     f[dihedrals[i].atom3-1] -=fl;
//     f[dihedrals[i].atom3-1] -= (temp1*fi);
//     f[dihedrals[i].atom3-1] += (temp2*fl);
    
  }
}
 /**===========================================================*/
 
void mdpack::md::computeDihedrals(void)  {
  Edihedrals =0.0;
  for (int i=0; i<nDihedrals; i++) {
    idx1 = dihedrals[i].atom1-1;
    idx2 = dihedrals[i].atom2-1;
    idx3 = dihedrals[i].atom3-1;
    idx4 = dihedrals[i].atom4-1;
//     const Vector *pos0 = coords + dihedral->atom1;
//     const Vector *pos1 = coords + dihedral->atom2;
//     const Vector *pos2 = coords + dihedral->atom3;
//     const Vector *pos3 = coords + dihedral->atom4;
    const mdpack::vector3 r12 = atoms[idx1].r - atoms[idx2].r;
    const mdpack::vector3 r23 = atoms[idx2].r - atoms[idx3].r;    //*pos1 - *pos2;
    const mdpack::vector3 r34 = atoms[idx3].r - atoms[idx4].r; //*pos2 - *pos3;

    mdpack::vector3 dcosdA;
    mdpack::vector3 dcosdB;
    mdpack::vector3 dsindC;
    mdpack::vector3 dsindB;
    mdpack:: vector3 f1, f2, f3;

    mdpack::vector3 A, B, C;
    A =crossProd(r12, r23);
    B =crossProd(r23, r34);
    C =crossProd(r23, A);

    double rA = norm(A); //A.length();
    double rB =  norm(B); //B.length(); 
    double rC = norm(C);  //C.length();

    double cos_phi = dotProd(A,B)/(rA*rB);  //(A*B)/(rA*rB);
    double sin_phi = dotProd(C,B)/(rC*rB);  //(C*B)/(rC*rB);

    // Normalize B
    rB = 1.0/rB;
    B *= rB;

     phi = -atan2(sin_phi, cos_phi);

    if (abs(sin_phi) > 0.1) {
      // Normalize A
      rA = 1.0/rA;
      A *= rA;
      dcosdA = rA*(cos_phi*A-B);
      dcosdB = rB*(cos_phi*B-A);
    }
    else {
      // Normalize C
      rC = 1.0/rC;
      C *= rC;
      dsindC = rC*(sin_phi*C-B);
      dsindB = rB*(sin_phi*B-C);
    }
 
    //int mult = dihedral->multiplicity; 
  //  for (int j=0; j<mult; j++) {
      double k = dihedrals[i].k;  //dihedral->k[j];
      double n = dihedrals[i].n;  //dihedral->n[j];
      double delta = dihedrals[i].delta;  //dihedral->delta[j];
      double K, K1;
      //if (n) {
        K = k * (1.0+cos(n*phi + degTorad*delta)); 
        K1 = -n*k*sin(n*phi + degTorad*delta);
      //}
//       else {
//         double diff = phi-delta;
//         if (diff < -M_PI) diff += 2.0*M_PI;
//         else if (diff > M_PI) diff -= 2.0*M_PI;
//         K = k*diff*diff;
//         K1 = 2.0*k*diff;
//       }
      Edihedrals += K;

      // forces
      if (abs(sin_phi) > 0.1) {
        K1 = K1/sin_phi;
        //f1.x += K1*(r23.y*dcosdA.z - r23.z*dcosdA.y);
	f1.x += K1*(r23.y*dcosdA.z - r23.z*dcosdA.y);
        //f1.y += K1*(r23.z*dcosdA.x - r23.x*dcosdA.z);
	f1.y += K1*(r23.z*dcosdA.x - r23.x*dcosdA.z);
        //f1.z += K1*(r23.x*dcosdA.y - r23.y*dcosdA.x);
	f1.z += K1*(r23.x*dcosdA.y - r23.y*dcosdA.x);

        //f3.x += K1*(r23.z*dcosdB.y - r23.y*dcosdB.z);
	f3.x += K1*(r23.z*dcosdB.y - r23.y*dcosdB.z);
        //f3.y += K1*(r23.x*dcosdB.z - r23.z*dcosdB.x);
	f3.y += K1*(r23.x*dcosdB.z - r23.z*dcosdB.x);
        //f3.z += K1*(r23.y*dcosdB.x - r23.x*dcosdB.y);
	f3.z += K1*(r23.y*dcosdB.x - r23.x*dcosdB.y);

//         f2.x += K1*(r12.z*dcosdA.y - r12.y*dcosdA.z
//                  + r34.y*dcosdB.z - r34.z*dcosdB.y);
	f2.x += K1*(r12.z*dcosdA.y - r12.y*dcosdA.z
                 + r34.y*dcosdB.z - r34.z*dcosdB.y);
//         f2.y += K1*(r12.x*dcosdA.z - r12.z*dcosdA.x
//                  + r34.z*dcosdB.x - r34.x*dcosdB.z);
	f2.y += K1*(r12.x*dcosdA.z - r12.z*dcosdA.x
                 + r34.z*dcosdB.x - r34.x*dcosdB.z);
//         f2.z += K1*(r12.y*dcosdA.x - r12.x*dcosdA.y
//                  + r34.x*dcosdB.y - r34.y*dcosdB.x);
	f2.z += K1*(r12.y*dcosdA.x - r12.x*dcosdA.y
                 + r34.x*dcosdB.y - r34.y*dcosdB.x);
      }
      else {
        //  This angle is closer to 0 or 180 than it is to
        //  90, so use the cos version to avoid 1/sin terms
        K1 = -K1/cos_phi;

//         f1.x += K1*((r23.y*r23.y + r23.z*r23.z)*dsindC.x
//                 - r23.x*r23.y*dsindC.y
//                 - r23.x*r23.z*dsindC.z);
	f1.x += K1*((r23.y*r23.y + r23.z*r23.z)*dsindC.x
                - r23.x*r23.y*dsindC.y
                - r23.x*r23.z*dsindC.z);
//         f1.y += K1*((r23.z*r23.z + r23.x*r23.x)*dsindC.y
//                 - r23.y*r23.z*dsindC.z
//                 - r23.y*r23.x*dsindC.x);
	f1.y += K1*((r23.z*r23.z + r23.x*r23.x)*dsindC.y
                - r23.y*r23.z*dsindC.z
                - r23.y*r23.x*dsindC.x);
//         f1.z += K1*((r23.x*r23.x + r23.y*r23.y)*dsindC.z
//                 - r23.z*r23.x*dsindC.x
//                 - r23.z*r23.y*dsindC.y);
	f1.z += K1*((r23.x*r23.x + r23.y*r23.y)*dsindC.z
                - r23.z*r23.x*dsindC.x
                - r23.z*r23.y*dsindC.y);

        f3 += K1*crossProd(dsindB,r23);

//         f2.x += K1*(-(r23.y*r12.y + r23.z*r12.z)*dsindC.x
//                +(2.0*r23.x*r12.y - r12.x*r23.y)*dsindC.y
//                +(2.0*r23.x*r12.z - r12.x*r23.z)*dsindC.z
//                +dsindB.z*r34.y - dsindB.y*r34.z);
	f2.x += K1*(-(r23.y*r12.y + r23.z*r12.z)*dsindC.x
               +(2.0*r23.x*r12.y - r12.x*r23.y)*dsindC.y
               +(2.0*r23.x*r12.z - r12.x*r23.z)*dsindC.z
               +dsindB.z*r34.y - dsindB.y*r34.z);
//         f2.y += K1*(-(r23.z*r12.z + r23.x*r12.x)*dsindC.y
//                +(2.0*r23.y*r12.z - r12.y*r23.z)*dsindC.z
//                +(2.0*r23.y*r12.x - r12.y*r23.x)*dsindC.x
//                +dsindB.x*r34.z - dsindB.z*r34.x);
	f2.y += K1*(-(r23.z*r12.z + r23.x*r12.x)*dsindC.y
               +(2.0*r23.y*r12.z - r12.y*r23.z)*dsindC.z
               +(2.0*r23.y*r12.x - r12.y*r23.x)*dsindC.x
               +dsindB.x*r34.z - dsindB.z*r34.x);
//         f2.z += K1*(-(r23.x*r12.x + r23.y*r12.y)*dsindC.z
//                +(2.0*r23.z*r12.x - r12.z*r23.x)*dsindC.x
//                +(2.0*r23.z*r12.y - r12.z*r23.y)*dsindC.y
//                +dsindB.y*r34.x - dsindB.x*r34.y);
	f2.z += K1*(-(r23.x*r12.x + r23.y*r12.y)*dsindC.z
               +(2.0*r23.z*r12.x - r12.z*r23.x)*dsindC.x
               +(2.0*r23.z*r12.y - r12.z*r23.y)*dsindC.y
               +dsindB.y*r34.x - dsindB.x*r34.y);
      }
   // }    // end loop over multiplicity
//     f[dihedral->atom1] += f1;
//     f[dihedral->atom2] += f2-f1;
//     f[dihedral->atom3] += f3-f2;
//     f[dihedral->atom4] += -f3;
    f[idx1] += f1;
    f[idx2] += f2-f1;
    f[idx3] += f3-f2;
    f[idx4] -= f3;

  }
 
}


void mdpack::md::computeImpDihedrals(void){
  EimpDihedrals =0.0;
  for (int i=0; i<nImpDihedrals; i++) {
//     const Vector *pos0 = coords + improper->atom1;
//     const Vector *pos1 = coords + improper->atom2;
//     const Vector *pos2 = coords + improper->atom3;
//     const Vector *pos3 = coords + improper->atom4;
//     const Vector r12 = *pos0 - *pos1;
//     const Vector r23 = *pos1 - *pos2;
//     const Vector r34 = *pos2 - *pos3;

    idx1 = impDihedrals[i].atom1-1;
    idx2 = impDihedrals[i].atom2-1;
    idx3 = impDihedrals[i].atom3-1;
    idx4 = impDihedrals[i].atom4-1;

    const mdpack::vector3 r12 = atoms[idx1].r - atoms[idx2].r;
    const mdpack::vector3 r23 = atoms[idx2].r - atoms[idx3].r;    //*pos1 - *pos2;
    const mdpack::vector3 r34 = atoms[idx3].r - atoms[idx4].r; //*pos2 - *pos3;
    
    mdpack::vector3 dcosdA;
    mdpack::vector3 dcosdB;
    mdpack::vector3 dsindC;
    mdpack::vector3 dsindB;
    mdpack::vector3 f1, f2, f3;

    mdpack::vector3 A, B, C;
    A = crossProd(r12, r23);
    B = crossProd(r23, r34);
    C = crossProd(r23, A);

    double rA = norm(A) ; //A.length();
    double rB = norm(B);  //B.length(); 
    double rC = norm(C);  //C.length();

    double cos_phi = dotProd(A,B)/(rA*rB); //(A*B)/(rA*rB);
    double sin_phi = dotProd(C,B)/(rC*rB); //(C*B)/(rC*rB);

    // Normalize B
    rB = 1.0/rB;
    B *= rB;

    phi = -atan2(sin_phi, cos_phi);

    if (abs(sin_phi) > 0.1) {
      // Normalize A
      rA = 1.0/rA;
      A *= rA;
      dcosdA = rA*(cos_phi*A-B);
      dcosdB = rB*(cos_phi*B-A);
    }
    else {
      // Normalize C
      rC = 1.0/rC;
      C *= rC;
      dsindC = rC*(sin_phi*C-B);
      dsindB = rB*(sin_phi*B-C);
    }
 
    //int mult = improper->multiplicity; 
   // for (int j=0; j<mult; j++) {
      double k = impDihedrals[i].k;  //improper->k[j];
     // double n = impDihedrals[i].n ;  //improper->n[j];
      double delta = impDihedrals[i].delta; //improper->delta[j];
      double K, K1;
//      // if (n) {
//         K = k * (1.0+cos(n*phi + delta)); 
//         K1 = -n*k*sin(n*phi + delta);
//      // }
      //else {
        double diff = phi-degTorad*delta;
        if (diff < -M_PI) diff += 2.0*M_PI;
        else if (diff > M_PI) diff -= 2.0*M_PI;
        K = k*diff*diff;
        K1 = 2.0*k*diff;
     // }
      EimpDihedrals += K;

      // forces
      if (abs(sin_phi) > 0.1) {
        K1 = K1/sin_phi;
//         f1.x += K1*(r23.y*dcosdA.z - r23.z*dcosdA.y);
//         f1.y += K1*(r23.z*dcosdA.x - r23.x*dcosdA.z);
//         f1.z += K1*(r23.x*dcosdA.y - r23.y*dcosdA.x);
// 
//         f3.x += K1*(r23.z*dcosdB.y - r23.y*dcosdB.z);
//         f3.y += K1*(r23.x*dcosdB.z - r23.z*dcosdB.x);
//         f3.z += K1*(r23.y*dcosdB.x - r23.x*dcosdB.y);
// 
//         f2.x += K1*(r12.z*dcosdA.y - r12.y*dcosdA.z
//                  + r34.y*dcosdB.z - r34.z*dcosdB.y);
//         f2.y += K1*(r12.x*dcosdA.z - r12.z*dcosdA.x
//                  + r34.z*dcosdB.x - r34.x*dcosdB.z);
//         f2.z += K1*(r12.y*dcosdA.x - r12.x*dcosdA.y
//                  + r34.x*dcosdB.y - r34.y*dcosdB.x);
	
	f1.x += K1*(r23.y*dcosdA.z - r23.z*dcosdA.y);
        f1.y += K1*(r23.z*dcosdA.x - r23.x*dcosdA.z);
        f1.z += K1*(r23.x*dcosdA.y - r23.y*dcosdA.x);

        f3.x += K1*(r23.z*dcosdB.y - r23.y*dcosdB.z);
        f3.y += K1*(r23.x*dcosdB.z - r23.z*dcosdB.x);
        f3.z += K1*(r23.y*dcosdB.x - r23.x*dcosdB.y);

        f2.x += K1*(r12.z*dcosdA.y - r12.y*dcosdA.z
                 + r34.y*dcosdB.z - r34.z*dcosdB.y);
        f2.y += K1*(r12.x*dcosdA.z - r12.z*dcosdA.x
                 + r34.z*dcosdB.x - r34.x*dcosdB.z);
        f2.z += K1*(r12.y*dcosdA.x - r12.x*dcosdA.y
                 + r34.x*dcosdB.y - r34.y*dcosdB.x);
	
	
      }
      else {
        //  This angle is closer to 0 or 180 than it is to
        //  90, so use the cos version to avoid 1/sin terms
        K1 = -K1/cos_phi;

//         f1.x += K1*((r23.y*r23.y + r23.z*r23.z)*dsindC.x
//                 - r23.x*r23.y*dsindC.y
//                 - r23.x*r23.z*dsindC.z);
//         f1.y += K1*((r23.z*r23.z + r23.x*r23.x)*dsindC.y
//                 - r23.y*r23.z*dsindC.z
//                 - r23.y*r23.x*dsindC.x);
//         f1.z += K1*((r23.x*r23.x + r23.y*r23.y)*dsindC.z
//                 - r23.z*r23.x*dsindC.x
//                 - r23.z*r23.y*dsindC.y);
	
	f1.x += K1*((r23.y*r23.y + r23.z*r23.z)*dsindC.x
                - r23.x*r23.y*dsindC.y
                - r23.x*r23.z*dsindC.z);
        f1.y += K1*((r23.z*r23.z + r23.x*r23.x)*dsindC.y
                - r23.y*r23.z*dsindC.z
                - r23.y*r23.x*dsindC.x);
        f1.z += K1*((r23.x*r23.x + r23.y*r23.y)*dsindC.z
                - r23.z*r23.x*dsindC.x
                - r23.z*r23.y*dsindC.y);

        f3 += K1*crossProd(dsindB,r23);

//         f2.x += K1*(-(r23.y*r12.y + r23.z*r12.z)*dsindC.x
//                +(2.0*r23.x*r12.y - r12.x*r23.y)*dsindC.y
//                +(2.0*r23.x*r12.z - r12.x*r23.z)*dsindC.z
//                +dsindB.z*r34.y - dsindB.y*r34.z);
//         f2.y += K1*(-(r23.z*r12.z + r23.x*r12.x)*dsindC.y
//                +(2.0*r23.y*r12.z - r12.y*r23.z)*dsindC.z
//                +(2.0*r23.y*r12.x - r12.y*r23.x)*dsindC.x
//                +dsindB.x*r34.z - dsindB.z*r34.x);
//         f2.z += K1*(-(r23.x*r12.x + r23.y*r12.y)*dsindC.z
//                +(2.0*r23.z*r12.x - r12.z*r23.x)*dsindC.x
//                +(2.0*r23.z*r12.y - r12.z*r23.y)*dsindC.y
//                +dsindB.y*r34.x - dsindB.x*r34.y);
	
	f2.x += K1*(-(r23.y*r12.y + r23.z*r12.z)*dsindC.x
               +(2.0*r23.x*r12.y - r12.x*r23.y)*dsindC.y
               +(2.0*r23.x*r12.z - r12.x*r23.z)*dsindC.z
               +dsindB.z*r34.y - dsindB.y*r34.z);
        f2.y += K1*(-(r23.z*r12.z + r23.x*r12.x)*dsindC.y
               +(2.0*r23.y*r12.z - r12.y*r23.z)*dsindC.z
               +(2.0*r23.y*r12.x - r12.y*r23.x)*dsindC.x
               +dsindB.x*r34.z - dsindB.z*r34.x);
        f2.z += K1*(-(r23.x*r12.x + r23.y*r12.y)*dsindC.z
               +(2.0*r23.z*r12.x - r12.z*r23.x)*dsindC.x
               +(2.0*r23.z*r12.y - r12.z*r23.y)*dsindC.y
               +dsindB.y*r34.x - dsindB.x*r34.y);
      }
    //}    // end loop over multiplicity
//     f[improper->atom1] += f1;
//     f[improper->atom2] += f2-f1;
//     f[improper->atom3] += f3-f2;
//     f[improper->atom4] += -f3;
    f[idx1] += f1;
    f[idx2] += f2-f1;
    f[idx3] += f3-f2;
    f[idx4] -= f3;

    
  }
 
}



 
 
/**============================================================
 * \c computeImpDihedrals: computes improper dihedral angle forces and  energies for all
 * improper dihedrals
 * (note this version is singularity free)
 */ 
void mdpack::md::computeImpDihedrals2(void){
   EimpDihedrals =0.0;

  for(int i=0; i< nImpDihedrals; i++){
    idx1 = impDihedrals[i].atom1-1;
    idx2 = impDihedrals[i].atom2-1;
    idx3 = impDihedrals[i].atom3-1;
    idx4 = impDihedrals[i].atom4-1;
    rij.x = atoms[idx1].r.x - atoms[idx2].r.x;
    rij.y = atoms[idx1].r.y - atoms[idx2].r.y;
    rij.z = atoms[idx1].r.z - atoms[idx2].r.z;
    rkj.x = atoms[idx3].r.x - atoms[idx2].r.x;
    rkj.y = atoms[idx3].r.y - atoms[idx2].r.y;
    rkj.z = atoms[idx3].r.z - atoms[idx2].r.z;
    rkl.x = atoms[idx3].r.x - atoms[idx4].r.x;
    rkl.y = atoms[idx3].r.y - atoms[idx4].r.y;
    rkl.z = atoms[idx3].r.z - atoms[idx4].r.z;
   
    mm = crossProd( rij, rkj);
    nn= crossProd(rkj, rkl);
    
    norm_rkj= norm(rkj);
    inv_norm_rkj = 1.0/norm_rkj;
    inv_n = 1.0/dotProd(nn,nn);
    inv_m = 1.0/dotProd(mm,mm);
    R= rij -( (pow2(inv_norm_rkj)*dotProd(rij,rkj))*rkj );
    //S= rlk -( (pow2(inv_norm_rkj)*dotProd(rlk,rkj))*rkj );
    S= ( (pow2(inv_norm_rkj)*dotProd(rkl,rkj))*rkj ) -rkl;
    normR= norm(R);
    normS= norm(S);
    if( dotProd(rij, nn)> 0.0)
      phi = acos( (1.0/(normR*normS))*dotProd(R,S) );
    else if (dotProd(rij, nn)==0.0)
      phi=0.0;
    else
      phi = -1.0*acos( (1.0/(normR*normS))*dotProd(R,S) );
    impDihedrals[i].phi = phi;
    
    KK = 2.0*impDihedrals[i].k*(phi -degTorad*impDihedrals[i].delta);
    EimpDihedrals += impDihedrals[i].k*pow2(phi -degTorad*impDihedrals[i].delta);
    
    fi.x = -(KK*norm_rkj*inv_m)*mm.x;
    fi.y = -(KK*norm_rkj*inv_m)*mm.y;
    fi.z = -(KK*norm_rkj*inv_m)*mm.z;
    fl.x =  (KK*norm_rkj*inv_n)*nn.x;
    fl.y =  (KK*norm_rkj*inv_n)*nn.y;
    fl.z =  (KK*norm_rkj*inv_n)*nn.z;
    f[idx1].x +=fi.x;
    f[idx1].y +=fi.y;
    f[idx1].z +=fi.z;
    f[idx4].x +=fl.x;
    f[idx4].y +=fl.y;
    f[idx4].z +=fl.z;
    
    factor1 = 1.0/(norm_rkj*norm_rkj);
    temp1 = factor1*dotProd(rij, rkj);
    temp2 = factor1*dotProd(rkl, rkj);
    
    
    
    f[idx2].x += -fi.x + (temp1*fi.x) -(temp2*fl.x);
    f[idx2].y += -fi.y + (temp1*fi.y) -(temp2*fl.y);
    f[idx2].z += -fi.z + (temp1*fi.z) -(temp2*fl.z);
//     
//   
    f[idx3].x += -fl.x - (temp1*fi.x) + (temp2*fl.x);
    f[idx3].y += -fl.y - (temp1*fi.y) + (temp2*fl.y);
    f[idx3].z += -fl.z - (temp1*fi.z) + (temp2*fl.z);
   
    
//     f[idx2] -= fi;
//     f[idx2] += (temp1*fi);
//     f[idx2] -= (temp2*fl);
//    
//     f[idx3] -=fl;
//     f[idx3] -= (temp1*fi);
//     f[idx3] += (temp2*fl);
    
  }
}
/**==========================================================*/  
      
/**============================================================
 * \c computeSphBC: computes improper dihedral angle forces and  energies for all
 * improper dihedrals
 * (note this version is singularity free)
 */   
void mdpack::md::computeSphBC(void){
 
  EsphBC=0.0;
  if(sphBC)
    for(int i=0; i<nAtoms ; i++){
      rr.x= sphBCCenter.x - atoms[i].r.x;
      rr.y= sphBCCenter.y - atoms[i].r.y;
      rr.z= sphBCCenter.z - atoms[i].r.z;
      mdpack::dreal dist2 = ( rr.x*rr.x + rr.y*rr.y+ rr.z*rr.z);
      if(dist2 > sphBCr1*sphBCr1){
	dist = sqrt(dist2);
	factor1=dist-sphBCr1;
	factor2 = sphBCk1*sphBCexp1/dist;
	factor3 = pow(factor1,sphBCexp1-1);
	f[i].x +=factor2*factor3*rr.x;
	f[i].y +=factor2*factor3*rr.y;
	f[i].z +=factor2*factor3*rr.z;
	EsphBC += sphBCk1*factor3*factor1;
      }
    }
}
/**===========================================================*/

/**============================================================
 *  computes forces
 */
void mdpack::md::computeForces(void){
  for (int i=0; i< nAtoms; i++)
    f[i].clear();
  EnonBondeds =0.0;
  Ebonds =0.0;
  Eangles =0.0;
  Edihedrals=0.0;
  EimpDihedrals =0.0;
  EsphBC =0.0;
  Evdw =0.0;
  Eelec=0.0;
  
  computeNonBondeds();
  computeBonds();
  computeAngles();
  computeSphBC();
  computeDihedrals();
  computeImpDihedrals();
  
  EnonBondeds= Evdw + Eelec;
  Epotential = EnonBondeds + Ebonds + Eangles +
                Edihedrals + EimpDihedrals + EsphBC;
}
/**===========================================================*/


/**============================================================
 *  verlet integrator 
 */
void mdpack::md::verlet(void){
  for(int i=0; i<nAtoms ; i++){
    atoms[i].v.x += 0.5*h*f[i].x*atoms[i].im;
    atoms[i].v.y += 0.5*h*f[i].y*atoms[i].im;
    atoms[i].v.z += 0.5*h*f[i].z*atoms[i].im;
    atoms[i].r.x += h* atoms[i].v.x;
    atoms[i].r.y += h* atoms[i].v.y;
    atoms[i].r.z += h* atoms[i].v.z;
  }
 // moveToCM(md1);

  Ekinetic =0.0;
  Etotal =0.0;
  
  computeForces();
  for (int i=0; i< nAtoms ; i++){
    atoms[i].v.x += 0.5*h*f[i].x*atoms[i].im;
    atoms[i].v.y += 0.5*h*f[i].y*atoms[i].im;
    atoms[i].v.z += 0.5*h*f[i].z*atoms[i].im;
    
    Ekinetic +=0.5*(atoms[i].m)*
		   (atoms[i].v.x*atoms[i].v.x +
		    atoms[i].v.y*atoms[i].v.y +
		    atoms[i].v.z*atoms[i].v.z);
  }
  if(useCellLists){
    for(int i=0; i<nCells ; i++)
      cells[i].update(atoms);
    for(int i=0; i<nCells; i++)
      cells[i].update(cells, atoms);
  }
  Etotal = Ekinetic + Epotential;
  tempInstant = 2.0*INV_BOLTZMANN*(1.0/(nDegsFreedom-4))*Ekinetic;
  //tempAvg += tempInstant;
  tempAvg = (stepn*tempAvg +tempInstant)/(stepn +1.0);
  simTime +=timeStep;
  stepn +=1;
}
/**===========================================================*/

/**===========================================================
 * Langevin method: based on spliting techniques 
 * ( weak second order and first order strong)
 * see Emad Noorizadeh thesis 2010
 * 
 * Performs constant temperature simulation
 */
void mdpack::md::langevin(void){
 
  for(int i=0;i< nAtoms; i++){
    atoms[i].v.x = expGamma_h*atoms[i].v.x + 0.5*atoms[i].im*h*f[i].x + langFactor[i]*W();
    atoms[i].v.y = expGamma_h*atoms[i].v.y + 0.5*atoms[i].im*h*f[i].y + langFactor[i]*W();
    atoms[i].v.z = expGamma_h*atoms[i].v.z + 0.5*atoms[i].im*h*f[i].z + langFactor[i]*W();
    
    atoms[i].r.x +=  h*atoms[i].v.x;
    atoms[i].r.y +=  h*atoms[i].v.y;
    atoms[i].r.z +=  h*atoms[i].v.z;
  }
  
  Ekinetic =0.0;
  Etotal =0.0;
  
  computeForces();
  for (int i=0; i< nAtoms ; i++){
    atoms[i].v.x += 0.5*h*f[i].x*atoms[i].im;
    atoms[i].v.y += 0.5*h*f[i].y*atoms[i].im;
    atoms[i].v.z += 0.5*h*f[i].z*atoms[i].im;
    
    Ekinetic +=0.5*(atoms[i].m)*
		   (atoms[i].v.x*atoms[i].v.x +
		    atoms[i].v.y*atoms[i].v.y +
		    atoms[i].v.z*atoms[i].v.z);
  }
  if(useCellLists){
    for(int i=0; i<nCells ; i++)
      cells[i].update(atoms);
    for(int i=0; i<nCells; i++)
      cells[i].update(cells, atoms);
  }
  Etotal = Ekinetic + Epotential;
  tempInstant = 2.0*INV_BOLTZMANN*(1.0/nDegsFreedom)*Ekinetic;
  //tempAvg += tempInstant;
  tempAvg = (stepn*tempAvg +tempInstant)/(stepn +1.0);
  simTime +=timeStep;
  stepn +=1;
}
/**=====================================================*/
 

/**===========================================================
 * Langevin method: based on spliting techniques 
 * see Emad Noorizadeh thesis 2010
 * 
 * Performs constant temperature simulation
 */
void mdpack::md::langevin1(void){
 
  for(int i=0;i< nAtoms; i++){
    atoms[i].v.x = expGamma_hh*atoms[i].v.x + 0.5*atoms[i].im*h*f[i].x + langFactorh[i]*W();
    atoms[i].v.y = expGamma_hh*atoms[i].v.y + 0.5*atoms[i].im*h*f[i].y + langFactorh[i]*W();
    atoms[i].v.z = expGamma_hh*atoms[i].v.z + 0.5*atoms[i].im*h*f[i].z + langFactorh[i]*W();
    
    atoms[i].r.x +=  h*atoms[i].v.x;
    atoms[i].r.y +=  h*atoms[i].v.y;
    atoms[i].r.z +=  h*atoms[i].v.z;
  }
  
  Ekinetic =0.0;
  Etotal =0.0;
  
  computeForces();
  for (int i=0; i< nAtoms ; i++){
    atoms[i].v.x = expGamma_hh*atoms[i].v.x + 0.5*atoms[i].im*h*f[i].x + langFactorh[i]*W();
    atoms[i].v.y = expGamma_hh*atoms[i].v.y + 0.5*atoms[i].im*h*f[i].y + langFactorh[i]*W();
    atoms[i].v.z = expGamma_hh*atoms[i].v.z + 0.5*atoms[i].im*h*f[i].z + langFactorh[i]*W();
    
    
    Ekinetic +=0.5*(atoms[i].m)*
		   (atoms[i].v.x*atoms[i].v.x +
		    atoms[i].v.y*atoms[i].v.y +
		    atoms[i].v.z*atoms[i].v.z);
  }
  if(useCellLists){
    for(int i=0; i<nCells ; i++)
      cells[i].update(atoms);
    for(int i=0; i<nCells; i++)
      cells[i].update(cells, atoms);
  }
  Etotal = Ekinetic + Epotential;
  tempInstant = 2.0*INV_BOLTZMANN*(1.0/nDegsFreedom)*Ekinetic;
  tempAvg = (stepn*tempAvg +tempInstant)/(stepn +1.0);
  simTime +=timeStep;
  stepn +=1;
}
/**=====================================================*/

/**===========================================================
 * Nose-Hoover-Langevin method: based on spliting techniques 
 * see Emad Noorizadeh thesis 2010
 * 
 * Performs constant temperature simulation
 */

void mdpack::md::noseHooverLangevin(void){

  Ekinetic =0.0;
  for(int i=0;i< nAtoms; i++){
    atoms[i].v.x =(atoms[i].v.x + 0.5*h*atoms[i].im*f[i].x)/(1.0+ 0.5*h*nhl_xi);
    atoms[i].v.y =(atoms[i].v.y + 0.5*h*atoms[i].im*f[i].y)/(1.0+ 0.5*h*nhl_xi);
    atoms[i].v.z =(atoms[i].v.z + 0.5*h*atoms[i].im*f[i].z)/(1.0+ 0.5*h*nhl_xi);
    
    atoms[i].r.x +=  h*atoms[i].v.x; 
    atoms[i].r.y +=  h*atoms[i].v.y;  
    atoms[i].r.z +=  h*atoms[i].v.z; 
    
    Ekinetic +=(atoms[i].m)*(atoms[i].v.x*atoms[i].v.x +
		             atoms[i].v.y*atoms[i].v.y +
		             atoms[i].v.z*atoms[i].v.z); 
  }
    
  nhl_xi= expNhl*nhl_xi + (nhl_imu*h)*(Ekinetic - (nDegsFreedom+1)*BOLTZMANN*targetTemp)
            +nhlFactor*W();
	    
  Ekinetic =0.0;
  Etotal =0.0;
  computeForces();
  for(int i=0;i<nAtoms; i++){
      atoms[i].v.x += 0.5*h*atoms[i].im*f[i].x - 0.5*h*nhl_xi*atoms[i].v.x; 
      atoms[i].v.y += 0.5*h*atoms[i].im*f[i].y - 0.5*h*nhl_xi*atoms[i].v.y; 
      atoms[i].v.z += 0.5*h*atoms[i].im*f[i].z - 0.5*h*nhl_xi*atoms[i].v.z; 
	
      Ekinetic +=0.5*(atoms[i].m)*(atoms[i].v.x*atoms[i].v.x +
		             atoms[i].v.y*atoms[i].v.y +
		             atoms[i].v.z*atoms[i].v.z); 
        
  }
  Etotal = Ekinetic + Epotential;
  tempInstant = 2.0*INV_BOLTZMANN*(1.0/(nDegsFreedom))*Ekinetic;
  tempAvg = (stepn*tempAvg +tempInstant)/(stepn +1.0);
  simTime +=timeStep;
  stepn +=1;
 
}
/**====================================================================*/

/**=====================================================================
 * over damped dynamics
 *  
 */

void mdpack::md::overDamped(void){
  for(int i=0;i<nAtoms; i++){
    atoms[i].r.x +=  h*f[i].x + overDampedFactor*W();
    atoms[i].r.y +=  h*f[i].y + overDampedFactor*W();
    atoms[i].r.z +=  h*f[i].z + overDampedFactor*W();
  }
  Ekinetic =0.0;
  Etotal =0.0;
  computeForces();
  Etotal = Ekinetic + Epotential;
  simTime +=timeStep;
  stepn +=1;
}
/**================================================================*/

/**=====================================================================
 * Minimization using gradient descent
 *  
 */
void mdpack::md::minimize(void){
  for(int i=0; i<nAtoms ; i++){
    atoms[i].v = 0.5*atoms[i].im*h*f[i];
    atoms[i].r += h*atoms[i].v;
  }
  Etotal =0.0;
  computeForces();
  for (int i=0; i< nAtoms ; i++)
    atoms[i].v.clear();
  Etotal =  Epotential;
}

/**=====================================================================*/
  
/**====================================================================
 * Moves atoms to the center of mass
 */
void mdpack::md::moveToCM(){
  mdpack::vector3 temp_cm;
  for(int i=0; i< nAtoms; i++)
    temp_cm += atoms[i].r;
  temp_cm *= (1.0/centerOfMass_m);
  for(int i=0; i<nAtoms; i++)
    atoms[i].r +=(centerOfMassCoords- temp_cm);
}
/**=================================================================*/

/**============================================================
 * \c openOutputFiles: opens necessary file for writing data
 * 
 */ 

void mdpack::md::openOutputFiles(void){
  char tempName [50];
  strcpy(tempName, outputName);
  energyLogFile.open(strcat(tempName,".log"), std::ofstream::out | std::ofstream::trunc );
	
  energyLogFile <<"ETITLE:"<<ESPACE<< "TS" <<ESPACE<< "BOND" <<ESPACE<<"ANGLE"<<ESPACE
                <<"DIHED"<<ESPACE<<"IMPRP"<<ESPACE<< "ELECT" <<ESPACE<<"VDW"<<ESPACE
                <<"BOUNDARY"<<ESPACE<<"MISC" <<ESPACE<<"KINETIC" <<ESPACE<<"TOTAL" << ESPACE
                <<"TEMP" <<ESPACE<<"POTENTIAL"<<ESPACE<<"TEMPAVG"<<"\n\n";
  if(isWritePosOn){
    strcpy(tempName, outputName);
    if(isBinaryOn)
      posFile.open(strcat(tempName,"Pos.xyz"),  
			 std::ofstream::out | std::ofstream::trunc | std::ios::binary);
    else
      posFile.open(strcat(tempName,"Pos.xyz"),  
			 std::ofstream::out | std::ofstream::trunc);
  }
	      
  if(isWriteVelOn){
    strcpy(tempName, outputName);
    if(isBinaryOn)
      velFile.open(strcat(tempName,"Vel.xyz"), 
			 std::ofstream::out | std::ofstream::trunc | std::ofstream::binary );
    else
      velFile.open(strcat(tempName,"Vel.xyz"), 
			 std::ofstream::out | std::ofstream::trunc);
  }
  if(isWritePhiOn){
    strcpy(tempName, outputName);
    if(isBinaryOn)
      phiFile.open(strcat(tempName,"Phi.dat"), 
			 std::ofstream::out | std::ofstream::trunc | std::ofstream::binary );
    else
      phiFile.open(strcat(tempName,"Phi.dat"), 
			 std::ofstream::out | std::ofstream::trunc | std::ofstream::binary );
  }
	    
}
/**=============================================================*/


/**============================================================
 * \c writeEnergie: writes energies to the output file named output.log
 * 
 */ 
/*
ETITLE:      TS           BOND          ANGLE          DIHED          IMPRP             
             ELECT        VDW       BOUNDARY         MISC        KINETIC          
             TOTAL         TEMP      POTENTIAL      TEMPAVG
             */

void mdpack::md::writeEnergies(void){
  energyLogFile.precision(4);
  energyLogFile <<"ETITLE:"<< std::fixed<< ESPACE<< stepn <<ESPACE<< Ebonds <<ESPACE<<Eangles<<ESPACE
                <<Edihedrals<<ESPACE<<EimpDihedrals<<ESPACE<< Eelec <<ESPACE<<Evdw<<ESPACE
                <<EsphBC<<ESPACE<< Emisc <<ESPACE<< Ekinetic<<ESPACE << Etotal << ESPACE
                <<tempInstant <<ESPACE<<Epotential<<ESPACE<< tempAvg<<"\n";
}
/**=============================================================*/

/**============================================================
 * \c printEnergies: prints energies to the std output 
 * 
 */ 
/*
ETITLE:      TS           BOND          ANGLE          DIHED          IMPRP             
             ELECT        VDW       BOUNDARY         MISC        KINETIC          
             TOTAL         TEMP      POTENTIAL      TEMPAVG
             */

void mdpack::md::printEtitles(void){
 std::cout<<"ETITLE:"<<ESPACE<< "TS" <<ESPACE<< "BOND" <<ESPACE<<"ANGLE"<<ESPACE
                <<"DIHED"<<ESPACE<<"IMPRP"<<ESPACE<< "ELECT" <<ESPACE<<"VDW"<<ESPACE
                <<"BOUNDARY"<<ESPACE<<"MISC" <<ESPACE<<"KINETIC" <<ESPACE<<"TOTAL" << ESPACE
                <<"TEMP" <<ESPACE<<"POTENTIAL"<<ESPACE<<"TEMPAVG"<<"\n";
 std::cout<<std::endl;
}
		
void mdpack::md::printEnergies(void){
  printf("ETITLE: %12d  %12.4f  %12.4f  %12.4f  %12.4f  %12.4f  %12.4f  %12.4f  %12.4f  %12.4f  ",
	 stepn,Ebonds,Eangles,Edihedrals,EimpDihedrals, Eelec, Evdw ,EsphBC, Emisc, Ekinetic);
  printf("%12.4f  %12.4f  %12.4f  %12.4f \n\n", Etotal, tempInstant, Epotential, tempAvg);
	 
//   std::cout <<"ETITLE:"<<ESPACE<< stepn <<ESPACE<< Ebonds <<ESPACE<<Eangles<<ESPACE
//                 <<Edihedrals<<ESPACE<<EimpDihedrals<<ESPACE<< Eelec <<ESPACE<<Evdw<<ESPACE
//                 <<EsphBC<<ESPACE<< Emisc <<ESPACE<< Ekinetic<<ESPACE << Etotal << ESPACE
//                 <<tempInstant <<ESPACE<<Epotential<<ESPACE<< tempAvg<<std::endl;
}
/**=============================================================*/


/**============================================================
 * \c printInfo: prints simulation details to the std output 
 * 
 */ 
void mdpack::md::printInfo(void){
  //------------------------------------------------------------
  std::cout<<"Info: ********************************************************************"<<"\n";
  std::cout<<std::left<<"Info:               SIMULATION PARAMETERS               "<<"\n";
  std::cout<<"Info: ********************************************************************"<<"\n";
  std::cout.width(40); std::cout<<"Info: TIMESTEP  "; std::cout.width(5); std::cout<<std::left<< timeStep<<"\n";
  std::cout.width(40); std::cout<<"Info: NUMBER OF STEPS  "; std::cout.width(5); 
  std::cout<< numSteps<<"\n";             
  std::cout.width(40); std::cout<<"Info: TEMPERATURE  "; std::cout.width(5); std::cout<<targetTemp<<"\n";
  std::cout.width(40); std::cout<<"Info: DIELECTRIC  "; std::cout.width(5); std::cout<<e0<<"\n";
  std::cout.width(40); std::cout<<"Info: EXCLUDE"; std::cout.width(5); std::cout<<exclude<<"\n";
  std::cout.width(40); std::cout<<"Info: OUTPUT FILENAME"; std::cout.width(5); 
  std::cout<<outputName<<"\n";
  std::cout.width(40); std::cout<<"Info: CUTOFF"; std::cout.width(5); std::cout<<cutoff<<"\n";
  std::cout.width(40); std::cout<<"Info: SWITCHDISTANCE"; std::cout.width(5); std::cout<<switchdist<<"\n";
  if(sphBC){
    std::cout.width(40); std::cout<<"Info: SPHERICAL BOUNDARY CONDITIONS"; std::cout.width(5); 
    std::cout<<"ACTIVE"<<"\n";
    std::cout.width(40); std::cout<<"Info: RADIUS #1 "; std::cout.width(5); std::cout<<sphBCr1<<"\n";
    std::cout.width(40); std::cout<<"Info: FORCE CONSTANT #1 "; std::cout.width(5); 
    std::cout<<sphBCk1<<"\n";
    std::cout.width(40); std::cout<<"Info: EXPONENT #1 "; std::cout.width(5); std::cout<<sphBCexp1<<"\n";
    std::cout.width(40); std::cout<<"Info: SPHERE BOUNDARY CENTER    "; 
    std::cout<<sphBCCenter.x<<"  "<<sphBCCenter.y<<"  "<<sphBCCenter.z<<"\n";
  }
  //-------------------------------------------------------------------
  std::cout<<"Info: *********************************************************************"<<"\n";
  std::cout<<"Info:               STRUCTURE SUMMARY   "<<"\n";
  std::cout<<"Info: *********************************************************************"<<"\n";
  std::cout<<"Info: "<<nAtoms<<" ATOMS"<<"\n";
  std::cout<<"Info: "<<nBonds<<" BONDS"<<"\n";
  std::cout<<"Info: "<<nAngles<<" ANGLES"<<"\n";
  std::cout<<"Info: "<<nDihedrals<<" DIHEDRALS"<<"\n";
  std::cout<<"Info: "<<nImpDihedrals<<" IMPROPERS"<<"\n";
  std::cout<<"Info: *********************************************************************"<<"\n";

}

/**================================================================*/

  
/**============================================================
 * \c writePositions: write positions to the output file named outputNamePos.xyz
 * The formatting of the .xyz file format is as follows:
 * <number of atoms>
 * comment line
 * atom_symbol1 x-coord1 y-coord1 z-coord1
 * atom_symbol2 x-coord2 y-coord2 z-coord2
 * ...
 * atom_symboln x-coordn y-coordn z-coordn
 * 
 * Example:
 * 5
 * methane molecule (in [[Ångström]]s)
 * C        0.000000        0.000000        0.000000
 * H        0.000000        0.000000        1.089000
 * H        1.026719        0.000000       -0.363000
 * H       -0.513360       -0.889165       -0.363000
 * H       -0.513360        0.889165       -0.363000
 * 
 * Note: if binary format is on then all position written in a line as follows
 * 
 * atom1.x_iatom1.y_iatom.z_iatom2.x_i.......atom1.x_i+1atom1.y_i+1
 * 
 * where i is the ith step, so all data will be written consequently and without space
 * the data whitten as type double so knowing thin and the number of atoms 
 * this file can be processed.
 */     
void mdpack::md::writePositions(void){
  if(isBinaryOn){
    double var;
    for(int i=0; i<nAtoms; i++){
      var= atoms[i].r.x;
      posFile.write((char *) &var, sizeof var);
      var= atoms[i].r.y;
      posFile.write((char *) &var, sizeof var);
      var= atoms[i].r.z;
      posFile.write((char *) &var, sizeof var);
    }
  }
  else{
    posFile << nAtoms <<"\n";
    posFile << simName <<"\n";
    posFile.precision(5);
    for(int i=0; i<nAtoms; i++){
      posFile <<atoms[i].name[0]<<"  "<<std::fixed<<atoms[i].r.y<<"  "<<atoms[i].r.y<<"  "
            <<atoms[i].r.z<<"\n";
    }
  }
}

/**==========================================================*/
/**============================================================
 * \c writeVelocities: writes velocities to the output file named outputNameVel.xyz
 * The formatting of the .xyz (see above method writePositions for description)
 * 
 */

void mdpack::md::writeVelocities(void){
  if(isBinaryOn){
    double var;
    for(int i=0; i<nAtoms; i++){
      var= atoms[i].v.x;
      velFile.write((char *) &var, sizeof var);
      var= atoms[i].v.y;
      velFile.write((char *) &var, sizeof var);
      var= atoms[i].v.z;
      velFile.write((char *) &var, sizeof var);
    }
  }
  else{
    posFile << nAtoms <<"\n";
    posFile << simName <<"\n";
    posFile.precision(5);
    for(int i=0; i<nAtoms; i++){
      velFile <<atoms[i].name[0]<<"  "<<std::fixed<<atoms[i].v.y<<"  "<<atoms[i].v.y<<"  "
            <<atoms[i].v.z<<"\n";
    }
  }
}
/**==============================================================*/

/**============================================================
 * \c writeRestart: writes two PDB files one containing all
 * positions and one containing all velocities. These file are used to 
 * restart the simulation at a given state and a time
 * step_i: is the time step
 * t: is the time
 */

void mdpack::md::writeRestart(const int& step_i,const mdpack::dreal& t){
  char tempName [50];
  strcpy(tempName, outputName);
  restartPosFile.open(strcat(tempName,"Restart.coor"), std::ofstream::out | std::ofstream::trunc );
  strcpy(tempName, outputName);
  restartVelFile.open(strcat(tempName,"Restart.vel"), std::ofstream::out | std::ofstream::trunc );
  
  restartPosFile << "REMARKS  positions at time step"
                 <<"("<<step_i<<")"<<" and time  "<<"("<<t<<")"<<"\n";
  restartVelFile << "REMARKS  velocities at time step"
                 <<"("<<step_i<<")"<<" and time  "<<"("<<t<<")"<<"\n";
		 
 
  for(int i=0; i<nAtoms; i++){
    restartPosFile<< "ATOM ";
    restartPosFile.width(7);
    restartPosFile <<std::right<<atoms[i].id<<"   ";
    restartPosFile.width(5);
    restartPosFile<<atoms[i].name<<"    ";
    restartPosFile.width(5);
    restartPosFile<<atoms[i].resName<<"   ";
    restartPosFile.width(4);
    restartPosFile <<atoms[i].resId<<"     ";
    restartPosFile.width(8);
    restartPosFile.precision(4);
    restartPosFile <<std::fixed<<atoms[i].r.x<<"  ";
    restartPosFile.width(8);
    restartPosFile.precision(4);
    restartPosFile <<std::right<<std::fixed<<atoms[i].r.y<<"  ";
    restartPosFile.width(8);
    restartPosFile.precision(4);
    restartPosFile <<std::right<<std::fixed<<atoms[i].r.z<<"  ";	   
    restartPosFile.precision(2);
    restartPosFile <<std::fixed<<atoms[i].occupancy<<"   ";
    restartPosFile <<std::fixed<<atoms[i].tempFactor<<"      ";
    restartPosFile.width(5);		   
    restartPosFile <<std::left<<atoms[i].segName<<"   ";
    restartPosFile.width(1);
    restartPosFile<<std::right<<atoms[i].name[0]<<"\n";
    
    restartVelFile<< "ATOM ";
    restartVelFile.width(7);
    restartVelFile <<std::right<<atoms[i].id<<"   ";
    restartVelFile.width(5);
    restartVelFile<<atoms[i].name<<"    ";
    restartVelFile.width(5);
    restartVelFile<<atoms[i].resName<<"   ";
    restartVelFile.width(4);
    restartVelFile <<atoms[i].resId<<"     ";
    restartVelFile.width(8);
    restartVelFile.precision(4);
    restartVelFile <<std::fixed<<atoms[i].v.x<<"  ";
    restartVelFile.width(8);
    restartVelFile.precision(4);
    restartVelFile <<std::right<<std::fixed<<atoms[i].v.y<<"  ";
    restartVelFile.width(8);
    restartVelFile.precision(4);
    restartVelFile <<std::right<<std::fixed<<atoms[i].v.z<<"  ";	   
    restartVelFile.precision(2);
    restartVelFile <<std::fixed<<atoms[i].occupancy<<"   ";
    restartVelFile <<std::fixed<<atoms[i].tempFactor<<"      ";
    restartVelFile.width(5);		   
    restartVelFile <<std::left<<atoms[i].segName<<"   ";
    restartVelFile.width(1);
    restartVelFile<<std::right<<atoms[i].name[0]<<"\n";
  }
}
/**=============================================================*/
  
/**============================================================
 * \c writepdb: writes pdb file
 */  
  
void mdpack::md::writepdb(const int& step_n, const mdpack::dreal& t ){
  char tempName [50];
  char step [20];
  sprintf(step, "%d", step_n);

  //itoa(step_n,step,10);
  strcpy(tempName, outputName);
  strcat(tempName,step);
  outpdbFile.open(strcat(tempName,".pdb"), std::ofstream::out | std::ofstream::trunc );
  
  outpdbFile << "REMARKS  positions at time step"
                 <<"("<<step_n<<")"<<" and time  "<<"("<<t<<")"<<"\n";
  for(int i=0; i<nAtoms; i++){
    outpdbFile<< "ATOM ";
    outpdbFile.width(7);
    outpdbFile <<std::right<<atoms[i].id<<"   ";
    outpdbFile.width(5);
    outpdbFile<<atoms[i].name<<"    ";
    outpdbFile.width(5);
    outpdbFile<<atoms[i].resName<<"   ";
    outpdbFile.width(4);
    outpdbFile <<atoms[i].resId<<"     ";
    outpdbFile.width(8);
    outpdbFile.precision(4);
    outpdbFile <<std::fixed<<atoms[i].r.x<<"  ";
    outpdbFile.width(8);
    outpdbFile.precision(4);
    outpdbFile <<std::right<<std::fixed<<atoms[i].r.y<<"  ";
    outpdbFile.width(8);
    outpdbFile.precision(4);
    outpdbFile <<std::right<<std::fixed<<atoms[i].r.z<<"  ";	   
    outpdbFile.precision(2);
    outpdbFile <<std::fixed<<atoms[i].occupancy<<"   ";
    outpdbFile <<std::fixed<<atoms[i].tempFactor<<"      ";
    outpdbFile.width(5);		   
    outpdbFile <<std::left<<atoms[i].segName<<"   ";
    outpdbFile.width(1);
    outpdbFile<<std::right<<atoms[i].name[0]<<"\n";
  }
  
  outpdbFile.close();
}
/**===============================================================*/

/**===============================================================
 * \c writePhi(int dihed):  write the dihedral angle
 * dihed: is the id of the dihedral angle (the id start at 1, 2, ..., n)
 * 
 * The output format for Binary is (phi_t1 phi_t2 phi_t3.......)
 * where subscript ti represent the time at which the angle is printed and there is
 * no space between successive print, so we need to read it in binary using size of 
 * double.
 * 
 * For non-binary, all phis are printed in a line separated bu space, hence it can be 
 * read easily in application such as matlab.
 * 
 *
 */
void mdpack::md::writePhi(int dihed){
  if(isWritePhiOn){
    double var = dihedrals[dihed -1].phi;
    if(isBinaryOn)
      phiFile.write((char *) &var, sizeof var);
    else
      phiFile << var <<" ";
  }
}
/**=================================================================*/

/**===============================================================
 * \c writePhi():  write set of dihedral angles, ginve by array of ids
 * diheds: is integer array ids of the dihedral angles (the id start at 1, 2, ..., n)
 * 
 * The output format for Binary is on one line. For example for two dihedrals we will 
 * have (phi1_t1 phi2_t1 phi1_t2 phi2_t2.......)
 * where subscript ti represent the time at which the angle is printed and there is
 * no space between successive print, so we need to read it in binary using size of 
 * double.
 * 
 * For non-binary, all phis are printed in a line separated bu space, and each line
 * is a different print time, for example for two dihedrals we will have:
 * phi1  phi2
 * phi1  phi2
 * ..........
 * .........
 * phi1  phi2
 * 
 * Hence it can be read easily in application such as matlab.
 * 
 *
 */
void mdpack::md::writePhi(const int * diheds, const int& num){
  if(isWritePhiOn){
    double var;
    if(isBinaryOn){
      for (int i=0; i<num; i++){
	var =dihedrals[diheds[i] -1].phi;
	phiFile.write((char *) &var, sizeof var);
      }
    }
    else{
      for(int i=0; i<num; i++){
	var =dihedrals[diheds[i] -1].phi;
	phiFile << var <<" ";
      }
      phiFile <<"\n";
    }
  }
}
/**=================================================================*/

// computing kinetic energy
void mdpack::md::computeKinetic(void){
  Ekinetic=0.0;
  for(int i=0; i<nAtoms; i++){
    Ekinetic += 0.5*atoms[i].m*(atoms[i].v.x*atoms[i].v.x +
                      atoms[i].v.y*atoms[i].v.y +
                      atoms[i].v.z*atoms[i].v.z);
  }
}

// computing kinetic energy and instantaneous temperature
void mdpack::md::computeTemperature(void){
  tempInstant =0.0;
  computeKinetic();
  tempInstant = 2.0*INV_BOLTZMANN*(1.0/nDegsFreedom)*Ekinetic;
}
  
  
  
  
/**=================  initVelocities =========================
 * Initializing velocities to the target temperature
 */      
void mdpack::md::initVelocities(void){
  mdpack::dreal mx=0.0;
  mdpack::dreal my=0.0;
  mdpack::dreal mz=0.0;
 
  // loop 1 and loop 2 are setting tptal linear momentum to be zeros
  for(int i=0; i<nAtoms; i++){
    atoms[i].v.x= W();
    atoms[i].v.y= W();
    atoms[i].v.z= W();
    mx += atoms[i].m*atoms[i].v.x;
    my += atoms[i].m*atoms[i].v.y;
    mz += atoms[i].m*atoms[i].v.z; 
  }
  mdpack::dreal inv_n= 1.0/nAtoms;
  mx =inv_n*mx; my =inv_n*my; mz= inv_n*mz;
  for(int i=0; i<nAtoms; i++) {
    atoms[i].v.x -= atoms[i].im*mx;
    atoms[i].v.y -= atoms[i].im*my;
    atoms[i].v.z -= atoms[i].im*mz;
  }

  // setting up initial temperature to target temperature
  computeTemperature(); // setup instantaneous temperature (tempInstant)
  double fct = sqrt(targetTemp/tempInstant);
  for(int i=0; i<nAtoms; i++) {
    atoms[i].v.x *= fct;
    atoms[i].v.y *= fct;
    atoms[i].v.z *= fct;
  }
  // setup initial values of kinetic energy and temperature
  computeTemperature();
  tempAvg = tempInstant;
}
/**==============================================================*/
  

/**=================   init =========================
 * Initializing velocities and position from pdb files
 */ 
void mdpack::md::init(const char * pos, const char * vel){
  std::ifstream posFile;
  posFile.open(pos, std::ios::in);
  if(!posFile)
    terminate(" cant open restart file for position ");
  std::ifstream velFile; 
  velFile.open(vel, std::ios::in); 
  if(!velFile)
    terminate(" cant open restart file for velocities ");
  
  mdpack::dreal x, y, z;
  char part1 [50];
  char part2 [50];
  int nParts, id;
  std::string line;
  int checkAtoms =0;
  // reading in positions
  while(!posFile.eof()){
    getline(posFile,line);
    if(getWord(line,part1, 1) && strcmp(part1,"ATOM")==0){
      nParts = countWords(line);
      if(nParts == 13){
	getWord(line,part2, 2);
	id= atoi(part2);
	getWord(line,part2, 7);
	x= atof(part2);
	getWord(line,part2, 8);
	y= atof(part2);
	getWord(line,part2, 9);
	z= atof(part2);
      }
      else if (nParts == 12){
	getWord(line,part2, 2);
	id= atoi(part2);
	getWord(line,part2, 6);
	x= atof(part2);
	getWord(line,part2, 7);
	y= atof(part2);
	getWord(line,part2, 8);
	z= atof(part2);
      }
      else
	terminate(" cannot read restart positions pdb file ");
      
      atoms[id-1].r.x=x; 
      atoms[id-1].r.y=y;
      atoms[id-1].r.z=z;
      checkAtoms +=1;
    }
  }
  
  if(checkAtoms != nAtoms)
    terminate("  cannot find coordinates for atoms from restart positions pdb file ");
  posFile.close();
  
  checkAtoms =0;
  // reading in velocities
  while(!velFile.eof()){
    getline(velFile,line);
    if(getWord(line,part1, 1) && strcmp(part1,"ATOM")==0){
      nParts = countWords(line);
      if(nParts == 13){
	getWord(line,part2, 2);
	id= atoi(part2);
	getWord(line,part2, 7);
	x= atof(part2);
	getWord(line,part2, 8);
	y= atof(part2);
	getWord(line,part2, 9);
	z= atof(part2);
      }
      else if (nParts == 12){
	getWord(line,part2, 2);
	id= atoi(part2);
	getWord(line,part2, 6);
	x= atof(part2);
	getWord(line,part2, 7);
	y= atof(part2);
	getWord(line,part2, 8);
	z= atof(part2);
      }
      else
	terminate(" cannot read restart velocities pdb file ");
      
      atoms[id-1].v.x=x; 
      atoms[id-1].v.y=y;
      atoms[id-1].v.z=z;
      checkAtoms +=1;
    }
  }
  
  if(checkAtoms != nAtoms)
    terminate("  cannot find coordinates for atoms from restart velocities pdb file ");
  velFile.close();
  
  computeTemperature();
  tempAvg = tempInstant;
  
}
/**=========================================================*/
  
  
  
  
  
  
  
  
  
  
  
  
    
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
    