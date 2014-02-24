/*
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

#ifndef MD_H
#define MD_H

#include <math.h>
#include <list>
//#include <algorithm>

#include "definitions.h"
#include "utilities.h"
#include "constants.h"
#include "atom.h"
#include "bond.h"
#include "angle.h"
#include "dihedral.h"
#include "impDihedral.h"
#include "ljTable.h"
#include "cell.h"
#include "random.h"



namespace mdpack {
  
  
  
  /**
   * class md: this is a core class for setting up and running 
   * molecular dynamics (MD) simulation
   * 
   */
  
  class md {
    
    /**==============================================================
    *                  Constructors and Destructor
    * ===============================================================*/
  public:
    // empty constructor
    md(){ initVars();}
	
    md(const char * pdb, const char * psf, const char * pars){
      initVars();
      readinpsf(psf);
      readinprm(pars);
      readinpdb(pdb);
      setup();
      initVelocities();
    }
    
    void initVars(void);
      
    ~md(){ 
      if (atoms != NULL) delete [] atoms; 
      if (bonds != NULL) delete [] bonds; 
      if (angles != NULL) delete [] angles;
      if(dihedrals != NULL) delete [] dihedrals; 
      if(impDihedrals !=NULL) delete [] impDihedrals;
      if(f != NULL) delete [] f;

      if(isNonBondedPairs){
	for (int i = 0; i < nAtoms; ++i)
	  if(nonBondedPairs[i] != NULL) 
	    delete [] nonBondedPairs[i];
	if(nonBondedPairs !=NULL)
	  delete [] nonBondedPairs;
      }
	     
     if(isCells)  delete [] cells;
     if(langFactor !=NULL) delete [] langFactor;
     if(langFactorh !=NULL) delete [] langFactorh;
   }
    
      
    /**==============================================================
    *                       System variables
    * ===============================================================*/
  private:
    
    // variables for electorstatic forces
    mdpack::dreal e0;
    mdpack::dreal inv_e0;
    mdpack::dreal e14;
      
    // random generator class
    randn W;
    // seed for initializing the random generator
    int seed;
  public:
    /** array of atoms, e.g., atoms[i] is the atom with id i+1,
      *  (see class atoms for more details)
      */
    mdpack::atom * atoms;
    /** array of bonds, e.g., bonds[i] contains bond parameters and information
      *  for bond i+1 (see class bond for more details)
      */
    mdpack::bond * bonds;
    /** array of angles, e.g., angles[i] contains angle parameters and information
      *  for angle i+1  (see class angle) 
      */
    mdpack::angle * angles;
    // array of dihedrals angle (see class dihedral for more details)
    mdpack::dihedral * dihedrals;
    // array of improper dihedral angles
    mdpack::impDihedral * impDihedrals;
    // array of ljTable for all lennard Jones interation (see class ljTable)
    mdpack::ljTable ** nonBondedPairs;
    // a list of ljTable for all lennard Jones interation
    // std::list<mdpack::ljTable> nonBondedPairsList;
    mdpack::cell * cells;
      
     
//     // array of position vector3s of all atoms
//     mdpack::vector * r;
//     // array of velocitie vectors of all atoms
//     mdpack::vector * v;
//     // array of momentum vectors of all atoms
//     mdpack::vector * p;
    // array of force vectors of all atoms
    mdpack::vector3 * f;
      
    int nAtoms;                   // number of atoms
    int nBonds;                   // number of bonds
    int nAngles;                  // number of angles
    int nDihedrals;               // number od dihedrals
    int nImpDihedrals;            // number of improper diherals
    int nDegsFreedom;            // number of degrees of freedom
      
      
    int nCells;  // number of cells used in cell lists for nonBonded computation
    mdpack::dreal boxSide;   // the length of simulation box, the box is divided into cells
      
    // the initial max values of coordinates
    mdpack::dreal max_x, max_y, max_z, min_x, min_y, min_z;
      
      
    // mass of centre of mass ( sum of all masses)
    mdpack::dreal centerOfMass_m; 
    // coordinates of center of mass
    mdpack::vector3 centerOfMassCoords;
      
    // energies
    mdpack::dreal EnonBondeds;
    mdpack::dreal Ebonds;
    mdpack::dreal Eangles;
    mdpack::dreal Edihedrals ;
    mdpack::dreal EimpDihedrals;
    mdpack::dreal EnonBonded ;
    mdpack::dreal Ekinetic;
    mdpack::dreal Epotential;
    mdpack::dreal Etotal;
    mdpack::dreal EsphBC;
    mdpack::dreal Evdw;
    mdpack::dreal Eelec;
    mdpack::dreal Emisc;
     
      
    /**==============================================================
      *             System and simulation conditional variables
      * ===============================================================*/
  public:  
    /** indicate which atoms are excluded from non-bonded interactions. Values are:
      * none: means all atoms are included in nonBonded interaction
      * 1-2: atoms that are connected by a linear bond will be excluded
      * 1-3: atoms that are connected by two linear bonds will be excluded
      * 1-4: atoms that are connected by three linear bonds will be excluded.
      * This parameter can be set in config file as:
      * exclude   on or off
      */
    char exclude [5];
      
    /** if true we  will compute electorstatic and lennard-JOnes forces using 
      * cutoff and smoothing function. This parameter can be set in config file as:
      * switching   on or off
      */
    bool switching; 
    /**
      * If true, then we use cells list to compute nonBonded forces.
      * if system size (i.e., number of atoms) is less than MINATOMS
      * then this will be ignored and cellslist would not be used. 
      * This parameter can be set in config file as:
      * cellsList   on or off
      */
    bool useCellLists;
      
    /** this is used to delete alocated memory of nonBondedPairs, true if nonBondedPairs 
      *is alocated, false otherwise.
      */
    bool isNonBondedPairs;
    
    /** this is used to delete alocated memory of cells (for cell lists), 
     * true if cells is alocated, false otherwise.
     */
    bool isCells;
    
      
    /**==============================================================
      *                      Simulation variables
      * ===============================================================*/
  public:
    // cutoff radious  for nonbonded interaction, default value 14.0
    mdpack::dreal cutoff;  
    // distance at which to activate switching function for nonbonded interaction
    //default value 13.0
    mdpack::dreal switchdist;
    // number of steps for the simulation run, default value 0
    int numSteps; 
    // the nth step
    int stepn;
    // time step-size for integration, default value TIMEFACTOR*0.001
    mdpack::dreal h;
    // time step
    mdpack::dreal timeStep;
    // initial time of simulation, default value 0.0
    mdpack::dreal initialTime;
    // simulation time
    mdpack::dreal simTime;
      
    // temperature at each step
    mdpack::dreal tempInstant;
    // average temperature
    mdpack::dreal tempAvg;
      
    // perscribed temperature of the system given by user
    mdpack::dreal targetTemp;
    
    // simulation name given by user
    char simName[100];
      
    /**==============================================================
      *             System and simulation setup functions
      * ===============================================================*/
  public:
    void initVelocities(void);
    void init(void){ initVelocities();}
    void init(const char * pos, const char * vel);
    
    /**
     * Functions isIn1_2, isIn1_3, isIn1_4, are used to exclude pairs that
     * are connected by bonds from non bonded interaction.
     * \c isIn1_2: Returns true if atoms with id1 and id2 are connected by a linear bonds.
     * \c isIn1_3: Returns true if atoms with id1 and id2 are connected by two linear bonds.
     * \c isIn1_4: Returns true if atoms with id1 and id2 are connected by three linear bonds.
     */
    bool isIn1_2(const int& id1, const int& id2);
    bool isIn1_3(const int& id1, const int& id2);
    bool isIn1_4(const int& id1, const int& id2);
    
    void setup();
    void setLJtable(void);
    void setupCellLists(void);
      
      
    /**==============================================================
     *           Temporary variables for force computation
    * ===============================================================*/
  private:
    mdpack::vector3 rij, rjk, rkj, rkl, rik, mm, nn, R, S, rr, fi, fl, fk, rji;
    mdpack::dreal dist_ij, dist_jk, dist_kj, dist_kl, dist_ik, KK;
    mdpack::dreal inv_dist_ij, inv_dist_jk, inv_dist_kj, inv_dist_kl, inv_dist_ik;
    mdpack::dreal inv_dist_R, inv_dist_S;
    mdpack::dreal dist, inv_dist ,inv_dist6, inv_dist12, dV, qi, qj ;
    mdpack::dreal factor1, factor2, factor3, factor4, factor5;
    mdpack::dreal factor6, factor7, factor8, factor9, factor10;
    mdpack::dreal temp1, temp2;
    mdpack::dreal inv_Rc2;
    mdpack::dreal fct;
    mdpack::dreal fix, fiy, fiz;
    mdpack::dreal fkx, fky, fkz;
    mdpack::dreal flx, fly, flz;
    int idx1, idx2, idx3, idx4;
    mdpack::dreal phi;
    mdpack::dreal cos_theta, sin_theta, atheta, ;
    mdpack::dreal inv_n, inv_m, norm_rkj,normR, normS,dot_R_S,dot_kj_ij,dot_kj, dot_ij_kj;
    mdpack::dreal inv_norm_rkj_ij,inv_norm_rkj,dot_kl_kj,inv_normR, inv_normS, drij,drkj,drik;
    mdpack::dreal dV_dtheta,inv_dist_rjk, inv_dist_rji;
      
      
    /**==============================================================
      *                Functions for force computation
      * ===============================================================*/
  public:
    void computeNonBondeds(void);
    void computeBonds(void);
    void computeAngles2(void);
    void computeAngles(void);
    void computeDihedrals2(void);
    void computeDihedrals(void); 
    void computeImpDihedrals(void);
    void computeImpDihedrals2(void);
    void computeSphBC(void);
    void computeForces(void);
    /**==============================================================
      *                   Computational functions
      * ===============================================================*/
  public:
    // computes kinetic energy
    void computeKinetic(void);
    //computes kinetic and instantaneous temperature
    void computeTemperature(void);
    // moves atoms to the center of mass
    void moveToCM(void);
      
      
    /**==============================================================
      *                      Integration methods
      * ===============================================================*/
  private:
    // constant factors for Langevin methods
    mdpack::dreal expGamma_h ;
    mdpack::dreal * langFactor;
    mdpack::dreal expGamma_hh ;
    mdpack::dreal * langFactorh; 
    
    // constant factors for NHL method
    mdpack::dreal nhl_imu; // inverse of mu
    mdpack::dreal expNhl;  // constat factor used in NHL method
    mdpack::dreal nhlFactor; // constant factor used in NHL
    
    //constat factor for over-damped dynamics
    mdpack::dreal overDampedFactor;
    
  public:
    void verlet(void);
    void langevin(void);
    void langevin1(void);
    // gamma parameter for Langevin
    mdpack::dreal gamma;
    
    // Nose-Hoover-Langevin (NHL) thermostat
    void noseHooverLangevin(void);
    // gamma parameter for (NHL)
    mdpack::dreal nhl_gamma;
    // mu parameter for NHL
    mdpack::dreal nhl_mu;
    // NHL xi variable
    mdpack::dreal nhl_xi;
    
    //over damped dynamics
    void overDamped(void);
    // Gradient descent minimization
    void minimize(void);
      
      
    /**==============================================================
      *         Parameters for shperical boundary condition
      * ===============================================================*/
  public:
      
    /**
      * \c sphBC: if ture the system uses shperical boundary codition that is
      * a potential of 
      * Vsph (ri) = sphBCk1(|ri- sphBCCenter|- sphBCr1)^sphBCexp1
      * acts on all atoms
      * ri: is the position vector of atom i
      * Default value: false, it can be changed using config file as
      * sphericalBC   on or off
      */
    bool sphBC;
      
    /**
      * \c sphBCCenter: Center of the shperical boundary condition, it is 
      * set by program to the center of mass of the system
      */
    mdpack::vector3 sphBCCenter;
      
    /**
      * \c sphBCr1: Radius of the shperical boundary condition, must be specified
      * by user in config file if spherical boundary condition is on. This can be 
      * in config file as
      * sphericalBCr1  value  ( in A )
      */
    mdpack::dreal sphBCr1;
      
    /**
      * \c sphBCk1: Force constant, must be specified by user in config file 
      * if spherical boundary condition is on. This can be 
      * in config file as
      * sphericalBCk1  value  (in A )
      */
    mdpack::dreal sphBCk1;
      
    /**
      * \c sphBCexp1: Exponent of shperical boundary condition potential. This is
      * optional, aceptable values are 2 or 4, can be set in config file as
      * sphericalBCexp1  2 or 4
      * Default is 2
      */
    mdpack::dreal sphBCexp1;
      
    /* these parameters are for spherical boundary condition with
      * tension effect and are not implemented in this version
      * Vsph (ri) = sphBCk1(|ri- sphBCCenter|- sphBCr1)^sphBCexp1 -
      *             sphBCk2(|ri- sphBCCenter|- sphBCr2)^sphBCexp2
      * in this case usually  sphBCexp1=4,  sphBCexp2 =2
      * and sphBCr1 > sphBCr2
      */
    /*
    mdpack::dreal sphBCr2;
    mdpack::dreal sphBCk2;
    mdpack::dreal sphBCexp2;
    */
      
      
    /**==============================================================
      *                   Parameters for Input Files 
      * ===============================================================*/
  public:
      
    /**
      * \c pdbFile: this is a PDB file containing initial position of all
      *atoms and must be specified by the user in the config file. It is specified by user in
      * the config file as:
      * coordinates   filename 
      * (note that the PDB file must be in NAMD PDB format, if you have a PDB in protein data bank 
      * format you must first use psfgen VMD to convert it to NAMD PDB format:
      *  http://www.ks.uiuc.edu/Training/Tutorials/namd/namd-tutorial-unix-html/node20.html
      * )
      * Acceptable Values: UNIX filename
      */
    char pdbFile [100];
      
    /**
      * \c structure: this is a PSF file that should be specified by user in the 
      * config file as:
      * structure   nameOfYourPSFfile
      * Acceptable Values: UNIX filename
      */
    char psfFile [100];
      
    /**
      * \c parameters: this is a parameter file that must be specified by user
      * in the config file as
      * parameters   nameOfParametersFile
      * (note that this version only accept the current Chamrm parameters format  CHARMM22 or
      * CHARMM36, the cross terms CAMP is not implemented in this version v 0.01)
      * Acceptable Values: UNIX filename
      */
    char prmFile [100];
      
    /**
      * \c velocities: this is a  PDB file containing initial velocities of all atoms. This is 
      * optional and can be specified by user to restart a simulation at given positions and velocities.
      * This can be specified in the config file as
      * velocities   filename
      * (note that as the for coordinates the PDB file is in NAMD PDB format
      * http://www.ks.uiuc.edu/Training/Tutorials/namd/namd-tutorial-unix-html/node20.html
      * ) 
      * Acceptable Values: UNIX filename
      */
    std::string velocities;
         
      
    /**==============================================================
      *                  Functions for reading input files 
      * ===============================================================*/
  public:
    // reading in psf file and seting up atoms, bonds, angles, dihedrals and impDihedrals
    void readinpsf(const char * psf);
    // reading in prm file and setting up parameters
    void readinprm(const char * prm);
    // reading in pdb file and setting up coordinates for atoms
    void readinpdb(const char * pdb);
    // reading config file and setting up variables and parameters from config file
    void readinConfig(const char * conf);
    
      
    /**==============================================================
      *                   Parameters for Output Files 
      * ===============================================================*/
  public:
    /**
      * \c outputName: This is the name of a PDB file to write the final coordinates,
      * velocities and momentum of all atoms.
      * This file then can be used as an input coordinates file to start a simulation. This can
      * be specified in config file as
      * outputName  filename
      * (note that the filename will be appended by .coor)
      */
    char outputName [50];
      
    std::ofstream energyLogFile;
    std::ofstream posFile;
    std::ofstream velFile;
    std::ofstream phiFile;
    std::ofstream outpdbFile;
    std::ofstream restartPosFile;
    std::ofstream restartVelFile;
    // if true we will write positions into a file
    bool isWritePosOn;  
    // if true we will write velocities into a file
    bool isWriteVelOn;
     // if true we will write dihedral angle or angles into a file
    bool isWritePhiOn;
    // if true we will write output in binary format
    bool isBinaryOn;
      
    /**==============================================================
      *               Functions for writing output files and std output 
      * ===============================================================*/
  public:
    void openOutputFiles(void);
	
    void closeOutputFiles(void){
      if(energyLogFile.is_open())
	energyLogFile.close();
      if(posFile.is_open())
	posFile.close();
      if(velFile.is_open())
	  velFile.close();
      if(phiFile.is_open())
	  phiFile.close();
    }
    void writeEnergies(void);
    void writePositions(void);
    void writeVelocities(void);
    void writePhi(int dihed);
    void writePhi(const int * diheds, const int& num);
    void writeRestart(const int& step_i,const mdpack::dreal& t);
    void writepdb(const int& step_n, const mdpack::dreal& t );
    void printEnergies(void);
    void printEtitles(void);
    void printInfo(void);
      
     
     
  }; // end of md
    
  
}; //end of namespace mdpack



#endif  //end of MD_H