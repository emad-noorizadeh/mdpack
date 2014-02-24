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

#ifndef ATOM_H
#define ATOM_H

#include <cstring>
#include <list>
//#include <iostream>

#include "definitions.h"

#include "vector3.h"
namespace mdpack {
  
  /**
   * class atom: provide parameters and information for each individual atom 
   * 
   */
  
  class atom {
    public:
      // type of the atom
      mdpack::atomType type;  
      // name of the atom
      mdpack::atomType name;     
      // name of the residue that the atom belong to
      mdpack::atomType resName;  
      // name of the segment that the atom belong to
      mdpack::segmentName segName; 
      // residue id in PDB file
      int resId;      
      // atoms id in PDB or PSF file
      int id;    
      // vector of positions and velocities
      mdpack::vector3 r; 
      mdpack::vector3 v;
      
      // occupancy ( from PDB file)
      mdpack::dreal occupancy;
      // temperature factor ( from PDB file)
      mdpack::dreal tempFactor;
      
      // parameters for Lennard Jones interaction
      mdpack::dreal halfSigma;      
      mdpack::dreal eps;
      
      // charg of the atom from PDB or PSF file
      mdpack::dreal q;      
      // mass of the atom
      mdpack::dreal m; 
      // inverse mass of the atom 1/m
      mdpack::dreal im;
      // true if LJ parameters has been set
      bool isSet;
      
      // list of atoms that are connected to this atom by a linear bond
      std::list<int> list1_2;
      // list of atoms that are connected to this atoms by two linear bonds
      std::list<int> list1_3;
      // list of atoms that are connected to this atoms by three linear bonds
      std::list<int> list1_4;
      
      
      
      atom();
      
      atom(int aId, const char * sname, int rId, const char * rname, const char * n,
	   const char * t, const mdpack::dreal& charge, const mdpack::dreal& mass,
	   mdpack::dreal& e,mdpack::dreal& sig, mdpack::vector3& coord ,
	   mdpack::dreal& oc, mdpack::dreal& temp);
      
      atom(const atom& a);
      
      atom& operator=(const atom& a);
      
      void set(int aId, const char * sname, int rId, const char * rname, const char * n,
	   const char * t, const mdpack::dreal& charge, const mdpack::dreal& mass);
      
      void setPrm(const mdpack::dreal& e,const mdpack::dreal& sig ){eps=e; halfSigma=sig;isSet = true;}
      
      void setCoords(const mdpack::dreal& x, const mdpack::dreal& y,
		  const mdpack::dreal& z ){ r.x=x; r.y=y; r.z=z;}
      
      bool operator==(const atom& a);
      
      bool operator!=(const atom& a);
      
      bool isIn1_2(const int& k);
      bool isIn1_3(const int& k);
      bool isIn1_4(const int& k);
      
      void setTempFactor(mdpack::dreal& d){ tempFactor = d;}
      void setOccupancy(mdpack::dreal& d){ occupancy = d;}
      
      void printToPSF() { 
	printf("%8d  %5s  %4d  %5s  %4s  %4s  %12.6f  %12.4f    %1d\n", 
				 id, segName, resId, resName, name, type, q, m,0);}                  
                         
      void printToPDB(){
	printf("ATOM  %8d  %-4s  %-5s  %4d      %+4.4f  %+4.4f  %+4.4f  %3.1f  %3.1f       %-5s  %-1c \n",
	id, name,resName,resId, r.x,r.y, r.z,occupancy,tempFactor, segName,name[0]);} 
			   
			   
			   
  }; // end of atom
    
  
}; //end of namespace mdpack



#endif  //end of ATOM_H
