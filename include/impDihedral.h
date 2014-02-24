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

#ifndef IMPDIHEDRAL_H
#define IMPDIHEDRAL_H

#include <cstring>
#include <cstdio>

#include "definitions.h"
//#include "operations.h"
namespace mdpack {
  
  /**
   * class impDihedralParams: provide parameters and information needed to 
   * compute improper dihedral potential V(improper) and forces:
   * 
   * V(improper) = Kpsi(psi - psi0)**2
   * Kpsi: kcal/mole/rad**2
   * psi0: degrees
   * 
   */
  
  class impDihedral {
    public:
      //index for first atom: atom1 correspond to atom id in the atom array
      int atom1;  
      //index for second atom: atom1 correspond to atom id in the atom array
      int atom2;  
      //index for third atom: atom1 correspond to atom id in the atom array
      int atom3;  
      //index for fourth atom: atom1 correspond to atom id in the atom array
      int atom4;
      //------------------------------------------
      /* type1, type2, type3 and type4 define three consecutive bonds: 
       * type1--type2--type3--type4
       */
      mdpack::atomType type1;
      mdpack::atomType type2;
      mdpack::atomType type3;
      mdpack::atomType type4;
      //-------------------------------------------
      
      mdpack::dreal phi;
      
      //-----------------------------------------
      /* parameters for computing angle potential:
       * V(improper) = K(psi - delta)^2
       */
      mdpack::dreal k;
      mdpack::dreal delta;
      //-----------------------------------------
      
      // true if parameters are set 
      bool isSet;  
      
      // constructors
      impDihedral();
      impDihedral(int a1, int a2,int a3, int a4, const char * str1, const char * str2,
	    const char * str3, const char * str4, const mdpack::dreal& k_phi, 
	    const mdpack::dreal& d);
      
      //copy constructor
      impDihedral(const mdpack::impDihedral& d);
      
      // assignment operator
      impDihedral& operator=(const mdpack::impDihedral& d);
      
      void set(int a1, int a2, int a3, int a4, const char * str1, 
	       const char * str2, const char * str3, const char * str4);
      
      void setPrm(const mdpack::dreal& k_phi,const mdpack::dreal& d);
      
      bool compare(const char * str1, const char * str2, const char * str3,
		   const char * str4);
      
      bool operator==(const impDihedral& d);
      bool operator!=(const impDihedral& d);
      
      void print(int id){
	printf("IMPDIHEDRAL  %-7d  atom1 = %-7d  atom2 = %-7d  atom3 = %-7d  atom4 = %-7d  ",
	       id,atom1, atom2,atom3, atom4);
	printf("%-4s  %-4s  %-4s  %-4s  k = %-6.3f  delta = %-6.3f\n",
	       type1, type2, type3, type4, k, delta );
      }
      
  };
    
  
}; //end of namespace mdpack



#endif  //end of IMPDIHEDRAL_H