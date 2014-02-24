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

#ifndef BOND_H
#define BOND_H

#include <cstring>
#include <iostream>
#include <cstdio>


#include "definitions.h"
//#include "operations.h"
namespace mdpack {
  /**
  * class bond: provide all information and functions that is needed to
  * compute bonded potential V(bond) and bonded forces:
  * 
  * V(bond) = K(b - r0)^2
  * K: kcal/mole/A^2
  * 
  * \param type1 it corresponds to a bond defined by type1--type2
  * \param type2
  * type1 and type2 are used to define the bond type1--type2
  * the parameters for the bond type1--type2 are obtained from the 
  * user defined input parameter file (this version only accept charmm parameters)
  * 
  */
  
  class bond{
    public:
      //index for first atom: atom1 correspond to atom id in the atom array
      int atom1;  
      //index for second atom: atom1 correspond to atom id in the atom array
      int atom2;
      // atom types for bond type1--type2
      mdpack::atomType type1;
      mdpack::atomType type2;
      // bond parameters for the bond (type1--type2)
      mdpack::dreal k;
      mdpack::dreal r0;
      // boolian to check if parameters k, and r0 are set 
      bool isSet;
      
      bond();
      bond(int a1, int a2,const char * str1, const char * str2, mdpack::dreal& kapa,
                  mdpack::dreal& b0);
      bond(const bond& b);
      bond& operator=(const bond& b);
      void set(int a1, int a2, const char * str1, const char * str2 );
      void setPrm(const mdpack::dreal& kapa, const mdpack::dreal& b0){k=kapa; r0=b0, isSet= true;}
      
      bool compare(const char * str1, const char * str2);
      bool operator==(const bond& b);
      bool operator!=(const bond& b);  
      void print(int id){ 
	printf("BOND  %-7d  atom1 = %-7d  atom2 = %-7d    %-4s  %-4s  k = %-6.3f  r0 = %-5.4f\n",
	       id,atom1, atom2,type1, type2, k, r0);
      }
    
  };
  
}; //end of namespace mdpack
    



#endif  //end of BOND_H






