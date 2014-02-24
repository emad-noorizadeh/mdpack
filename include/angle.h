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

#ifndef ANGLE_H
#define ANGLE_H


#include <cstring>
#include <cstdio>

#include "definitions.h"
//#include "operations.h"
namespace mdpack {
  
  /**
   * class angle: provide parameters and information needed to 
   * compute angle potential V(angle) and forces:
   * V(angle) = Ktheta(Theta - Theta0)^2
   * V(Urey-Bradley) = Kub(S - S0)^2
   * Ktheta: kcal/mole/rad^2
   * Theta0: degrees
   * Kub: kcal/mole/A^2 (Urey-Bradley)
   * S0: A
   */
  class angle {
    public:
      //index for first atom: atom1 correspond to atom id in the atom array
      int atom1;  
      //index for second atom: atom1 correspond to atom id in the atom array
      int atom2;  
      //index for third atom: atom1 correspond to atom id in the atom array
      int atom3;  
      //------------------------------------------
      /* type1, type2 and type3 define two consecutive bonds: 
       * type1--type2--type3
       */
      mdpack::atomType type1;
      mdpack::atomType type2;
      mdpack::atomType type3;
      //-------------------------------------------
      
      //-----------------------------------------
      /* parameters for computing angle potential:
       * V(angle) = Ktheta(Theta - Theta0)^2
       */
      mdpack::dreal ktheta;
      mdpack::dreal theta0;
      mdpack::dreal kub;
      mdpack::dreal rub;
      //-----------------------------------------
      
      // true if parameters are set 
      bool isSet;  
      
      // constructors
      angle();
      angle(int a1, int a2,int a3, const char * str1, const char * str2,
	    const char * str3, const mdpack::dreal& k_theta, const mdpack::dreal& theta_0,
	    const mdpack::dreal& k_ub, const mdpack::dreal& r_ub);
      // copy constructor
      angle(const angle& a);
      
      // assignment operator
      angle& operator=(const angle& a);
      
      void set(int a1, int a2, int a3,const char * str1, const char * str2, const char * str3);
      void setPrm(const mdpack::dreal& k_theta, const mdpack::dreal& theta_0,
			   const mdpack::dreal& k_ub, const mdpack::dreal& r_ub);
      
      bool compare(const char * str1, const char * str2, const char * str3);
      
      bool operator==(const angle& a);
      bool operator!=(const angle& a); 
      
      void print(int id){
	printf("ANGLE  %-7d  atom1 = %-7d  atom2 = %-7d  atom3 = %-7d  ",id,atom1, atom2,atom3);
	printf("%-4s  %-4s  %-4s  ktheta = %-6.3f  theta0 = %-7.4f  kub = %-6.3f rub = %-7.4f\n",
	       type1, type2, type3, ktheta, theta0, kub, rub );
      }
      
      
  };
    
  
}; //end of namespace mdpack



#endif  //end of ANGLE_H
