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

#ifndef LJTABLE_H
#define LJTABLE_H

#include <cstring>
#include <cstdio>

#include "definitions.h"

namespace mdpack {
  
  /**
   * class ljTable provides parameters for Lennard-Jones potential
   * between pairs of atom i j
   * V(|rij|) = A/|rij|^12 - B/|rij|^6
   * rij = rj-ri  distance from i to j
   */
  class ljTable {
    public:
      /**
       * \c inc
       * 1 if the interaction is included, sero otherwise
       * it value depend on the exclude option
      */
      int inc;
      // atom id for atom i
      int id1;
      // atom id for atom j
      int id2;
      mdpack::dreal A;
      mdpack::dreal B;
      ljTable(){ inc=0; A=0.0; B=0.0; id1=0; id2=0;}
  
      // copy constructor
      ljTable(const mdpack::ljTable& LJ);
      ljTable& operator=(const ljTable& LJ);
      
      
      
      // print information into std output
      void print();
    
  };
  
};

#endif // LJTABLE_H