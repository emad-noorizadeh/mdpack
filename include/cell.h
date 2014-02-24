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

#ifndef CELL_H
#define CELL_H

#include <cstring>
#include <iostream>
#include <cstdio>
#include <list>

#include "atom.h"
#include "vector3.h"
#include "definitions.h"
//#include "operations.h"
namespace mdpack {
  /**
  * class cell: the simulation domain is divided into computation cells
  * each cell has list of neighbour that are within cutoff distance from its centre
  * the computation of non-bonded forces is done over the cells instead of atoms,
  * hence reducing the computational cost.
  
  */
  
  class cell{
    public:
//       //list of atoms in the cell
//       mdpack::atom * atoms;
//       // list of atoms that left the cell
//       mdpack::atom * exitAtoms;
//       // list of atoms that entered the cell
//       mdpack::atom * enterAtoms;
//       
//       // list of neighbouring cells
//       mdpack::cell * neighbours;

      //list of atoms in the cell
      std::list<int>  atomIds;
     // std::list<int> listOfAtomIds;
      // list of atoms that left the cell
//       std::list<int> * exitAtomIds;
      std::list<int> exitAtomIds;
      // list of atoms that entered the cell
//       int * enterAtomIds;
//       std::list<int> listOfEnterAtomIds;
      // list of neighbouring cells
      std::list<int> neighbours;
      //std::list<int> listOfNeighbour;
      int nNeighbours;
      
      // number of atoms in the cell
      int size;
      // number of atoms that left the cell
      //int exitSize;
      // number of atoms that entered the cell 
     // int enterSize;
      
      // center of the cell
      mdpack::dreal x, y, z;
      // dimension of the cell
      mdpack::dreal lx, ly, lz;
      // corners 
      mdpack::dreal leftx, rightx, lefty, righty, leftz, rightz;
      
      // array of random variables for future MPI implementation
      // mdpack:dreal * w;
      
      cell(){ size= 0;  x=0.0; y=0.0; z=0.0; lx=0.0; ly =0.0; lz=0.0;nNeighbours =0;
	     leftx=x-0.5*lx; rightx=x+0.5*lx; lefty=y-0.5*ly; righty=y+0.5*ly; 
	     leftz=z -0.5*lz; rightz=z+ 0.5*lz;
      }
      cell(const mdpack::dreal& cx,const mdpack::dreal& cy, const mdpack::dreal& cz, const mdpack::dreal& L){
	     size =0; x=cx; y=cy; z=cz; lx=L; ly =L; lz=L;
	     leftx=x-0.5*lx; rightx=x+0.5*lx; lefty=y-0.5*ly; righty=y+0.5*ly; 
	     leftz=z -0.5*lz; rightz=z+ 0.5*lz;
	     
	     nNeighbours =0;
      }
	   
      cell(const mdpack::dreal& cx,const mdpack::dreal& cy, const mdpack::dreal& cz, const mdpack::dreal& Lx,
	   const mdpack::dreal& Ly,const mdpack::dreal& Lz){
	     size =0;  x=cx; y=cy; z=cz; lx=Lx; ly =Ly; lz=Lz;
	     leftx=x-0.5*lx; rightx=x+0.5*lx; lefty=y-0.5*ly; righty=y+0.5*ly; 
	     leftz=z -0.5*lz; rightz=z+ 0.5*lz;
	     nNeighbours =0;
	  }
             
//       ~cell(){ if(atomIds !=NULL) delete [] atomIds; if(exitAtomIds != NULL) delete exitAtomIds;
//                if(enterAtomIds != NULL) delete [] enterAtomIds; }
               
      void addAtom(const int& id){ atomIds.push_back(id);}
      
      void setNeighbours(const std::list<int>& L){
	neighbours = L;
	nNeighbours = neighbours.size();
      }
      
      int getSize(){ return atomIds.size();}
      
      bool isIn(mdpack::vector3& pos){
	if( (pos.x>= leftx && pos.x < rightx) && (pos.y>= lefty && pos.y < righty)
	    && (pos.z>= leftz && pos.z < rightz))
	  return true;
	else
	  return false;
      }
	
      
      // update the list of atomIds and exitAtomIds 
      void update(mdpack::atom * atoms){
	exitAtomIds.clear();
	std::list<int>::iterator it;
	for(it=atomIds.begin(); it!=atomIds.end(); ++it){
	  if(!isIn(atoms[*it-1].r)){
	    atomIds.erase(it);
	    exitAtomIds.push_back(*it);
	  }
	}
	size = atomIds.size();
      }
      
      // update the list of atomIds by checking the list of exitAtomIds in the
      // neighbouring cells
      void update(mdpack::cell * cells, mdpack::atom * atoms){
	std::list<int>::iterator it1;
	std::list<int>::iterator it2;
	for(it1=neighbours.begin(); it1!=neighbours.end(); ++it1){
	  for(it2=cells[*it1].exitAtomIds.begin(); it2!=cells[*it1].exitAtomIds.end(); ++it2){
	    if(isIn(atoms[*it2-1].r)){
	      atomIds.push_back(*it2);
	      cells[*it1].exitAtomIds.erase(it2);
	    }
	  }
	}
	size= atomIds.size();
      }
	 
	  
	   
    
  };
  
}; //end of namespace mdpack
    



#endif  //end of CELL_H


