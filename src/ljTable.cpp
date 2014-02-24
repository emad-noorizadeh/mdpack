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

#include "ljTable.h"


mdpack::ljTable::ljTable(const mdpack::ljTable& LJ){
  this->inc = LJ.inc; 
  this->A   = LJ.A; 
  this->B   = LJ.B; 
  this->id1 = LJ.id1; 
  this->id2 = LJ.id2;
}

mdpack::ljTable& mdpack::ljTable::operator=(const mdpack::ljTable& LJ){ 
  this->inc = LJ.inc; 
  this->A   = LJ.A; 
  this->B   = LJ.B; 
  this->id1 = LJ.id1; 
  this->id2 = LJ.id2; 
  return(*this);
  
}   

// print information into std output
void mdpack::ljTable::print(){ 
  printf("nonBondedPair  atom1 = %-7d  atom2 = %-7d  A = %-6.3f  B = %-6.3f  inc = %d\n",
	       id1, id2,A,B, inc );
}