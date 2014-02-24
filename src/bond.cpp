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

#include "bond.h"
      
mdpack::bond::bond(){
  atom1=0;
  atom2=0;
  strcpy(type1, mdpack::emptystr);
  strcpy(type2, mdpack::emptystr);
  isSet = false;
  k=0.0;
  r0=0.0;
}
      
mdpack::bond::bond(int a1, int a2,const char * str1, const char * str2, mdpack::dreal& kapa,
                  mdpack::dreal& b0){
  atom1=a1;
  atom2=a2;
  strcpy(type1, str1);
  strcpy(type2, str2);
  isSet = true;
  k=kapa;
  r0=b0;
}

mdpack::bond::bond(const mdpack::bond& b){
  atom1=b.atom1;
  atom2=b.atom2;
  strcpy(type1, b.type1);
  strcpy(type2, b.type2);
  isSet = b.isSet;
  k= b.k;
  r0=b.r0;
}

mdpack::bond& mdpack::bond::operator=(const mdpack::bond& b){
  if(this !=&b){
    atom1=b.atom1;
    atom2=b.atom2;
    strcpy(type1, b.type1);
    strcpy(type2, b.type2);
    isSet = b.isSet;
    k=b.k;
    r0=b.r0;
  }
  return (*this);
}
    
void mdpack::bond::set(int a1, int a2, const char * str1, const char * str2 ){
  atom1=a1;
  atom2=a2;
  strcpy(type1, str1);
  strcpy(type2, str2);
}
  
      
bool mdpack::bond::operator==(const mdpack::bond& b){
  if((atom1 == b.atom1 && atom2== b.atom2) || (atom1 == b.atom2 && atom2== b.atom1))
    return true;
  else
    return false;
}

bool mdpack::bond::operator!=(const mdpack::bond& b){ 
  return !(*this == b);
}

bool mdpack::bond::compare(const char * str1, const char * str2){
  if((strcmp(type1, str1)==0 && strcmp(type2, str2)==0) || 
     (strcmp(type1, str2)==0 && strcmp(type2, str1)==0))
    return true;
  else
    return false;
}
  
      
      
      
      
      
      
      