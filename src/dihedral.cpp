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

#include "dihedral.h"

mdpack::dihedral::dihedral(){
  atom1=0;
  atom2=0;
  atom3=0;
  atom4=0;
  strcpy(type1, mdpack::emptystr);
  strcpy(type2, mdpack::emptystr);
  strcpy(type3, mdpack::emptystr);
  strcpy(type4, mdpack::emptystr);
  isSet= false;
  k=0.0;
  n=0.0;
  delta=0.0;
}

mdpack::dihedral::dihedral(int a1, int a2,int a3, int a4, const char * str1, const char * str2,
	    const char * str3, const char * str4, const mdpack::dreal& k_phi, 
	    const mdpack::dreal& mult,const mdpack::dreal& d){
  atom1=a1;
  atom2=a2;
  atom3=a3;
  atom4=a4;
  strcpy(type1, str1);
  strcpy(type2, str2);
  strcpy(type3, str3);
  strcpy(type4, str4);
  isSet= true;
  k=k_phi;
  n=mult;
  delta=d;
}

mdpack::dihedral::dihedral(const dihedral& d){
  atom1=d.atom1;
  atom2=d.atom2;
  atom3=d.atom3;
  atom4=d.atom4;
  strcpy(type1, d.type1);
  strcpy(type2, d.type2);
  strcpy(type3, d.type3);
  strcpy(type4, d.type4);
  isSet= d.isSet;
  k=d.k;
  n=d.n;
  delta=d.delta;
}

// assignment operator
mdpack::dihedral& mdpack::dihedral::operator=(const mdpack::dihedral& d){
  if(this != &d){
    atom1=d.atom1;
    atom2=d.atom2;
    atom3=d.atom3;
    atom4=d.atom4;
    strcpy(type1, d.type1);
    strcpy(type2, d.type2);
    strcpy(type3, d.type3);
    strcpy(type4, d.type4);
    isSet= d.isSet;
    k=d.k;
    n=d.n;
    delta=d.delta;
  }
  return (*this);
}

void mdpack::dihedral::set(int a1, int a2, int a3, int a4, const char * str1, 
	       const char * str2, const char * str3, const char * str4){
    atom1=a1;
    atom2=a2;
    atom3=a3;
    atom4=a4;
    strcpy(type1, str1);
    strcpy(type2, str2);
    strcpy(type3, str3);
    strcpy(type4, str4);
}

void mdpack::dihedral::setPrm(const mdpack::dreal& k_phi,const mdpack::dreal& mult,
		  const mdpack::dreal& d){
  isSet = true;
  k= k_phi;
  n= mult;
  delta = d;
}

bool mdpack::dihedral::compare(const char * str1, const char * str2, const char * str3,
		   const char * str4){
  if( strcmp(type1, str1)==0 && strcmp(type2, str2)==0 &&
     strcmp(type3, str3)==0 && strcmp(type4, str4)==0    )
    return true;
  else if( strcmp(type1, str4)==0 && strcmp(type2, str3)==0 &&
           strcmp(type3, str2)==0 && strcmp(type4, str1)==0   )
    return true;
  else if( strcmp(str1,"X")==0 && strcmp(type2, str2)==0 &&
           strcmp(type3, str3)==0 && strcmp(str4,"X")==0    )
    return true;
  else if( strcmp(str1,"X")==0 && strcmp(type2, str3)==0 &&
           strcmp(type3, str2)==0 && strcmp(str4,"X")==0    )
    return true;
  else
    return false;
}

bool mdpack::dihedral::operator==(const dihedral& d){
  return (this->compare(d.type1, d.type2, d.type3, d.type4));
}

bool mdpack::dihedral::operator!=(const dihedral& d){
  return !(*this ==d);
}



// bool mdpack::dihedralParams::operator==(const mdpack::dihedralParams& d){
//   if( strcmp(atype1, d.atype1)==0 && strcmp(atype2, d.atype2)==0 &&
//      strcmp(atype3, d.atype3)==0 && strcmp(atype4, d.atype4)==0    )
//     return true;
//   else if( strcmp(atype1, d.atype4)==0 && strcmp(atype2, d.atype3)==0 &&
//            strcmp(atype3, d.atype2)==0 && strcmp(atype4, d.atype1)==0   )
//     return true;
//   else if( strcmp(d.atype1,"X")==0 && strcmp(atype2, d.atype2)==0 &&
//            strcmp(atype3, d.atype3)==0 && strcmp(d.atype4,"X")==0    )
//     return true;
//   else if( strcmp(d.atype1,"X")==0 && strcmp(atype2, d.atype3)==0 &&
//            strcmp(atype3, d.atype2)==0 && strcmp(d.atype4,"X")==0    )
//     return true;
//   else
//     return false;
// }
//     
// 
// bool mdpack::dihedralParams::operator!=(const mdpack::dihedralParams& d){
//   return !(*this==d);
// }



