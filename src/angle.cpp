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

#include "angle.h"
      
mdpack::angle::angle(){
  atom1=0;
  atom2=0;
  atom3=0;
  strcpy(type1, mdpack::emptystr);
  strcpy(type2, mdpack::emptystr);
  strcpy(type3, mdpack::emptystr);
  isSet= false;
  ktheta=0.0;
  theta0=0.0;
  kub=0.0;
  rub=0.0;
}

mdpack::angle::angle(int a1, int a2,int a3, const char * str1, const char * str2,
	    const char * str3, const mdpack::dreal& k_theta, const mdpack::dreal& theta_0,
	    const mdpack::dreal& k_ub, const mdpack::dreal& r_ub){
  atom1=a1;
  atom2=a2;
  atom3=a3;
  strcpy(type1,str1);
  strcpy(type2,str2);
  strcpy(type3,str3);
  isSet= true;
  ktheta=k_theta;
  theta0=theta_0;
  kub=k_ub;
  rub=r_ub;
}

mdpack::angle::angle(const mdpack::angle& a){
  atom1=a.atom1;
  atom2=a.atom2;
  atom3=a.atom3;
  strcpy(type1,a.type1);
  strcpy(type2,a.type2);
  strcpy(type3,a.type3);
  isSet= a.isSet;
  ktheta= a.ktheta;
  theta0= a.theta0;
  kub=a.kub;
  rub=a.rub;
}
  
      
mdpack::angle& mdpack::angle::operator=(const mdpack::angle& a){
  if(this !=&a){
    atom1=a.atom1;
    atom2=a.atom2;
    atom3=a.atom3;
    strcpy(type1,a.type1);
    strcpy(type2,a.type2);
    strcpy(type3,a.type3);
    isSet= a.isSet;
    ktheta= a.ktheta;
    theta0= a.theta0;
    kub=a.kub;
    rub=a.rub;
  }
  return (*this);
}
    
      
void mdpack::angle::set(int a1, int a2, int a3, const char * str1,
			const char * str2, const char * str3){
  atom1=a1;
  atom2=a2;
  atom3=a3;
  strcpy(type1, str1);
  strcpy(type2, str2);
  strcpy(type3, str3);
}

void mdpack::angle::setPrm(const mdpack::dreal& k_theta, const mdpack::dreal& theta_0,
			   const mdpack::dreal& k_ub, const mdpack::dreal& r_ub){
  isSet = true;
  ktheta = k_theta;
  theta0= theta_0;
  kub= k_ub;
  rub= r_ub;
}
      
bool mdpack::angle::compare(const char * str1, const char * str2, const char * str3){
  if( ((strcmp(type1, str1)==0 && strcmp(type2, str2)==0) &&
      (strcmp(type3, str3)==0)) ||
      ((strcmp(type1, str3)==0 && strcmp(type2, str2)==0) &&
      (strcmp(type3, str1)==0)) )
    return true;
  else
    return false;
}
  
bool mdpack::angle::operator==(const mdpack::angle& a){
  return this->compare(a.type1, a.type2, a.type3);
}
  
bool mdpack::angle::operator!=(const mdpack::angle& a){
  return !(*this==a);
}

  
  






