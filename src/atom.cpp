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

#include "atom.h"

mdpack::atom::atom(){
  id=0;
  strcpy(segName,mdpack::emptystr);
  resId= 0;
  strcpy(resName,mdpack::emptystr );
  strcpy(name, mdpack::emptystr);
  strcpy(type, mdpack::emptystr);
  q= 1.0;
  m= 1.0;
  im =1.0/m;
  isSet= false;
  eps=0.0;
  halfSigma =0.0;
  tempFactor =0.0;
  occupancy =1.0;
  r.clear();
}

mdpack::atom::atom(int aid, const char * sname, int rId, const char * rname, const char * n,
	   const char * t, const mdpack::dreal& charge, const mdpack::dreal& mass,
	   mdpack::dreal& e,mdpack::dreal& sig, mdpack::vector3& coord, mdpack::dreal& oc,
		  mdpack::dreal& tempf){
  id= aid;
  strcpy(segName, sname);
  resId= rId;
  strcpy(resName, rname);
  strcpy(name, n);
  strcpy(type, t);
  q= charge;
  m= mass;
  im = 1.0/mass;
  isSet= true;
  eps=e;
  halfSigma =sig;
  r = coord;
  occupancy = oc;
  tempFactor = tempf;
}

mdpack::atom::atom(const mdpack::atom& a){
  id= a.id;
  strcpy(segName, a.segName);
  resId= a.resId;
  strcpy(resName, a.resName);
  strcpy(name, a.name);
  strcpy(type, a.type);
  q= a.q;
  m= a.m;
  im = a.im;
  r = a.r;
  eps = a.eps;
  halfSigma = a.halfSigma;
  isSet = a.isSet;
  tempFactor = a.tempFactor;
  occupancy = a.occupancy;
}
      
mdpack::atom& mdpack::atom::operator=(const mdpack::atom& a){
  if(this != &a){
    id= a.id;
    strcpy(segName, a.segName);
    resId= a.resId;
    strcpy(resName, a.resName);
    strcpy(name, a.name);
    strcpy(type, a.type);
    q= a.q;
    m= a.m;
    im = a.im;
    r = a.r;
    eps = a.eps;
    halfSigma = a.halfSigma;
    isSet = a.isSet;
    tempFactor = a.tempFactor;
    occupancy = a.occupancy;
  }
  return (*this);
}
  
void mdpack::atom::set(int aid, const char * sname, int rId, const char * rname, const char * n,
	   const char * t, const mdpack::dreal& charge, const mdpack::dreal& mass){
  id= aid;
  strcpy(segName, sname);
  resId= rId;
  strcpy(resName, rname);
  strcpy(name, n);
  strcpy(type, t);
  q= charge;
  m= mass;
  im = 1.0/mass;
}

bool mdpack::atom::operator==(const atom& a){
  if(id==a.id && resId == a.resId && strcmp(name, a.name)==0 &&
     strcmp(type, a.type)==0)
    return true;
  else
    return false;
}
      
bool mdpack::atom::operator!=(const atom& a){
  return !(*this == a);
}

bool mdpack::atom::isIn1_2(const int& k){
  std::list<int>::iterator it = list1_2.begin();
  bool out = false;
  while(it!=list1_2.end() && (*it) <= k){
    if((*it) == k)
      out = true;
    ++it;
  }
  return out;
}

bool mdpack::atom::isIn1_3(const int& k){
  std::list<int>::iterator it = list1_3.begin();
  bool out = false;
  while(it!=list1_3.end() && (*it) <= k){
    if((*it) == k)
      out = true;
    ++it;
  }
  return out;
}

bool mdpack::atom::isIn1_4(const int& k){
  std::list<int>::iterator it = list1_4.begin();
  bool out = false;
  while(it!=list1_4.end() && (*it) <= k){
    if((*it) == k)
      out = true;
    ++it;
  }
  return out;
}
 
    





  
  
