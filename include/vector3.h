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


#ifndef VECTOR3_H
#define VECTOR3_H

#include <stdio.h>
//#include <stdlib.h>
#include <math.h>
#include <iostream>

//#include "utilities.h"
#include "constants.h"
#include "definitions.h"


// inline float Q_rsqrt( float number )
// {
//   long i;
//   float x2, y;
//   const float threehalfs = 1.5F;
//  
//   x2 = number * 0.5F;
//   y  = number;
//   i  = * ( long * ) &y;                       // evil floating point bit level hacking
//   i  = 0x5f3759df - ( i >> 1 );               // what the fuck?
//   y  = * ( float * ) &i;
//   y  = y * ( threehalfs - ( x2 * y * y ) );   // 1st iteration
//   y  = y * ( threehalfs - ( x2 * y * y ) );   // 2nd iteration, this can be removed
// //y  = y * ( threehalfs - ( x2 * y * y ) );
//  
//   return y;
// }
// 
// //sqrt(1/x) = 2^((1/2)*log2(1/x)) = 2^(-(1/2)*log2(x))
// inline const double invSqrt(const double x){
//   //it is 7 times slower than 1/sqrt
//  // return pow(2.0, -(0.5)*log2(x));
//   return 1.0/sqrt(x);
// }

namespace mdpack{
  
  class vector3 {
  public:
    mdpack::dreal x,y,z;
    vector3(){x=0.0; y=0.0; z=0.0;}
    vector3(const mdpack::dreal& x1, const mdpack::dreal& y1,const  mdpack::dreal& z1){ x=x1; y=y1; z=z1;}
      
    // copy constructor 
    vector3(const vector3& v);
      
    // constructor with array argument
    vector3(const mdpack::dreal a[]);
    
    ~vector3( void) { }
      
    void set(mdpack::dreal x1, mdpack::dreal y1, mdpack::dreal z1){x=x1; y=y1; z=z1;}
      
      
    vector3& operator=(const vector3& v);
      
      
    inline vector3& operator+=(const vector3& v){
      x += v.x; y += v.y; z += v.z; return (*this);
    }
    inline vector3& operator+=(const mdpack::dreal& alpha){
      x +=alpha; y +=alpha; z +=alpha; return (*this);
    }
	
    const vector3 operator+(const vector3& v) const;
    const vector3 operator+(const mdpack::dreal& alpha) const;
      
    vector3& operator-=(const vector3& v){
      x -= v.x; y -= v.y; z -= v.z; return (*this);
    }
	
    vector3& operator-=(const mdpack::dreal& alpha){
      x -=alpha; y -=alpha; z -=alpha; return (*this);
    }
	
    const vector3 operator-(const vector3& v)const;
    const vector3 operator-(const mdpack::dreal& alpha) const;
      
    vector3& operator*=(const mdpack::dreal& alpha){
      x *=alpha; y *=alpha; z *=alpha; return (*this);
    }

    const vector3 operator*(const mdpack::dreal& alpha)const;
    inline vector3& operator/=(const mdpack::dreal& alpha){
      mdpack::dreal temp = 1.0/alpha;
      x *=temp; y *=temp; z *=temp; return (*this);
    }
      

	
    const vector3 operator/(const mdpack::dreal& alpha) const;
      
    const mdpack::dreal dotProd(const vector3& v);
    const vector3 crossProd(const vector3& v);
      
      
    void clear(){x=0.0; y=0.0; z=0.0;}
    const mdpack::dreal sum(){ return (x+y+z);}
      
    const mdpack::dreal norm();
    void print(){ printf("x=%f, y=%f, z=%f \n",x,y,z);}
      
    const mdpack::dreal min(){ 
      mdpack::dreal ans; 
      if(x<y && x<z) ans=x;
      else if(y<x && y<z) ans=y;
      else ans =z;
      return ans;
    }
    const mdpack::dreal max(){ 
      mdpack::dreal ans; 
      if(x>y && x>z) ans=x;
      else if(y>x && y>z) ans=y;
      else ans =z;
      return ans;
    }
};// vector

}; // mdpack



inline const mdpack::dreal dotProd(const mdpack::vector3& v1, const mdpack::vector3& v2){
  return (v1.x*v2.x + v1.y*v2.y + v1.z*v2.z);
}

inline const mdpack::dreal norm(mdpack::vector3& v){
  return sqrt(v.x*v.x + v.y*v.y + v.z*v.z);
}

inline const mdpack::dreal inorm(mdpack::vector3& v){
  return 1.0/sqrt(v.x*v.x + v.y*v.y + v.z*v.z);
}

inline const mdpack::dreal sum(mdpack::vector3& v){
  return (v.x + v.y + v.z);
}

inline const mdpack::vector3 operator*(const mdpack::dreal& alpha, 
				      const mdpack::vector3& v){
  return(v*alpha);
}

inline const mdpack::vector3 crossProd( const mdpack::vector3& v1,
				       const mdpack::vector3& v2){
  return mdpack::vector3( (v1.y*v2.z-v1.z*v2.y), (v1.z*v2.x-v1.x*v2.z), (v1.x*v2.y-v1.y*v2.x) );
 /* 
  mdpack::vector3 temp;
  temp.x=v1.y*v2.z-v1.z*v2.y;
  temp.y=v1.z*v2.x-v1.x*v2.z;
  temp.z=v1.x*v2.y-v1.y*v2.x;
  return temp;*/
}


/** vector max(v): returns the maximum element of v
 */
inline const mdpack::dreal max(mdpack::vector3& v){
  return v.max();
}

/** vector min(v): returns the minimum element of v
*/
inline const mdpack::dreal min(mdpack::vector3& v){
  return v.min();
}
      


#endif
