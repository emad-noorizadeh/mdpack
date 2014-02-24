#include "vector3.h"



mdpack::vector3::vector3(const mdpack::vector3& v)
{
  x=v.x; y=v.y; z=v.z;
}

mdpack::vector3& mdpack::vector3::operator=(const mdpack::vector3& v){
  if(this != &v){
   x=v.x; y=v.y; z=v.z;
  }
  return (*this);
}

const mdpack::vector3 mdpack::vector3::operator+(const vector3& v) const{
  //return (mdpack::vector3(*this) += v);
  return (mdpack::vector3(x+v.x, y+v.y, z+v.z));
}


const mdpack::vector3 mdpack::vector3::operator+(const mdpack::dreal& alpha) const{
  //return (mdpack::vector3(*this) += alpha);
  return (mdpack::vector3(x+alpha, y+alpha, z+alpha));
}

const mdpack::vector3 mdpack::vector3::operator-(const vector3& v) const{
  return mdpack::vector3( x-v.x, y-v.y, z-v.z);
}

const mdpack::vector3 mdpack::vector3::operator-(const mdpack::dreal& alpha) const{
  return mdpack::vector3( x-alpha, y-alpha, z-alpha);
}


const mdpack::vector3 mdpack::vector3::operator/(const mdpack::dreal& alpha) const{
  mdpack::dreal temp= 1.0/alpha;
  //return (mdpack::vector3(*this) /= alpha);
  return mdpack::vector3( temp*x, temp*y, temp*z);
  
}

const mdpack::vector3 mdpack::vector3::operator*(const mdpack::dreal& alpha) const{
  return mdpack::vector3( alpha*x, alpha*y, alpha*z);
  //return (mdpack::vector3(*this) *= alpha);
}


const mdpack::dreal mdpack::vector3::norm(){
  return sqrt(x*x + y*y + z*z);
}

const mdpack::dreal mdpack::vector3::dotProd(const mdpack::vector3& v){
  return (x*x + y*y + z*z);
}


const mdpack::vector3 mdpack::vector3::crossProd(const mdpack::vector3& v){
  return mdpack::vector3( (y*v.z-z*v.y), (z*v.x-x*v.z), (x*v.y-y*v.x) );
//   mdpack::vector3 temp;
//   temp.x=y*v.z-z*v.y;
//   temp.y=z*v.x-x*v.z;
//   temp.z=x*v.y-y*v.x;
//   return temp;
}





 




