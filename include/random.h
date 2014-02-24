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
    
    this file uses modified third party codes for Mersenne twister, with copyright
    Copyright (c) 2006, 2007 by Mutsuo Saito, Makoto Matsumoto and Hiroshima University.
    Copyright (c) 2008 by Agner Fog.
  
 */

/** 
 * This file provide functions for generating normal random variables based on 
 * uniform generators that are based on Mersenne twister, Park and Miller with Bays-Durham
 * shuffle. In our test, randn1 was the weakest performer and randn and randn3 they were 
 * very simular but randn (Mersenne twister) has much longer period, hence we chose randn
 * as a default generator for normal random generator.
 * 
 */

#ifndef RANDOM_H
#define RANDOM_H

#include <math.h>

// Define integer types with known size: int32_t, uint32_t, int64_t, uint64_t.
// If this doesn't work then insert compiler-specific definitions here:
#if defined(__GNUC__)
  // Compilers supporting C99 or C++0x have inttypes.h defining these integer types
  #include <inttypes.h>
  #define INT64_SUPPORTED // Remove this if the compiler doesn't support 64-bit integers
#elif defined(_WIN16) || defined(__MSDOS__) || defined(_MSDOS) 
   // 16 bit systems use long int for 32 bit integer
  typedef   signed long int int32_t;
  typedef unsigned long int uint32_t;
#elif defined(_MSC_VER)
  // Microsoft have their own definition
  typedef   signed __int32  int32_t;
  typedef unsigned __int32 uint32_t;
  typedef   signed __int64  int64_t;
  typedef unsigned __int64 uint64_t;
  #define INT64_SUPPORTED // Remove this if the compiler doesn't support 64-bit integers
#else
  // This works with most compilers
  typedef signed int          int32_t;
  typedef unsigned int       uint32_t;
  typedef long long           int64_t;
  typedef unsigned long long uint64_t;
  #define INT64_SUPPORTED // Remove this if the compiler doesn't support 64-bit integers
#endif



// //#include "randomc.h"
// #if defined(__GNUC__)
//   // Compilers supporting C99 or C++0x have inttypes.h defining these integer types
//   #include <inttypes.h>
//   #define INT64_SUPPORTED // Remove this if the compiler doesn't support 64-bit integers
// #else
//  typedef signed int          int32_t;
//  typedef unsigned int       uint32_t;
//  typedef long long           int64_t;
//  typedef unsigned long long uint64_t;
// #define INT64_SUPPORTED // Remove this if the compiler doesn't support 64-bit integers
//  
 // constants for type MT19937:
#define MERS_N   624
#define MERS_M   397
#define MERS_R   31
#define MERS_U   11
#define MERS_S   7
#define MERS_T   15
#define MERS_L   18
#define MERS_A   0x9908B0DF
#define MERS_B   0x9D2C5680
#define MERS_C   0xEFC60000
 
 
/**======================================================================*/
 /** class CRandomMersenne:
  * Random number generator of type Mersenne twister.
  */
 /*
 Copyright (c) 2006, 2007 by Mutsuo Saito, Makoto Matsumoto and Hiroshima University.
  Copyright (c) 2008 by Agner Fog.
 */
  

class CRandomMersenne {  
  public:
   CRandomMersenne() {         // Constructor
      RandomInit(0); LastInterval = 0;}
   CRandomMersenne(int seed) {         // Constructor
      RandomInit(seed); LastInterval = 0;}
   void RandomInit(int seed);          // Re-seed
   int IRandom (int min, int max);     // Output random integer
   double Random();                    // Output random float
   uint32_t BRandom();                 // Output random bits
private:
   void Init0(int seed);               // Basic initialization procedure
   uint32_t mt[MERS_N];                // State vector
   int mti;                            // Index into mt
   uint32_t LastInterval;              // Last interval length for IRandomX
   uint32_t RLimit;                    // Rejection limit used by IRandomX
   
};
/**==================================================================*/
  
/**=================Random Number Generator=========================*/
inline double ran3(int &seed)
// Random number generator (from Numerical Recipes)
// seed must be negative integer
{ const int MBIG = 1000000000;
  const int MSEED = 161803398;
  const int MZ = 0;
  const double FAC = 1.0/MBIG;
  static int inext,inextp;
  static long ma[56];
  static int iff=0;
  long mj,mk;
  int i,ii,k;

  if (seed < 0 || iff == 0) 
  { iff = 1;
    mj = MSEED-(seed < 0 ? -seed : seed);
    mj %= MBIG;   ma[55]=mj;   mk=1;
    for (i=1;i<=54;i++) 
    { ii=(21*i) % 55;   ma[ii]=mk;   mk=mj-mk;
      if (mk < MZ) mk += MBIG;   mj=ma[ii]; }
    for (k=1;k<=4;k++)
      for (i=1;i<=55;i++) 
      { ma[i] -= ma[1+(i+30) % 55];
        if (ma[i] < MZ) ma[i] += MBIG; }
    inext=0;   inextp=31;   seed=1; }
  if (++inext == 56) inext=1;   
  if (++inextp == 56) inextp=1;
  mj = ma[inext] - ma[inextp];
  if (mj < MZ) mj += MBIG;
  ma[inext] = mj;

  return(mj*FAC); 
}
/**=============================================================================*/

/**=======================from Numerical Recipes================================*/
/** random number generator of Park and Miller with Bays-Durham shuffle and 
 * added safegaurds. Returns a uniform random deviate between 0.0 and 1.0
 * (exclusive of the the endpoint values). Call with seed a negative integer
 * to initialize; thereafter do not alter seed between succesive deviates in a 
 * sequence. RNMX should approximate the largest floating value that is less that 1.
 */
inline double ran1(int &seed)
{
const int IA=16007, IM=2147483647, IQ=127773, IR=2836, NTAB=32;
const int NDIV=(1+(IM-1)/NTAB);
const double EPS=3.0e-16, AM=1.0/IM, RNMX = (1.0-EPS);
static int iy =0;
static int iv [NTAB];
int j,k;
double temp;
if(seed <=0 || !iy){
  if(-seed < 1) 
     seed=1;
  else 
     seed = -seed;
  for(j= NTAB +7; j>=0;j--){
     k=seed/IQ;
     seed= IA*(seed-k*IQ)-IR*k;
     if(seed < 0)
       seed +=IM;
     if(j < NTAB)
       iv[j]= seed;
   }
   iy=iv[0];
}
k=seed/IQ;
seed=IA*(seed-k*IQ)-IR*k;
if(seed<0)
  seed +=IM;
j=iy/NDIV;
iy=iv[j];
iv[j]= seed;
if((temp=AM*iy) > RNMX)
  return RNMX;
else
  return temp;
}  
/**=============================================================================*/


/**=============================================================================*/
  /**
  * class randn, randn1, randn3: provide normal random number generator using box-muller method,
  * based on uniform random number generators using Mersenne twister, ran1, ran3.
  * 
  * 
  */

//============== randn ==============
class randn {
private:
  int seed;
  int count;
  double gaussianVar1;
  double gaussianVar2;
  double v1;
  double v2;
  //CRandomSFMT RanGen(1);       // make instance of random number generator
  CRandomMersenne RanGen;
  //temp t(0);
 
public: 
  
  randn(){seed =0; count =0; gaussianVar1=0.0; gaussianVar2=0.0;v1=0.0; v2=0.0;}
  randn(int s);
  const double operator()();
  const int getSeed(void){ return seed;}
  void setSeed(int s){seed =s; RanGen.RandomInit(seed);}
    
};

//=============== randn1  ====================
class randn1 {
protected:
  int seed;
  int count;
  double gaussianVar1;
  double gaussianVar2;
  double v1;
  double v2;
 
public: 
  
  randn1(){seed =-1; count =0; gaussianVar1=0.0; gaussianVar2=0.0;v1=0.0; v2=0.0;}
  randn1(int s){ if(s <0) seed =s; else if(s >0) seed =-s; else seed =-1;}
  virtual const double operator()();
  const int getSeed(void){ return seed;}
  void setSeed(int s){ seed =s;}
    
};

//=============== randn3  ====================
class randn3 : public randn1{
 
public: 
  randn3(int s) : randn1(s){}
  virtual const double operator()();
  
};
/**==============================================================================*/





#endif  //end of RANDOM_H

