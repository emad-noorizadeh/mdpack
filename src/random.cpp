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


#include "random.h"

void CRandomMersenne::Init0(int seed) {
   // Seed generator
   const uint32_t factor = 1812433253UL;
   mt[0]= seed;
   for (mti=1; mti < MERS_N; mti++) {
      mt[mti] = (factor * (mt[mti-1] ^ (mt[mti-1] >> 30)) + mti);
   }
}

void CRandomMersenne::RandomInit(int seed) {
   // Initialize and seed
   Init0(seed);

   // Randomize some more
   for (int i = 0; i < 37; i++) BRandom();
}



uint32_t CRandomMersenne::BRandom() {
   // Generate 32 random bits
   uint32_t y;

   if (mti >= MERS_N) {
      // Generate MERS_N words at one time
      const uint32_t LOWER_MASK = (1LU << MERS_R) - 1;       // Lower MERS_R bits
      const uint32_t UPPER_MASK = 0xFFFFFFFF << MERS_R;      // Upper (32 - MERS_R) bits
      static const uint32_t mag01[2] = {0, MERS_A};

      int kk;
      for (kk=0; kk < MERS_N-MERS_M; kk++) {    
         y = (mt[kk] & UPPER_MASK) | (mt[kk+1] & LOWER_MASK);
         mt[kk] = mt[kk+MERS_M] ^ (y >> 1) ^ mag01[y & 1];}

      for (; kk < MERS_N-1; kk++) {    
         y = (mt[kk] & UPPER_MASK) | (mt[kk+1] & LOWER_MASK);
         mt[kk] = mt[kk+(MERS_M-MERS_N)] ^ (y >> 1) ^ mag01[y & 1];}      

      y = (mt[MERS_N-1] & UPPER_MASK) | (mt[0] & LOWER_MASK);
      mt[MERS_N-1] = mt[MERS_M-1] ^ (y >> 1) ^ mag01[y & 1];
      mti = 0;
   }
   y = mt[mti++];

   // Tempering (May be omitted):
   y ^=  y >> MERS_U;
   y ^= (y << MERS_S) & MERS_B;
   y ^= (y << MERS_T) & MERS_C;
   y ^=  y >> MERS_L;

   return y;
}


double CRandomMersenne::Random() {
   // Output random float number in the interval 0 <= x < 1
   // Multiply by 2^(-32)
   return (double)BRandom() * (1./(65536.*65536.));
}


int CRandomMersenne::IRandom(int min, int max) {
   // Output random integer in the interval min <= x <= max
   // Relative error on frequencies < 2^-32
   if (max <= min) {
      if (max == min) return min; else return 0x80000000;
   }
   // Multiply interval with random and truncate
   int r = int((double)(uint32_t)(max - min + 1) * Random() + min); 
   if (r > max) r = max;
   return r;
}





/** Returns a normally distributed deviate with sero mean and unit variance.
 *  This uses Box-Muller algorithm
 * 
 */
randn::randn(int s){
  seed=s; 
  count =0;
  gaussianVar1=0.0; 
  gaussianVar2=0.0; 
  v1=0.0;
  v2=0.0;
  RanGen.RandomInit(seed);
}


const double randn::operator()(void){
  double s, ans;

  if (count == 0) {
      do { 
	  v1 = 2.0 * RanGen.Random() - 1;   // between -1.0 and 1.0
	  v2 = 2.0 * RanGen.Random() - 1;   // between -1.0 and 1.0
	  s = v1 * v1 + v2 * v2;
      } while (s >= 1 || s == 0);
      double multiplier = sqrt(-2 * log(s)/s);
      gaussianVar2 = v2 * multiplier;
      gaussianVar1 = v1 * multiplier;
      count = 1;
      ans= gaussianVar1;
  }
  else{
  count = 0;
  ans= gaussianVar2;
  }
  return ans;
}


const double randn1::operator()(void){
  double s, ans;

  if (count == 0) {
      do { 
	  v1 = 2.0 * ran1(seed) - 1;   // between -1.0 and 1.0
	  v2 = 2.0 * ran1(seed) - 1;   // between -1.0 and 1.0
	  s = v1 * v1 + v2 * v2;
      } while (s >= 1 || s == 0);
      double multiplier = sqrt(-2 * log(s)/s);
      gaussianVar2 = v2 * multiplier;
      gaussianVar1 = v1 * multiplier;
      count = 1;
      ans= gaussianVar1;
  }
  else{
  count = 0;
  ans= gaussianVar2;
  }
  return ans;
}

const double randn3::operator()(void){
  double s, ans;

  if (count == 0) {
      do { 
	  v1 = 2.0 * ran3(seed) - 1;   // between -1.0 and 1.0
	  v2 = 2.0 * ran3(seed) - 1;   // between -1.0 and 1.0
	  s = v1 * v1 + v2 * v2;
      } while (s >= 1 || s == 0);
      double multiplier = sqrt(-2 * log(s)/s);
      gaussianVar2 = v2 * multiplier;
      gaussianVar1 = v1 * multiplier;
      count = 1;
      ans= gaussianVar1;
  }
  else{
  count = 0;
  ans= gaussianVar2;
  }
  return ans;
}
  
  
