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

#ifndef UTILITIES_H
#define UTILITIES_H

#include <stdio.h>
#include <stdlib.h>
#include <cctype>
#include <math.h>
#include <iostream>
#include <fstream>
#include <sstream>
//#include <limits>

#include <cstring>

#include "definitions.h"


// ========== Max ===============
template<class TYPE>
TYPE max(TYPE x, TYPE y)
{
 return ((x>=y)?x:y);
}

// ========== Min ===============
template<class TYPE>
TYPE min(TYPE x, TYPE y)
{
 return ((x<=y)?x:y);
}

// ========== square ===============
template<class TYPE>
TYPE sqr(TYPE x)
{
  return x*x;
}
// ========== Abs ===============
template<class TYPE>
TYPE abs(TYPE x)
{
  return ((x>=0)?x:-x);
}

// ========== sig ===============
// template<class TYPE>
// double sig(TYPE x)
// {
//   return ((x>=0)?1.0:-1.0);
// }
#define SIG(x) ((x>=0)?1.0:-1.0)

// ========== isnan ===============
// template<typename T>
// bool isNaN(T value)
// {
//   return value != value;
// 
// }

// ========== isinf ===============
// template<typename T>
// bool isinf(T value)
// {
//   return  ( (std::numeric_limits<T>::has_infinity) &&
//           (value == std::numeric_limits<T>::infinity()) );
// }

// // function to write in binary 
// template<typename T>
// std::ostream& binary_write(std::ostream* outStream, const T& value){
//   return outStream->write(reinterpret_cast<const char*>(&value), sizeof(T));
// }
// 
// // function to read in binary 
// template<typename T>
// std::istream & binary_read(std::istream* stream, T& value){
//   return stream->read(reinterpret_cast<char*>(&value), sizeof(T));
// }



inline const mdpack::dreal pow2(const mdpack::dreal& d){
  return (d*d);
}
inline const mdpack::dreal pow3(const mdpack::dreal& d){
  return (d*d*d);
}
inline const mdpack::dreal pow6(const mdpack::dreal& d){
  return pow2(d*d*d);
}
inline const mdpack::dreal pow7(const mdpack::dreal& d){
  return pow2(d*d*d)*d;
}
inline const mdpack::dreal pow12(const mdpack::dreal& d){
  return (pow3(d*d*d*d));
}
inline const mdpack::dreal pow13(const mdpack::dreal& d){
  return (pow12(d)*d);
}

  


  
inline void terminate(const char * msg){
  std::cerr <<"\n"<< msg << ".....exiting.."<<"\n"<<std::endl;
  exit(1);
}

inline int countParts (const std::string& str){
   int numwords=0;
   int size= (int) (str.length());
   int i=0;
   int count =0;
   if(size >0 && str.at(0) =='!'){
     return 0;
   }
   else{
     while(i<size-1 && isspace(str.at(i)) && str.at(i) !='!' )
      i++;
     for (; i<size; i++){
       if(str.at(i) =='!')
	 return numwords;
       else if (isspace(str[i])){
	 numwords +=1;
	 while(i<size-1 && isspace(str.at(i+1)))
	   i++;
	}
	else if(i== size -1)
	  numwords +=1;
	else
	  count +=1;
     }
     if(count == size && str[0] !='\0')
       return 1;
     else
       return numwords;
  }
}

int countWords(const std::string& line);
//  int countWords(const std::string& line){
//   std::istringstream iss(line,std::istringstream::in);
//   char temp[100];
//   int count =0;
//   while( iss >> temp){
//     count +=1;
//   }
//   return count;
// }
  
bool getWord(const std::string& line,char * str, const int& n);
// bool getWord(const std::string& line,char * str, const int& n){
//   std::istringstream iss(line,std::istringstream::in);
//   //char temp [100];
//   int k=0;
//   for(int i=0; i<n; i++){
//     if(iss >> str)
//       k +=1;
//   }
//   if(k == n)
//     return true;
//   else
//     return false;
// }




/** \c getWord: 
  * given a string (line) it extract the words from it and returns a 1D array of words:
  * \c inputs:
  * std::string line 
  * int nwords: will be replaced with the number of words in the line
  * char *words[maxWords]: must be initilized array of char pointer. it will be changed so
  * that word[i] will hold the i^th word in the line
  * const int maxWords: it is the size of words
  * 
 */
void getWords(const std::string& line, int& nwords, char *words[], const int& maxWords);
// void getWords(const std::string& line, int& nwords, char *words[], const int& maxWords){
//   std::istringstream iss(line,std::istringstream::in);
//   for(int i=0; i<maxWords; i++){
//     if(words[i]) delete [] words[i];
//   }
//   char temp [100];
//   int k=0;
//   int i;
//   while(iss >> temp){
//     i=0;
//     while(temp[i] !='\0') i +=1;
//     words[k] = new char[i+1];
//     strcpy(words[k], temp);
//     k +=1;
//   }
//   
//   nwords =k;
// }
   


inline void readBond(const std::string& line, char * str1,
		     char * str2, mdpack::dreal& d1, 
		     mdpack::dreal& d2){
  std::istringstream iss(line,std::istringstream::in);
  iss >> str1 >> str2 >> d1 >> d2;
}

inline void readAngle( const std::string& line, char * str1,
		       char * str2, char * str3,
		       mdpack::dreal& d1,mdpack::dreal& d2){
  std::istringstream iss(line,std::istringstream::in);
  iss >> str1 >> str2 >> str3 >> d1 >> d2;
}
  




inline void readAngle( const std::string& line, char * str1,
		       char * str2, char * str3,
		       mdpack::dreal& d1,mdpack::dreal& d2,
		       mdpack::dreal& d3,mdpack::dreal& d4){
  
  std::istringstream iss(line,std::istringstream::in);
  iss >> str1 >> str2 >> str3 >> d1 >> d2 >> d3 >> d4;
}

inline void readDihedral(const std::string& line, char * str1, char * str2, 
		  char * str3, char * str4,  mdpack::dreal& k_chi,
		  mdpack::dreal& multip, mdpack::dreal& delta){
  
  std::istringstream iss(line,std::istringstream::in);
  iss >> str1 >> str2 >> str3 >> str4>> k_chi >> multip >> delta;
}


inline void readImproper(const std::string& line, char * str1, char * str2, char * str3, 
		  char * str4, mdpack::dreal& k_psi,mdpack::dreal& unused, mdpack::dreal& psi_0){
  
  std::istringstream iss(line,std::istringstream::in);
  iss >> str1 >> str2 >> str3 >> str4>> k_psi >> unused >> psi_0;
}

inline void readNonBonded(const std::string& line, char * str1, mdpack::dreal& unused, 
		   mdpack::dreal& epsilon, mdpack::dreal& rmin){
  
  std::istringstream iss(line,std::istringstream::in);
  iss >> str1 >> unused >> epsilon >> rmin ;
}
  




#endif  //UTILITIES_H