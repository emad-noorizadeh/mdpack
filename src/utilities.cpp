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


#include "utilities.h"




int countWords(const std::string& line){
  std::istringstream iss(line,std::istringstream::in);
  char temp[100];
  int count =0;
  while( iss >> temp){
    count +=1;
  }
  return count;
}
  

bool getWord(const std::string& line,char * str, const int& n){
  std::istringstream iss(line,std::istringstream::in);
  //char temp [100];
  int k=0;
  for(int i=0; i<n; i++){
    if(iss >> str)
      k +=1;
  }
  if(k == n)
    return true;
  else
    return false;
}




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

void getWords(const std::string& line, int& nwords, char *words[], const int& maxWords){
  std::istringstream iss(line,std::istringstream::in);
  for(int i=0; i<maxWords; i++){
    if(words[i]) delete [] words[i];
  }
  char temp [100];
  int k=0;
  int i;
  while(iss >> temp){
    i=0;
    while(temp[i] !='\0') i +=1;
    words[k] = new char[i+1];
    strcpy(words[k], temp);
    k +=1;
  }
  
  nwords =k;
}
   


// void readBond(const std::string& line, char * str1,
// 		     char * str2, mdpack::dreal& d1, 
// 		     mdpack::dreal& d2){
//   std::istringstream iss(line,std::istringstream::in);
//   iss >> str1 >> str2 >> d1 >> d2;
// }
// 
// void readAngle( const std::string& line, char * str1,
// 		       char * str2, char * str3,
// 		       mdpack::dreal& d1,mdpack::dreal& d2){
//   std::istringstream iss(line,std::istringstream::in);
//   iss >> str1 >> str2 >> str3 >> d1 >> d2;
// }
//   
// 
// 
// 
// 
// void readAngle( const std::string& line, char * str1,
// 		       char * str2, char * str3,
// 		       mdpack::dreal& d1,mdpack::dreal& d2,
// 		       mdpack::dreal& d3,mdpack::dreal& d4){
//   
//   std::istringstream iss(line,std::istringstream::in);
//   iss >> str1 >> str2 >> str3 >> d1 >> d2 >> d3 >> d4;
// }
// 
// void readDihedral(const std::string& line, char * str1, char * str2, 
// 		  char * str3, char * str4,  mdpack::dreal& k_chi,
// 		  mdpack::dreal& multip, mdpack::dreal& delta){
//   
//   std::istringstream iss(line,std::istringstream::in);
//   iss >> str1 >> str2 >> str3 >> str4>> k_chi >> multip >> delta;
// }
// 
// 
// void readImproper(const std::string& line, char * str1, char * str2, char * str3, 
// 		  char * str4, mdpack::dreal& k_psi,mdpack::dreal& unused, mdpack::dreal& psi_0){
//   
//   std::istringstream iss(line,std::istringstream::in);
//   iss >> str1 >> str2 >> str3 >> str4>> k_psi >> unused >> psi_0;
// }
// 
// void readNonBonded(const std::string& line, char * str1, mdpack::dreal& unused, 
// 		   mdpack::dreal& epsilon, mdpack::dreal& rmin){
//   
//   std::istringstream iss(line,std::istringstream::in);
//   iss >> str1 >> unused >> epsilon >> rmin ;
// }
  