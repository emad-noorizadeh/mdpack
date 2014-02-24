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

#ifndef CONFIGOBJECT_H
#define CONFIGOBJECT_H

#include <cstring>
//#include <string>
#include "definitions.h"

namespace mdpack {
  
  /**
   * class configObject: provides easy way to add entry into config file 
   * and read them from config file
   * 
   */
  
  template <class Type1>
  class configObject {
    private:
      Type1 * field;
      std::string keyword;
      Type1 defaultValue;
    public:
      
      configObject(){ field = NULL; }
      
      configObject(Type1& T1, const char * strKeyword, const Type1& value){
	field = &T1;
	keyword= strKeyword;
	defaultValue = value;
      }
      
      bool check(const char * str){
	if( keyword == str)
	  return true;
	else
	  return false;
      }
      
      void setValue(const Type1& value){ *field= value;}
	
      
			   
  }; // end of configObject
    
  
}; //end of namespace mdpack



#endif  //end of CONFIGOBJECT_H
