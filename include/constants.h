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

#ifndef CONSTANTS_H
#define CONSTANTS_H 

//const double TIMEFACTOR=48.88821;
const double TIMEFACTOR =20.459716518075222;
const double ITIMEFACTOR = 0.048876532532430;

const double COULOMB=332.0636;
const double ELECTRON = 1.6021892;
const double BOLTZMANN = 0.001987191;  //1.987191 * 10^-3 kcAL/(MOL*K)
const double INV_BOLTZMANN = 5.032228910054444e+02;

const double MAX =99999999999999999999999999.0;

const double HALFCUTOFF = 5.0;


const double pi=3.141592653589793238462643383279502884197;

// conversion factor from radian to degrees 180/pi
const double radTodeg =  57.295779513082323;
// conversion factor from degrees to radians pi/180
const double degTorad = 0.017453292519943;

// minimum atoms that required to gain efficiency
//by using cell lists
const int MINATOMS = 100;

const char ESPACE [] = "       ";





#endif // CONSTANTS_H
