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


/**
 * Users can defined their functions or methods here. The methods can 
 * be defined as inline functions here or they can be defined in methods.cpp
 */

#ifndef METHODS_H
#define METHODS_H

#include "md.h"

/**==============================================================
 * To illustrate how to defined new method, we write here
 * the verlet method.
 */
inline void verlet(mdpack::md& md1){
  // first half of verlet
  for(int i=0; i<md1.nAtoms ; i++){
    // im is the inverse of mass, v is the velocity and r is the position]
    // h is the step-size, f is the force
    
    // updating velocities at half step
    md1.atoms[i].v += 0.5*md1.h*md1.f[i]*md1.atoms[i].im;
    // updating positions at full step
    md1.atoms[i].r += md1.h* md1.atoms[i].v; 
  }
  // kinetic energy and total energy, if you write integrator 
  // you must update kinetic and total energy
  md1.Ekinetic =0.0;
  md1.Etotal =0.0;
  // computing forces
  md1.computeForces();
  //second half of verlet
  for (int i=0; i< md1.nAtoms ; i++){
    md1.atoms[i].v += 0.5*md1.h*md1.f[i]*md1.atoms[i].im;
    // updating kinetic energy
    md1.Ekinetic +=0.5*(md1.atoms[i].m)*
		   (md1.atoms[i].v.x*md1.atoms[i].v.x +
		    md1.atoms[i].v.y*md1.atoms[i].v.y +
		    md1.atoms[i].v.z*md1.atoms[i].v.z);
  }
  md1.Etotal = md1.Ekinetic + md1.Epotential;
  // updating instantaneous temperature and avrage temperature
  // these must be updated by any integrator
  
  // updating instantaneous temperature
  md1.tempInstant = 2.0*INV_BOLTZMANN*(1.0/(md1.nDegsFreedom-4))*md1.Ekinetic;
  // updating avrage temperature 
  md1.tempAvg = (md1.stepn*md1.tempAvg +md1.tempInstant)/(md1.stepn +1.0);
  // updating simulation time
  md1.simTime +=md1.timeStep;
  // updating the step number
  md1.stepn +=1;
}
/**================================================================*/
  








#endif
