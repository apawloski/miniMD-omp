/* ----------------------------------------------------------------------
   miniMD is a simple, parallel molecular dynamics (MD) code.   miniMD is
   an MD microapplication in the Mantevo project at Sandia National 
   Laboratories ( http://software.sandia.gov/mantevo/ ). The primary 
   authors of miniMD are Steve Plimpton and Paul Crozier 
   (pscrozi@sandia.gov).

   Copyright (2008) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This library is free software; you 
   can redistribute it and/or modify it under the terms of the GNU Lesser 
   General Public License as published by the Free Software Foundation; 
   either version 3 of the License, or (at your option) any later 
   version.
  
   This library is distributed in the hope that it will be useful, but
   WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU 
   Lesser General Public License for more details.
    
   You should have received a copy of the GNU Lesser General Public 
   License along with this software; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
   USA.  See also: http://www.gnu.org/licenses/lgpl.txt .

   For questions, contact Paul S. Crozier (pscrozi@sandia.gov). 

   Please read the accompanying README and LICENSE files.
---------------------------------------------------------------------- */

#include "stdio.h"
#include "integrate.h"

Integrate::Integrate() {}
Integrate::~Integrate() {}

void Integrate::setup()
{
  dtforce = 48.0*dt;
}

void Integrate::run(Atom &atom, Force &force, Neighbor &neighbor,
		    Comm &comm, Thermo &thermo, Timer &timer)
{
  int i,n;
  double *x,*y,*z,*v,*f,*vold;

  for (n = 0; n < ntimes; n++) {

    x = atom.x;
    y = atom.y;
    z = atom.z;
    v = &(atom.v[0][0]);

    const int n3local = 3*atom.nlocal;
    const double dt1 = dt;

#if defined(_OPENMP)
#pragma omp parallel for private(i) schedule(static) default(none) shared(x,y,z,v)
#endif
    for (i = 0; i < n3local; i++) {
      //We loop over each atom's 3 dimensions
      if ( i%3 == 0 ) {	
	x[i/3] += dt1*v[i];
      }
      else if ( i%3 == 1 ) {
	y[i/3] += dt1*v[i];
      }
      else {
	z[i/3] += dt1*v[i];      
      }
    }
    
    timer.stamp();

    if ((n+1) % neighbor.every) {
      comm.communicate(atom);
      timer.stamp(TIME_COMM);
    } else {
      comm.exchange(atom);
      comm.borders(atom);
      timer.stamp(TIME_COMM);
      neighbor.build(atom);
      timer.stamp(TIME_NEIGH);
    }      

    force.compute(atom,neighbor,comm.me);
    timer.stamp(TIME_FORCE);

    comm.reverse_communicate(atom);
    timer.stamp(TIME_COMM);

    vold = &(atom.vold[0][0]);
    f = &(atom.f[0][0]);

    const int n3local2 = 3*atom.nlocal;
    const double dtforce1 = dtforce;

#if defined(_OPENMP)
#pragma omp parallel for private(i) schedule(static) default(none) shared(vold,v,f)
#endif
    for (i = 0; i < n3local2; i++) {
      vold[i] = v[i];
      v[i] += dtforce1*f[i];
    }

    if (thermo.nstat) thermo.compute(n+1,atom,neighbor,force);
  }
}
