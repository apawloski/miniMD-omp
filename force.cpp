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
#include "math.h"
#include "force.h"
#include "string.h"

#if defined(_OPENMP)
#include <omp.h>
#endif

Force::Force() {}
Force::~Force() {}

void Force::setup()
{
  cutforcesq = cutforce*cutforce;
}

void Force::compute(Atom &atom, Neighbor &neighbor, int me)
{

#if defined(_OPENMP)
#pragma omp parallel default(none) shared(atom,neighbor,me)
#endif
{
  int i,j,k,numneigh;
  double xtmp,ytmp,ztmp,delx,dely,delz,rsq;
  double fxtmp,fytmp,fztmp,sr2,sr6,force;
  int *neighs;

  const int nlocal = atom.nlocal;
  const int nall = atom.nlocal + atom.nghost;

#if defined(_OPENMP)
  const int tid = omp_get_thread_num();
  const int nthreads = atom.nthreads;

  // each thread works on a fixed chunk of atoms.
  const int idelta = 1 + nlocal/nthreads;
  const int ifrom = tid*idelta;
  const int imax  = ifrom + idelta;
  const int ito   = (imax > nlocal) ? nlocal : imax;

  double **f = atom.f + nall*tid;
#else
  const int tid = 0;
  const int ifrom = 0;
  const int ito = nlocal;
  double **f = atom.f;
#endif

  double **x = atom.x;

  // clear force on own and ghost atoms
  memset(&(f[0][0]),0,3*nall*sizeof(double));

  // loop over all neighbors of my atoms
  // store force on both atoms i and j
  
  for (i = ifrom; i < ito; i++) {
    neighs = neighbor.firstneigh[i];
    numneigh = neighbor.numneigh[i];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    fxtmp = fytmp = fztmp = 0.0;
    for (k = 0; k < numneigh; k++) {
      j = neighs[k];
      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;
      if (rsq < cutforcesq) {
	sr2 = 1.0/rsq;
	sr6 = sr2*sr2*sr2;
	force = sr6*(sr6-0.5)*sr2;
	fxtmp += delx*force;
	fytmp += dely*force;
	fztmp += delz*force;
	f[j][0] -= delx*force;
	f[j][1] -= dely*force;
	f[j][2] -= delz*force;
      }
    }
    f[i][0] += fxtmp;
    f[i][1] += fytmp;
    f[i][2] += fztmp;
  }
#if defined(_OPENMP)
  // NO-OP for one thread
#pragma omp barrier
  if (nthreads > 1) {
    // reduction of per-thread forces into first part of array
    const int nvals = 3*nall;
    const int idelta = nvals/nthreads + 1;
    const int ifrom = tid*idelta;
    const int imax  = ifrom + idelta;
    const int ito   = (imax > nvals) ? nvals : imax;
    double *fall = &(atom.f[0][0]);

    for (int m = ifrom; m < ito; ++m)
      for (int n = 1; n < nthreads; ++n)
        fall[m] += fall[n*nvals + m];
  }
#else
  // NO-OP in non-threaded compile
#endif
}
}
