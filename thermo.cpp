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
#include "stdlib.h"
#include "mpi.h"

#include "thermo.h"

Thermo::Thermo() {}
Thermo::~Thermo() {}

void Thermo::setup(double rho_in, int ntimes_in)
{
  rho = rho_in;
  ntimes = ntimes_in;

  int maxstat;
  if (nstat == 0) maxstat = 2;
  else maxstat = ntimes/nstat + 1;

  steparr = (int *) malloc(maxstat*sizeof(int));
  tmparr = (double *) malloc(maxstat*sizeof(double));
  engarr = (double *) malloc(maxstat*sizeof(double));
  prsarr = (double *) malloc(maxstat*sizeof(double));
}

void Thermo::compute(int iflag, Atom &atom, Neighbor &neighbor, Force &force)
{
  double t,eng,p;

  if (iflag > 0 && iflag % nstat) return;
  if (iflag == -1 && nstat > 0 && ntimes % nstat == 0) return;

  t = temperature(atom);
  eng = energy(atom,neighbor,force);

  force.compute(atom,neighbor,-1);
  p = pressure(t,atom);

  int istep = iflag;
  if (iflag == -1) istep = ntimes;

  if (iflag == 0) mstat = 0;

  steparr[mstat] = istep;
  tmparr[mstat] = t;
  engarr[mstat] = eng;
  prsarr[mstat] = p;

  mstat++;
}

/* reduced potential energy */

double Thermo::energy(Atom &atom, Neighbor &neighbor, Force &force)
{
  int i;
  double eng = 0.0;
#if defined(_OPENMP)
#pragma omp parallel for private(i) default(none) reduction(+:eng) shared(atom,neighbor,force)
#endif
  for (i = 0; i < atom.nlocal; i++) {
    int j,k;
    double delx,dely,delz,rsq,sr2,sr6,phi;

    const int * const neighs = neighbor.firstneigh[i];
    const int numneigh = neighbor.numneigh[i];
    for (k = 0; k < numneigh; k++) {
      j = neighs[k];
      delx = atom.x[i] - atom.x[j];
      dely = atom.y[i] - atom.y[j];
      delz = atom.z[i] - atom.z[j];
      rsq = delx*delx + dely*dely + delz*delz;
      if (rsq < force.cutforcesq) {
	sr2 = 1.0/rsq;
	sr6 = sr2*sr2*sr2;
	phi = sr6*(sr6-1.0);
	eng += phi;
      }
    }
  }

  double engtmp = 4.0*eng;
  MPI_Allreduce(&engtmp,&eng,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  return eng/atom.natoms;
}

/*  reduced temperature */

double Thermo::temperature(Atom &atom)
{
  const double * const v = &(atom.v[0][0]);
  const double * const vold = &(atom.vold[0][0]);
  const int n3local = 3*atom.nlocal;

  int i;
  double t = 0.0;
#if defined(_OPENMP)
#pragma omp parallel for private(i) default(none) reduction(+:t)
#endif
  for (i = 0; i < n3local; i++) {
    const double vx = 0.5 * (v[i] + vold[i]);
    t += vx*vx;
  }

  double t1;
  MPI_Allreduce(&t,&t1,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  return t1/(3.0*atom.natoms);
}
      
      
/* reduced pressure from virial
   virial = Fi dot Ri summed over own and ghost atoms, since PBC info is
   stored correctly in force array before reverse_communicate is performed */

double Thermo::pressure(double t, Atom &atom)
{
  const double * const x = atom.x;
  const double * const f = &(atom.f[0][0]);
  const int n3all = 3*(atom.nlocal+atom.nghost);

  int i;
  double virial = 0.0;
#if defined(_OPENMP)
#pragma omp parallel for private(i) default(none) reduction(+:virial)
#endif
  for (i = 0; i < n3all; i++)
    virial += f[i]*x[i];

  double virtmp = 48.0*virial;
  MPI_Allreduce(&virtmp,&virial,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  return (t*rho + rho/3.0/atom.natoms * virial);
}
