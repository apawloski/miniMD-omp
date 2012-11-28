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
#include "string.h"
#include "stdlib.h"
#include "mpi.h"
#include "atom.h"

#if defined(_OPENMP)
#include <omp.h>
#endif

#define DELTA 20000

Atom::Atom()
{
  natoms = 0;
  nlocal = 0;
  nghost = 0;
  nmax = 0;
  nthreads = 1;
#if defined(_OPENMP)
#pragma omp parallel default(none)
  {
#pragma omp master 
    { nthreads = omp_get_num_threads(); }
  }
#endif

  v = f = vold = NULL;
  x = y = z = NULL;

  comm_size = 3;
  reverse_size = 3;
  border_size = 3;
}

Atom::~Atom()
{
  if (nmax) {
    destroy_double_array(x);
    destroy_double_array(y);
    destroy_double_array(z);
    destroy_2d_double_array(v);
    destroy_2d_double_array(f);
    destroy_2d_double_array(vold);
  }
}

void Atom::growarray()
{
  int nold = nmax;
  nmax += DELTA;
  x = (double *) realloc_double_array(x,nmax,3*nold);
  y = (double *) realloc_double_array(y,nmax,3*nold);
  z = (double *) realloc_double_array(z,nmax,3*nold);
  v = (double **) realloc_2d_double_array(v,nmax,3,3*nold);
  f = (double **) realloc_2d_double_array(f,nmax*nthreads,3,3*nold);
  vold = (double **) realloc_2d_double_array(vold,nmax,3,3*nold);
  if (x == NULL || y == NULL || z == NULL || v == NULL || f == NULL || vold == NULL) {
    printf("ERROR: No memory for atoms\n");
  }
}

void Atom::addatom(double x_in, double y_in, double z_in, 
		   double vx_in, double vy_in, double vz_in)
{
  if (nlocal == nmax) growarray();

  x[nlocal] = x_in;
  y[nlocal] = y_in;
  z[nlocal] = z_in;
  v[nlocal][0] = vx_in;
  v[nlocal][1] = vy_in;
  v[nlocal][2] = vz_in;

  nlocal++;
}

/* enforce PBC
   order of 2 tests is important to insure lo-bound <= coord < hi-bound
   even with round-off errors where (coord +/- epsilon) +/- period = bound */

void Atom::pbc()
{
  for (int i = 0; i < nlocal; i++) {
    if (x[i] < 0.0) x[i] += box.xprd;
    if (x[i] >= box.xprd) x[i] -= box.xprd;
    if (y[i] < 0.0) y[i] += box.yprd;
    if (y[i] >= box.yprd) y[i] -= box.yprd;
    if (z[i] < 0.0) z[i] += box.zprd;
    if (z[i] >= box.zprd) z[i] -= box.zprd;
  }
}

void Atom::copy(int i, int j)
{
  x[j] = x[i];
  y[j] = y[i];
  z[j] = z[i];
  v[j][0] = v[i][0];
  v[j][1] = v[i][1];
  v[j][2] = v[i][2];
}

void Atom::pack_comm(int n, int *list, double *buf, int *pbc_flags)
{
  int i,j,m;

  m = 0;
  if (pbc_flags[0] == 0) {
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = x[j];
      buf[m++] = y[j];
      buf[m++] = z[j];
    }
  } else {
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = x[j] + pbc_flags[1]*box.xprd;
      buf[m++] = y[j] + pbc_flags[2]*box.yprd;
      buf[m++] = z[j] + pbc_flags[3]*box.zprd;
    }
  }
}

void Atom::unpack_comm(int n, int first, double *buf)
{
  int i,j,m;
  m = 0;
  j = first;
  for (i = 0; i < n; i++, j++) {
    x[j] = buf[m++];
    y[j] = buf[m++];
    z[j] = buf[m++];
  }
}

void Atom::pack_reverse(int n, int first, double *buf)
{
  int i,j,m;
  m = 0;
  j = first;
  for (i = 0; i < n; i++, j++) {
    buf[m++] = f[j][0];
    buf[m++] = f[j][1];
    buf[m++] = f[j][2];
  }
}

void Atom::unpack_reverse(int n, int *list, double *buf)
{
  int i,j,m;
  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    f[j][0] += buf[m++];
    f[j][1] += buf[m++];
    f[j][2] += buf[m++];
  }
}

int Atom::pack_border(int i, double *buf, int *pbc_flags)
{
  int m = 0;
  if (pbc_flags[0] == 0) {
    buf[m++] = x[i];
    buf[m++] = y[i];
    buf[m++] = z[i];
  } else {
    buf[m++] = x[i] + pbc_flags[1]*box.xprd;
    buf[m++] = y[i] + pbc_flags[2]*box.yprd;
    buf[m++] = z[i] + pbc_flags[3]*box.zprd;
  }
  return m;
}

int Atom::unpack_border(int i, double *buf)
{
  if (i == nmax) growarray();

  int m = 0;
  x[i] = buf[m++];
  y[i] = buf[m++];
  z[i] = buf[m++];
  return m;
}

int Atom::pack_exchange(int i, double *buf)
{
  int m = 0;
  buf[m++] = x[i];
  buf[m++] = y[i];
  buf[m++] = z[i];
  buf[m++] = v[i][0];
  buf[m++] = v[i][1];
  buf[m++] = v[i][2];
  return m;
}

int Atom::unpack_exchange(int i, double *buf)
{
  if (i == nmax) growarray();

  int m = 0;
  x[i] = buf[m++];
  y[i] = buf[m++];
  z[i] = buf[m++];
  v[i][0] = buf[m++];
  v[i][1] = buf[m++];
  v[i][2] = buf[m++];
  return m;
}

int Atom::skip_exchange(double *buf)
{
  return 6;
}

/* realloc a 2-d double array */

double **Atom::realloc_2d_double_array(double **array, 
				       int n1, int n2, int nold)

{
  double **newarray;

  newarray = create_2d_double_array(n1,n2);
  if (nold) memcpy(newarray[0],array[0],nold*sizeof(double));
  destroy_2d_double_array(array);

  return newarray;
}

/* create a 2-d double array */

double **Atom::create_2d_double_array(int n1, int n2)

{
  double **array;
  double *data;
  int i,n;

  if (n1*n2 == 0) return NULL;

  array = (double **) malloc(n1*sizeof(double *));
  data = (double *) malloc(n1*n2*sizeof(double));

  n = 0;
  for (i = 0; i < n1; i++) {
    array[i] = &data[n];
    n += n2;
  }

  return array;
}

/* free memory of a 2-d double array */

void Atom::destroy_2d_double_array(double **array)

{
  if (array != NULL) {
    free(array[0]);
    free(array);
  }
}

/* The following 1-d double array methods are just to mimick 2-d array methods */

/* realloc a double array */

double *Atom::realloc_double_array(double *array, int n, int nold)

{
  double *newarray;

  newarray = create_double_array(n);
  if (nold) memcpy(newarray, array, nold*sizeof(double));
  destroy_double_array(array);

  return newarray;
}

/* create a double array */

double *Atom::create_double_array(int n)

{
  double *array;

  if (n == 0) return NULL;

  array = (double *) malloc(n*sizeof(double));

  return array;
}

/* free memory of a double array */

void Atom::destroy_double_array(double *array)

{
  if (array != NULL) {
    free(array);
  }
}
