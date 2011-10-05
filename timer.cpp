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

#include "stdlib.h"
#include "mpi.h"
#include "timer.h"

Timer::Timer()
{
  array = (double *) malloc(TIME_N*sizeof(double));
  reset();
}

Timer::~Timer()
{
  if (array) free(array);
}

void Timer::reset()
{
  for (int i = 0; i < TIME_N; i++) array[i] = 0.0;
}

void Timer::stamp()
{
  previous_time = MPI_Wtime();
}

void Timer::stamp(int which)
{
  double current_time = MPI_Wtime();
  array[which] += current_time - previous_time;
  previous_time = current_time;
}

void Timer::barrier_start(int which)
{
  MPI_Barrier(MPI_COMM_WORLD);
  array[which] = MPI_Wtime();
}

void Timer::barrier_stop(int which)
{
  MPI_Barrier(MPI_COMM_WORLD);
  double current_time = MPI_Wtime();
  array[which] = current_time - array[which];
}
