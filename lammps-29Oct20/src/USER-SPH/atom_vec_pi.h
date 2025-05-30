/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef ATOM_CLASS

AtomStyle(pi,AtomVecPi)

#else

#ifndef LMP_ATOM_VEC_PI_H
#define LMP_ATOM_VEC_PI_H

#include "atom_vec.h"

namespace LAMMPS_NS {

class AtomVecPi : public AtomVec {
 public:
  AtomVecPi(class LAMMPS *);

  void grow_pointers();
  void force_clear(int, size_t);
  void create_atom_post(int);
  void data_atom_post(int);
  int property_atom(char *);
  void pack_property_atom(int, double *, int, int);

 private:
  double *rho,*drho,*esph,*desph,*cv;
  double **vest; // estimated velocity during force computation
  //Modified by NMC
  double **stmic; // stress tensor microstructure per particle
  double **stfluid; // fluid stress tensor 
  double **gradv; // velocity gradient tensor per particle

};

}

#endif
#endif
