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

AtomStyle(difdif,AtomVecDIFDIF)

#else

#ifndef LMP_ATOM_VEC_DIFDIF_H
#define LMP_ATOM_VEC_DIFDIF_H

#include "atom_vec.h"

namespace LAMMPS_NS {

class AtomVecDIFDIF : public AtomVec {
 public:
  AtomVecDIFDIF(class LAMMPS *);
  void init();

  void grow_pointers();
  int property_atom(char *);
  void pack_property_atom(int, double *, int, int);

 private:
  double *gammai;
};

}

#endif
#endif
