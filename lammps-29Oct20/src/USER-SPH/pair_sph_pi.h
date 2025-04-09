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

#ifdef PAIR_CLASS

PairStyle(sph/pi,PairSPHPi)

#else

#ifndef LMP_PAIR_PI_H
#define LMP_PAIR_PI_H

#include "pair.h"

namespace LAMMPS_NS {

class PairSPHPi : public Pair {
 public:
  PairSPHPi(class LAMMPS *);
  virtual ~PairSPHPi();
  virtual void compute(int, int);
  void settings(int, char **);
  void coeff(int, char **);
  virtual double init_one(int, int);

 protected:
  double *rho0, *soundspeed, *B;
  double **cut;//,**viscosity;
  int first;
  
  double epsilonMM;  //Modified by NMC epsilon micro-macro. weighting parameter for micro-macro coupling
  double micro; // flag to use microscale info.
  double viscosity, pb, bulkViscosity; //pb is the background pressure, and visco and bulk visco correspond to

  void allocate();
};

}

#endif
#endif
