/* ----------------------------------------------------------------------
 LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
 https://lammps.sandia.gov/, Sandia National Laboratories
 Steve Plimpton, sjplimp@sandia.gov

 Copyright (2003) Sandia Corporation.  Under the terms of Contract
 DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
 certain rights in this software.  This software is distributed under
 the GNU General Public License.

 See the README file in the top-level LAMMPS directory.
 ------------------------------------------------------------------------- */

#include "fix_pi.h"
#include "atom.h"
#include "force.h"
#include "update.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixPI::FixPI(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg) {

  if ((atom->esph_flag != 1) || (atom->rho_flag != 1))
    error->all(FLERR,
        "Fix sph command requires atom_style with both energy and density");

  if (narg != 3)
    error->all(FLERR,"Illegal number of arguments for fix sph command");

  time_integrate = 1;
}

/* ---------------------------------------------------------------------- */

int FixPI::setmask() {
  int mask = 0;
  mask |= INITIAL_INTEGRATE;
  mask |= FINAL_INTEGRATE;
  mask |= PRE_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixPI::init() {
  dtv = update->dt;
  dtf = 0.5 * update->dt * force->ftm2v;
}

void FixPI::setup_pre_force(int /*vflag*/)
{
  // set vest equal to v
  double **v = atom->v;
  double **vest = atom->vest;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup)
    nlocal = atom->nfirst;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      vest[i][0] = v[i][0];
      vest[i][1] = v[i][1];
      vest[i][2] = v[i][2];
    }
  }
}

void FixPI::pre_force(int /*vflag*/)
{
  // set vest equal to v
 // double **stfluid = atom->stfluid;

//  double **gradv = atom->gradv;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup)
    nlocal = atom->nfirst;

  //for (int i = 0; i < nlocal; i++) {
  //  if (mask[i] & groupbit) {
     // gradv[i][0] = 0.0;
     // gradv[i][1] = 0.0;
     // gradv[i][2] = 0.0;
     // gradv[i][3] = 0.0;
     // gradv[i][4] = 0.0;
     // gradv[i][5] = 0.0;
     // gradv[i][6] = 0.0;
      //gradv[i][7] = 0.0;
     // gradv[i][8] = 0.0;

      
    //}
  //}
}

/* ----------------------------------------------------------------------
 allow for both per-type and per-atom mass
 ------------------------------------------------------------------------- */

void FixPI::initial_integrate(int /*vflag*/) {
  // update v and x and rho and e of atoms in group
  double **gradv = atom->gradv;
 // double **stfluid = atom->stfluid;

  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  double **vest = atom->vest;
  double *rho = atom->rho;
  double *drho = atom->drho;
  double *esph = atom->esph;
  double *desph = atom->desph;
  double *mass = atom->mass;
  double *rmass = atom->rmass;
  int rmass_flag = atom->rmass_flag;

  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  int i;
  double dtfm;

  if (igroup == atom->firstgroup)
    nlocal = atom->nfirst;

  for (i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      if (rmass_flag) {
        dtfm = dtf / rmass[i];
      } else {
        dtfm = dtf / mass[type[i]];
      }

          //Only 2d to test 
      // double trace;
       //trace = gradv[i][0]+gradv[i][1];//+gradv[i][0]    
      // stfluid[i][0] =  0;//(2*gradv[i][0]+trace/2);
      // stfluid[i][1] =  0;//(2*gradv[i][1]+trace/2);
      // stfluid[i][2] =  0.0;//(gradv[i][2]-trace/2);
      // stfluid[i][3] =  (gradv[i][3]+gradv[i][6]);
     // stfluid[i][4] =  0.0;//0.5*(gradv[i][4]+gradv[i][7]);
     //  stfluid[i][5] =  0.0;//0.5*(gradv[i][5]+gradv[i][8]);

      gradv[i][0] = 0.0;
      gradv[i][1] = 0.0;
      gradv[i][2] = 0.0;
      gradv[i][3] = 0.0;
      gradv[i][4] = 0.0;
      gradv[i][5] = 0.0;
      gradv[i][6] = 0.0;
      gradv[i][7] = 0.0;
      gradv[i][8] = 0.0;
      
      esph[i] += dtf * desph[i]; // half-step update of particle internal energy
      rho[i] += dtf * drho[i]; // ... and density

      // extrapolate velocity for use with velocity-dependent potentials, e.g. SPH
      vest[i][0] = v[i][0] + 2.0 * dtfm * f[i][0];
      vest[i][1] = v[i][1] + 2.0 * dtfm * f[i][1];
      vest[i][2] = v[i][2] + 2.0 * dtfm * f[i][2];

      v[i][0] += dtfm * f[i][0];
      v[i][1] += dtfm * f[i][1];
      v[i][2] += dtfm * f[i][2];

      x[i][0] += dtv * v[i][0];
      x[i][1] += dtv * v[i][1];
      x[i][2] += dtv * v[i][2];
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixPI::final_integrate() {
  double **stfluid = atom->stfluid;
  double **gradv = atom->gradv;
  // update v, rho, and e of atoms in group

  double **v = atom->v;
  double **f = atom->f;
  double *esph = atom->esph;
  double *desph = atom->desph;
  double *rho = atom->rho;
  double *drho = atom->drho;
  int *type = atom->type;
  int *mask = atom->mask;
  double *mass = atom->mass;
  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup)
    nlocal = atom->nfirst;
  double dtfm;
  double *rmass = atom->rmass;
  int rmass_flag = atom->rmass_flag;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {

      if (rmass_flag) {
        dtfm = dtf / rmass[i];
      } else {
        dtfm = dtf / mass[type[i]];
      }
      v[i][0] += dtfm * f[i][0];
      v[i][1] += dtfm * f[i][1];
      v[i][2] += dtfm * f[i][2];

      esph[i] += dtf * desph[i];
      rho[i] += dtf * drho[i]; // ****** MODIFIED TO NOT UPDTADE RHO
       
       //Only 2d to test 
       //double trace;
       //trace = gradv[i][0]+gradv[i][1];//+gradv[i][0]    
       ///stfluid[i][0] =  (2*gradv[i][0]);//-trace/2);    ///this iscretization including zeta recovers properly the profile
       //stfluid[i][1] =  (2*gradv[i][1]);//-trace/2);
       //stfluid[i][2] =  0.0;//(gradv[i][2]-trace/2);
       //stfluid[i][3] =  (gradv[i][3]+gradv[i][6]);
       //stfluid[i][4] =  0.0;//0.5*(gradv[i][4]+gradv[i][7]);
      // stfluid[i][5] =  0.0;//0.5*(gradv[i][5]+gradv[i][8]);


    }
  }
}

/* ---------------------------------------------------------------------- */

void FixPI::reset_dt() {
  dtv = update->dt;
  dtf = 0.5 * update->dt * force->ftm2v;
}
