/* ----------------------------------------------------------------------
 LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
 https://lammps.sandia.gov/, Sandia National Laboratories
 Steve Plimpton, sjplimp@sandia.gov

 Copyright (2003) Sandia Corporation.  Under the terms of Contract
 DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
 certain rights in this software.  This software is distributed under
 the GNU General Public License.

 See the README file in the top-level LAMMPS directory.

 Modifed by NMC 03/18/2021 from sph style atom. Including settings for stmic, stfluid and gradv

 ------------------------------------------------------------------------- */

#include "atom_vec_pi.h"

#include "atom.h"

#include <cstring>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

AtomVecPi::AtomVecPi(LAMMPS *lmp) : AtomVec(lmp)
{
  molecular = Atom::ATOMIC;
  mass_type = PER_TYPE;
  forceclearflag = 1;

  atom->esph_flag = 1;
  atom->rho_flag = 1;
  atom->cv_flag = 1;
  atom->vest_flag = 1;
  atom->stmic_flag =1; // Modifed by NMC
  atom->stfluid_flag =1; // Modifed by NMC
  atom->gradv_flag =1; // Modifed by NMC

  // strings with peratom variables to include in each AtomVec method
  // strings cannot contain fields in corresponding AtomVec default strings
  // order of fields in a string does not matter
  // except: fields_data_atom & fields_data_vel must match data file

  fields_grow = (char *) "rho drho esph desph cv vest stmic stfluid gradv";
  fields_copy = (char *) "rho drho esph desph cv vest stmic stfluid gradv";
  fields_comm = (char *) "rho esph vest stmic stfluid gradv";
  fields_comm_vel = (char *) "rho esph vest stmic stfluid gradv";
  fields_reverse = (char *) "drho desph gradv";  //gradv is needed to communicate for sph pi for macroscale to compute gradv
  fields_border = (char *) "rho esph cv vest stmic stfluid gradv";
  fields_border_vel = (char *) "rho esph cv vest stmic stfluid gradv";
  fields_exchange = (char *) "rho esph cv vest stmic stfluid gradv";
  fields_restart = (char * ) "rho esph cv vest stmic stfluid gradv";
  fields_create = (char *) "rho esph cv vest desph drho stmic stfluid gradv";
  fields_data_atom = (char *) "id type rho esph cv x";
  fields_data_vel = (char *) "id v";

  setup_fields();
}

/* ----------------------------------------------------------------------
   set local copies of all grow ptrs used by this class, except defaults
   needed in replicate when 2 atom classes exist and it calls pack_restart()
------------------------------------------------------------------------- */

void AtomVecPi::grow_pointers()
{
  rho = atom->rho;
  drho = atom->drho;
  esph = atom->esph;
  desph = atom->desph;
  cv = atom->cv;
  vest = atom->vest;
  //Modified by NMC
  stmic = atom->stmic;
  stfluid = atom->stfluid;
  gradv = atom->gradv;
}

/* ----------------------------------------------------------------------
   clear extra forces starting at atom N
   nbytes = # of bytes to clear for a per-atom vector
------------------------------------------------------------------------- */

void AtomVecPi::force_clear(int n, size_t nbytes)
{
  memset(&desph[n],0,nbytes);
  memset(&drho[n],0,nbytes);

  // Modified by NMC. Gradv can be reinitialized after force computation
  
  memset(&gradv[n][0],0,nbytes);
  memset(&gradv[n][1],0,nbytes);
  memset(&gradv[n][2],0,nbytes);
  memset(&gradv[n][3],0,nbytes);
  memset(&gradv[n][4],0,nbytes);
  memset(&gradv[n][5],0,nbytes);
  memset(&gradv[n][6],0,nbytes);
  memset(&gradv[n][7],0,nbytes);
  memset(&gradv[n][8],0,nbytes);

}

/* ----------------------------------------------------------------------
   initialize non-zero atom quantities
------------------------------------------------------------------------- */

void AtomVecPi::create_atom_post(int ilocal)
{
  cv[ilocal] = 1.0;
}

/* ----------------------------------------------------------------------
   modify what AtomVec::data_atom() just unpacked
   or initialize other atom quantities
------------------------------------------------------------------------- */

void AtomVecPi::data_atom_post(int ilocal)
{
  vest[ilocal][0] = 0.0;
  vest[ilocal][1] = 0.0;
  vest[ilocal][2] = 0.0;
  desph[ilocal] = 0.0;
  drho[ilocal] = 0.0;

  // Modifed by NMC
  stmic[ilocal][0] = 0.0;
  stmic[ilocal][1] = 0.0;
  stmic[ilocal][2] = 0.0;
  stmic[ilocal][3] = 0.0;
  stmic[ilocal][4] = 0.0;
  stmic[ilocal][5] = 0.0;

  stfluid[ilocal][0] = 0.0;
  stfluid[ilocal][1] = 0.0;
  stfluid[ilocal][2] = 0.0;
  stfluid[ilocal][3] = 0.0;
  stfluid[ilocal][4] = 0.0;
  stfluid[ilocal][5] = 0.0;

// Modifed by NMC
  
  gradv[ilocal][0] = 0.0;
  gradv[ilocal][1] = 0.0;
  gradv[ilocal][2] = 0.0;
  gradv[ilocal][3] = 0.0;
  gradv[ilocal][4] = 0.0;
  gradv[ilocal][5] = 0.0;
  gradv[ilocal][6] = 0.0;
  gradv[ilocal][7] = 0.0;
  gradv[ilocal][8] = 0.0;
  
}

/* ----------------------------------------------------------------------
   assign an index to named atom property and return index
   return -1 if name is unknown to this atom style
------------------------------------------------------------------------- */

int AtomVecPi::property_atom(char *name)
{
  if (strcmp(name,"rho") == 0) return 0;
  if (strcmp(name,"drho") == 0) return 1;
  if (strcmp(name,"esph") == 0) return 2;
  if (strcmp(name,"desph") == 0) return 3;
  if (strcmp(name,"cv") == 0) return 4;


  //////////////Return each component of the tensors separately
  if (strcmp(name,"stfxx") == 0) return 5; //Modified by NMC
  if (strcmp(name,"stfyy") == 0) return 6; //Modified by NMC
  if (strcmp(name,"stfzz") == 0) return 7; //Modified by NMC
  if (strcmp(name,"stfxy") == 0) return 8; //Modified by NMC
  if (strcmp(name,"stfxz") == 0) return 9; //Modified by NMC
  if (strcmp(name,"stfyz") == 0) return 10; //Modified by NMC

  if (strcmp(name,"stmxx") == 0) return 11; //Modified by NMC
  if (strcmp(name,"stmyy") == 0) return 12; //Modified by NMC
  if (strcmp(name,"stmzz") == 0) return 13; //Modified by NMC
  if (strcmp(name,"stmxy") == 0) return 14; //Modified by NMC
  if (strcmp(name,"stmxz") == 0) return 15; //Modified by NMC
  if (strcmp(name,"stmyz") == 0) return 16; //Modified by NMC

  if (strcmp(name,"gdvxx") == 0) return 17; //Modified by NMC
  if (strcmp(name,"gdvyy") == 0) return 18; //Modified by NMC
  if (strcmp(name,"gdvzz") == 0) return 19; //Modified by NMC
  if (strcmp(name,"gdvxy") == 0) return 20; //Modified by NMC
  if (strcmp(name,"gdvxz") == 0) return 21; //Modified by NMC
  if (strcmp(name,"gdvyz") == 0) return 22; //Modified by NMC
  if (strcmp(name,"gdvyx") == 0) return 23; //Modified by NMC
  if (strcmp(name,"gdvzx") == 0) return 24; //Modified by NMC
  if (strcmp(name,"gdvzy") == 0) return 25; //Modified by NMC

  return -1;
}

/* ----------------------------------------------------------------------
   pack per-atom data into buf for ComputePropertyAtom
   index maps to data specific to this atom style
------------------------------------------------------------------------- */

void AtomVecPi::pack_property_atom(int index, double *buf,
                                     int nvalues, int groupbit)
{
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  int n = 0;

  if (index == 0) {
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) buf[n] = rho[i];
      else buf[n] = 0.0;
      n += nvalues;
    }
  } else if (index == 1) {
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) buf[n] = drho[i];
      else buf[n] = 0.0;
      n += nvalues;
    }
  } else if (index == 2) {
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) buf[n] = esph[i];
      else buf[n] = 0.0;
      n += nvalues;
    }
  } else if (index == 3) {
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) buf[n] = desph[i];
      else buf[n] = 0.0;
      n += nvalues;
    }
  } else if (index == 4) {
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) buf[n] = cv[i];
      else buf[n] = 0.0;
      n += nvalues;
    }
  } else if (index == 5) {  //Modified from here to include the prop atome for the new tensors
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) buf[n] = stfluid[i][0];
      else buf[n] = 0.0;
      n += nvalues;
    }
  } else if (index == 6) {
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) buf[n] = stfluid[i][1];
      else buf[n] = 0.0;
      n += nvalues;
    }
  } else if (index == 7) {
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) buf[n] = stfluid[i][2];
      else buf[n] = 0.0;
      n += nvalues;
    }
  } else if (index == 8) {
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) buf[n] = stfluid[i][3];
      else buf[n] = 0.0;
      n += nvalues;
    }
  } else if (index == 9) {
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) buf[n] = stfluid[i][4];
      else buf[n] = 0.0;
      n += nvalues;
    }
  } else if (index == 10) {
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) buf[n] = stfluid[i][5];
      else buf[n] = 0.0;
      n += nvalues;
    }
  } else if (index == 11) {
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) buf[n] = stmic[i][0];
      else buf[n] = 0.0;
      n += nvalues;
    }
  } else if (index == 12) {
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) buf[n] = stmic[i][1];
      else buf[n] = 0.0;
      n += nvalues;
    }
  } else if (index == 13) {
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) buf[n] = stmic[i][2];
      else buf[n] = 0.0;
      n += nvalues;
    }
  } else if (index == 14) {
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) buf[n] = stmic[i][3];
      else buf[n] = 0.0;
      n += nvalues;
    }
  } else if (index == 15) {
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) buf[n] = stmic[i][4];
      else buf[n] = 0.0;
      n += nvalues;
    }
  } else if (index == 16) {
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) buf[n] = stmic[i][5];
      else buf[n] = 0.0;
      n += nvalues;
    }
  } else if (index == 17) {
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) buf[n] = gradv[i][0];
      else buf[n] = 0.0;
      n += nvalues;
    }
  } else if (index == 18) {
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) buf[n] = gradv[i][1];
      else buf[n] = 0.0;
      n += nvalues;
    }
  } else if (index == 19) {
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) buf[n] = gradv[i][2];
      else buf[n] = 0.0;
      n += nvalues;
    }
  } else if (index == 20) {
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) buf[n] = gradv[i][3];
      else buf[n] = 0.0;
      n += nvalues;
    }
  } else if (index == 21) {
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) buf[n] = gradv[i][4];
      else buf[n] = 0.0;
      n += nvalues;
    }
  } else if (index == 22) {
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) buf[n] = gradv[i][5];
      else buf[n] = 0.0;
      n += nvalues;
    }
  } else if (index == 23) {
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) buf[n] = gradv[i][6];
      else buf[n] = 0.0;
      n += nvalues;
    }
  } else if (index == 24) {
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) buf[n] = gradv[i][7];
      else buf[n] = 0.0;
      n += nvalues;
    }
  } else if (index == 25) {
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) buf[n] = gradv[i][8];
      else buf[n] = 0.0;
      n += nvalues;
    }
  }
}
