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

#include "atom_vec_difdif.h"
#include <cstring>
#include "atom.h"
#include "update.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

AtomVecDIFDIF::AtomVecDIFDIF(LAMMPS *lmp) : AtomVec(lmp)
{
  molecular = Atom::ATOMIC;
  mass_type = PER_TYPE;
  forceclearflag = 1;

  atom->gammai_flag = 1;
  
  // strings with peratom variables to include in each AtomVec method
  // strings cannot contain fields in corresponding AtomVec default strings
  // order of fields in a string does not matter
  // except: fields_data_atom & fields_data_vel must match data file

  fields_grow = (char *) "gammai";
  fields_copy = (char *) "gammai";
  fields_comm = (char *) "gammai";
  fields_comm_vel = (char *) "gammai";
  fields_reverse = (char *) "";
  fields_border = (char *) "";
  fields_border_vel = (char *) "";
  fields_exchange = (char *) "gammai";
  fields_restart = (char * ) "gammai";
  fields_create = (char *) "gammai";
  fields_data_atom = (char *) "id type x gammai";
  fields_data_vel = (char *) "id v";

  setup_fields();
}

/* ---------------------------------------------------------------------- */

void AtomVecDIFDIF::init()
{
  AtomVec::init();

  if (strcmp(update->unit_style,"lj") != 0)
    error->all(FLERR,"Atom style DIFDIF requires lj units");
}

/* ----------------------------------------------------------------------
   set local copies of all grow ptrs used by this class, except defaults
   needed in replicate when 2 atom classes exist and it calls pack_restart()
------------------------------------------------------------------------- */

void AtomVecDIFDIF::grow_pointers()
{
  gammai = atom->gammai;

}


/* ----------------------------------------------------------------------
   assign an index to named atom property and return index
   return -1 if name is unknown to this atom style
------------------------------------------------------------------------- */

int AtomVecDIFDIF::property_atom(char *name)
{
  if (strcmp(name,"gammai") == 0) return 0;
  return -1;
}

/* ----------------------------------------------------------------------
   pack per-atom data into buf for ComputePropertyAtom
   index maps to data specific to this atom style
------------------------------------------------------------------------- */

void AtomVecDIFDIF::pack_property_atom(int index, double *buf,
                                     int nvalues, int groupbit)
{
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  int n = 0;
  if (index == 0) {
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) buf[n] = gammai[i];
      else buf[n] = 0.0;
      n += nvalues;
    }
  } 
}
