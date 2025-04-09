/* ----------------------------------------------------------------------
 LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
 https://lammps.sandia.gov/, Sandia National Laboratories
 Steve Plimpton, sjplimp@sandia.gov

 Copyright (2003) Sandia Corporation.  Under the terms of Contract
 DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
 certain rights in this software.  This software is distributed under
 the GNU General Public License.

 See the README file in the top-level LAMMPS directory.

Modified version of SPH to include microscopic derived stress tensor component to the sph particles
Additionally it uses the viscosoti formulation of Espanol and Ravenga 2003
 ------------------------------------------------------------------------- */

#include "pair_sph_pi.h"
#include <cmath>
#include "atom.h"
#include "force.h"
#include "comm.h"
#include "neigh_list.h"
#include "memory.h"
#include "error.h"
#include "domain.h"


using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

PairSPHPi::PairSPHPi(LAMMPS *lmp) : Pair(lmp)
{
  restartinfo = 0;
  first = 1;
  single_enable = 0;
}

/* ---------------------------------------------------------------------- */

PairSPHPi::~PairSPHPi() {
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);

    memory->destroy(cut);
    memory->destroy(rho0);
    memory->destroy(soundspeed);
    memory->destroy(B);
   // memory->destroy(viscosity);
  }
}

/* ---------------------------------------------------------------------- */

void PairSPHPi::compute(int eflag, int vflag) {
  int i, j, ii, jj, inum, jnum, itype, jtype;
  double xtmp, ytmp, ztmp, delx, dely, delz, fpair;

  int *ilist, *jlist, *numneigh, **firstneigh;
  double vxtmp, vytmp, vztmp, imass, jmass, fi, fj, fvisc, h, ih, ihsq, velx, vely, velz;
  double rsq, tmp, wfd, delVdotDelR, deltaE;

  //Modified NMC
  double smxx, smyy, smzz, smxy, smxz, smyz;  //six components of the tensor microstrcuture
  double dsmxx, dsmyy, dsmzz, dsmxy, dsmxz, dsmyz;  //deltas for stress microstructucture
  double sflxx, sflyy, sflzz, sflxy, sflxz, sflyz;  //six components of the tensor fluid
  double dsflxx, dsflyy, dsflzz, dsflxy, dsflxz, dsflyz;  //deltas for stress fluid
  double fpi;

  double fmicx, fmicy, fmicz;  //force components from microscale
  double fflux, ffluy, ffluz;  //force components from microscale


  ev_init(eflag, vflag);

  //Modified NMC
  double **stmic = atom->stmic;
  double **stfluid = atom->stfluid;
  double **gradv = atom->gradv;
  double a,b,dim, fvisx, fvisy, fvisz;  //viscosity coefficients
  double prefKernel;
  int dimension = domain->dimension;

  double **v = atom->vest;
  double **x = atom->x;
  double **f = atom->f;
  double *rho = atom->rho;
  double *mass = atom->mass;
  double *desph = atom->desph;
  double *drho = atom->drho;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  int newton_pair = force->newton_pair;

  // check consistency of pair coefficients

  if (first) {
    for (i = 1; i <= atom->ntypes; i++) {
      for (j = 1; i <= atom->ntypes; i++) {
        if (cutsq[i][j] > 1.e-32) {
          if (!setflag[i][i] || !setflag[j][j]) {
            if (comm->me == 0) {
              printf(
                  "SPH-PI particle types %d and %d interact with cutoff=%g, but not all of their single particle properties are set.\n",
                  i, j, sqrt(cutsq[i][j]));
            }
          }
        }
      }
    }
    first = 0;
  }

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;
  
if(dimension==3){
  dim=3.0;
} else {
  dim=2.0;
}
//NMC added amplitude of thermal noise, accounting for different dimensionality
// to use with Espanol and Revenga model.

a = (2.0-1.0/dim)*viscosity - bulkViscosity;
b = (2.0+dim)/dim * viscosity +(2.0+dim)*bulkViscosity - a*(2.0*dim-4.0)/(2.0*dim); //b+1/3a  eqs 61 Espanol2003, //for 2D a contribution on b vanishes

prefKernel = -25.066903536973515383e0*(dim-2.0) - 19.098593171027440292e0*(3.0-dim); //first term for 3d second 2d

  // loop over neighbors of my atoms

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    vxtmp = v[i][0];
    vytmp = v[i][1];
    vztmp = v[i][2];
    itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    imass = mass[itype];

    //Modified by NMC
    smxx = stmic[i][0]/ (rho[i] * rho[i]);
    smyy = stmic[i][1]/ (rho[i] * rho[i]);
    smzz = stmic[i][2]/ (rho[i] * rho[i]);
    smxy = stmic[i][3]/ (rho[i] * rho[i]);
    smxz = stmic[i][4]/ (rho[i] * rho[i]);
    smyz = stmic[i][5]/ (rho[i] * rho[i]);

    sflxx = stfluid[i][0]/ (rho[i] * rho[i]);
    sflyy = stfluid[i][1]/ (rho[i] * rho[i]);
    sflzz = stfluid[i][2]/ (rho[i] * rho[i]);
    sflxy = stfluid[i][3]/ (rho[i] * rho[i]);
    sflxz = stfluid[i][4]/ (rho[i] * rho[i]);
    sflyz = stfluid[i][5]/ (rho[i] * rho[i]);

    // compute pressure of atom i with Tait EOS
    tmp = rho[i] / rho0[itype];
    fi = tmp * tmp * tmp;
    fi = (B[itype] * (fi * fi * tmp - 1.0)+pb) / (rho[i] * rho[i]);

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx * delx + dely * dely + delz * delz;
      jtype = type[j];
      jmass = mass[jtype];

      if (rsq < cutsq[itype][jtype]) {
        h = cut[itype][jtype];
        ih = 1.0 / h;
        ihsq = ih * ih;

        wfd = h - sqrt(rsq);
         //if (dimension == 3) {       Removing if loop since prefKernel and wfd account for that
          // Lucy Kernel, 3d
          // Note that wfd, the derivative of the weight function with respect to r,
          // is lacking a factor of r.
          // The missing factor of r is recovered by
          // (1) using delV . delX instead of delV . (delX/r) and
          // (2) using f[i][0] += delx * fpair instead of f[i][0] += (delx/r) * fpair
         // wfd = -25.066903536973515383e0 * wfd * wfd * ihsq * ihsq * ihsq * ih;
        //} else {
          // Lucy Kernel, 2d
         // wfd = -19.098593171027440292e0 * wfd * wfd * ihsq * ihsq * ihsq;
        //}

        wfd = prefKernel * wfd * wfd * ihsq * ihsq * ihsq * (1*(3.0-dim)+ih*(dim-2.0));

        // compute pressure  of atom j with Tait EOS
        tmp = rho[j] / rho0[jtype];
        fj = tmp * tmp * tmp;
//        fj = B[jtype] * (fj * fj * tmp - 1.0) / (rho[j] * rho[j]);
        fj = (B[jtype] * (fj * fj * tmp - 1.0)+pb) / (rho[j] * rho[j]); //term with background pressure

        velx=vxtmp - v[j][0];
        vely=vytmp - v[j][1];
        velz=vztmp - v[j][2];

///////////////////////////////////////////////////////////////
        //Modified by NMC
        ///Symmetric
        //Compute the difference on the stress between particles
        // dsmxx = smxx - stmic[j][0];
        //dsmyy = smyy - stmic[j][1];
        //dsmzz = smzz - stmic[j][2];
        //dsmxy = smxy - stmic[j][3];
        //dsmxz = smxz - stmic[j][4];
        //dsmyz = smyz - stmic[j][5];

       // dsflxx = sflxx - stfluid[j][0];
       // dsflyy = sflyy - stfluid[j][1];
       // dsflzz = sflzz - stfluid[j][2];
       // dsflxy = sflxy - stfluid[j][3];
       // dsflxz = sflxz - stfluid[j][4];
       // dsflyz = sflyz - stfluid[j][5]; 

        ///Asymmetric

        dsmxx = smxx + stmic[j][0]/(rho[j] * rho[j]);
        dsmyy = smyy + stmic[j][1]/(rho[j] * rho[j]);
        dsmzz = smzz + stmic[j][2]/(rho[j] * rho[j]);
        dsmxy = smxy + stmic[j][3]/(rho[j] * rho[j]);
        dsmxz = smxz + stmic[j][4]/(rho[j] * rho[j]);
        dsmyz = smyz + stmic[j][5]/(rho[j] * rho[j]);


        dsflxx = sflxx + stfluid[j][0]/(rho[j] * rho[j]); 
        dsflyy = sflyy + stfluid[j][1]/(rho[j] * rho[j]);
        dsflzz = sflzz + stfluid[j][2]/(rho[j] * rho[j]);
        dsflxy = sflxy + stfluid[j][3]/(rho[j] * rho[j]);
        dsflxz = sflxz + stfluid[j][4]/(rho[j] * rho[j]);
        dsflyz = sflyz + stfluid[j][5]/(rho[j] * rho[j]);        

        ///pi x X
        fmicx = dsmxx*delx + dsmxy*dely + dsmxz*delz; 
        fmicy = dsmxy*delx + dsmyy*dely + dsmyz*delz;
        fmicz = dsmxz*delx + dsmyz*dely + dsmzz*delz;

        fflux = dsflxx*delx + dsflxy*dely + dsflxz*delz; 
        ffluy = dsflxy*delx + dsflyy*dely + dsflyz*delz;
        ffluz = dsflxz*delx + dsflyz*dely + dsflzz*delz;
        
        // Temporal to test the validity of eta*gradient
        //fflux*=viscosity;
        //ffluy*=viscosity;
        //ffluz*=viscosity;
        

        //Approximating gradient of the velocity at that point using the extrapolated velocity vest
        gradv[i][0] += -(jmass  * wfd / (rho[j])) * velx*delx;
        gradv[i][1] += -(jmass  * wfd / (rho[j])) * vely*dely;
        gradv[i][2] += -(jmass  * wfd / (rho[j])) * velz*delz;
        gradv[i][3] += -(jmass  * wfd / (rho[j])) * velx*dely; 
        gradv[i][4] += -(jmass  * wfd / (rho[j])) * velx*delz;
        gradv[i][5] += -(jmass  * wfd / (rho[j])) * vely*delz;
        gradv[i][6] += -(jmass  * wfd / (rho[j])) * vely*delx; 
        gradv[i][7] += -(jmass  * wfd / (rho[j])) * velz*delx;
        gradv[i][8] += -(jmass  * wfd / (rho[j])) * velz*dely;

//////////////////////////////////////////////////////////

// total pair force & thermal energy increment
        fpair = -imass * jmass * (fi + fj) * wfd; 

        // dot product of velocity delta and distance vector
        delVdotDelR = delx * velx + dely * vely + delz * velz;

        // Espanol Viscosity (Espanol, 2003)
        fvisc = imass * jmass * wfd / (rho[i]*rho[j]);

       // printf("testvar= %f, %f \n", delx, dely);
        fvisx = (a*velx + (b+a*(2.0*dim-4.0)/(2.0*dim))*delx * delVdotDelR / rsq) * fvisc;
        fvisy = (a*vely + (b+a*(2.0*dim-4.0)/(2.0*dim))*dely * delVdotDelR / rsq) * fvisc;
        fvisz = (a*velz + (b+a*(2.0*dim-4.0)/(2.0*dim))*delz * delVdotDelR / rsq) * fvisc;

        // Morris Viscosity (Morris, 1996)
        //fvisc = 2 * viscosity / (rho[i] * rho[j]);
        //fvisc *= imass * jmass * wfd;
        deltaE = -0.5 *(fpair * delVdotDelR + fvisc * (velx*velx + vely*vely + velz*velz));  ////THIS NEED TO BE UPDATED

  /////////////// Need to double check the signs of the new terms 
        fpi = imass*jmass*wfd*micro;

        f[i][0] += delx * fpair + epsilonMM*(fvisx) + (fflux + fmicx)*fpi;
        f[i][1] += dely * fpair + epsilonMM*(fvisy) + (ffluy + fmicy)*fpi;
        f[i][2] += delz * fpair + epsilonMM*(fvisz) + (ffluz + fmicz)*fpi;

        // and change in density
        drho[i] += jmass * delVdotDelR * wfd;

        // change in thermal energy
        desph[i] += deltaE; 

        if (newton_pair || j < nlocal) {

          f[j][0] -= delx * fpair + epsilonMM*(fvisx) + (fflux + fmicx)*fpi;
          f[j][1] -= dely * fpair + epsilonMM*(fvisy) + (ffluy + fmicy)*fpi;
          f[j][2] -= delz * fpair + epsilonMM*(fvisz) + (ffluz + fmicz)*fpi;
          desph[j] += deltaE;
          drho[j] += imass * delVdotDelR * wfd;

          //Modified by NMC
          //Approximating gradient of the velocity at that point using the extrapolated velocity vest
          //needed to proper stimation of gradient with ghost atoms
          gradv[j][0] += -(imass * wfd / (rho[i])) * delx*velx;
          gradv[j][1] += -(imass * wfd / (rho[i])) * dely*vely;
          gradv[j][2] += -(imass * wfd / (rho[i])) * delz*velz;
          gradv[j][3] += -(imass * wfd / (rho[i])) * dely*velx;
          gradv[j][4] += -(imass * wfd / (rho[i])) * velx*delz;
          gradv[j][5] += -(imass * wfd / (rho[i])) * vely*delz;
          gradv[j][6] += -(imass * wfd / (rho[i])) * vely*delx; 
          gradv[j][7] += -(imass * wfd / (rho[i])) * velz*delx;
          gradv[j][8] += -(imass * wfd / (rho[i])) * velz*dely;
        }

        if (evflag)
          ev_tally(i, j, nlocal, newton_pair, 0.0, 0.0, fpair, delx, dely, delz);
      }
    }
  }

  if (vflag_fdotr) virial_fdotr_compute();
}

/* ----------------------------------------------------------------------
 allocate all arrays
 ------------------------------------------------------------------------- */

void PairSPHPi::allocate() {
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag, n + 1, n + 1, "pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  memory->create(cutsq, n + 1, n + 1, "pair:cutsq");

  memory->create(rho0, n + 1, "pair:rho0");
  memory->create(soundspeed, n + 1, "pair:soundspeed");
  memory->create(B, n + 1, "pair:B");
  memory->create(cut, n + 1, n + 1, "pair:cut");
  //memory->create(viscosity, n + 1, n + 1, "pair:viscosity");
}

/* ----------------------------------------------------------------------
 global settings
 ------------------------------------------------------------------------- */

void PairSPHPi::settings(int narg, char **arg) {
 //if (narg != 0)
   /// error->all(FLERR,
  //      "Illegal number of arguments for pair_style sph/pi");

//Modified by NMC adding viscosity, buklvisco, background press, and epsilon paramter from constructor of the pair style
 if (narg != 3 && narg != 4 )
    error->all (FLERR, "Illegal number of arguments for "
                "pair_style sph/pi ");
  
  viscosity = utils::numeric(FLERR, arg[0], false, lmp);

   /* NMC adding bulk viscosity and background pressure as variable*/
  bulkViscosity = utils::numeric(FLERR, arg[1], false, lmp);
  pb = utils::numeric(FLERR, arg[2], false, lmp);


  if (bulkViscosity < 0) error->all (FLERR, "Bulk viscosity must be positive");
  if (viscosity <= 0) error->all (FLERR, "Viscosity must be positive");
  if (pb <0) error->all (FLERR,"Background pressure must be positive");

  epsilonMM = 1.0;
  micro = 1.0;
  if (narg == 4) epsilonMM = utils::numeric(FLERR, arg[3], false, lmp);
  if (epsilonMM < 0){
    // flag to not use any microscale information
    micro = 0.0;
    epsilonMM = -1*epsilonMM;
  }
  if (epsilonMM > 1 ) error->all (FLERR, "weighting parameter should be between 0 and 1");

}

/* ----------------------------------------------------------------------
 set coeffs for one or more type pairs
 ------------------------------------------------------------------------- */

void PairSPHPi::coeff(int narg, char **arg) {
  if (narg != 5)
    error->all(FLERR,
        "Incorrect args for pair_style sph/pi coefficients");
  if (!allocated)
    allocate();

  int ilo, ihi, jlo, jhi;
  utils::bounds(FLERR,arg[0], 1, atom->ntypes, ilo, ihi, error);
  utils::bounds(FLERR,arg[1], 1, atom->ntypes, jlo, jhi, error);

  double rho0_one = utils::numeric(FLERR,arg[2],false,lmp);
  double soundspeed_one = utils::numeric(FLERR,arg[3],false,lmp);
  //double viscosity_one = utils::numeric(FLERR,arg[4],false,lmp);
  double cut_one = utils::numeric(FLERR,arg[4],false,lmp);
  double B_one = soundspeed_one * soundspeed_one * rho0_one / 7.0;

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    rho0[i] = rho0_one;
    soundspeed[i] = soundspeed_one;
    B[i] = B_one;
    for (int j = MAX(jlo,i); j <= jhi; j++) {
      //viscosity[i][j] = viscosity_one;
      //printf("setting cut[%d][%d] = %f\n", i, j, cut_one);
      cut[i][j] = cut_one;

      setflag[i][j] = 1;
      count++;
    }
  }

  if (count == 0)
    error->all(FLERR,"Incorrect args for pair coefficients");
}

/* ----------------------------------------------------------------------
 init for one type pair i,j and corresponding j,i
 ------------------------------------------------------------------------- */

double PairSPHPi::init_one(int i, int j) {

  if (setflag[i][j] == 0) {
    error->all(FLERR,"All pair sph/pi coeffs are not set");
  }

  cut[j][i] = cut[i][j];
  //viscosity[j][i] = viscosity[i][j];

  return cut[i][j];
}

