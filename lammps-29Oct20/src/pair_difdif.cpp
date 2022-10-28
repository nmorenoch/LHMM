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

/* ----------------------------------------------------------------------
   Contributing author: Kurt Smith (U Pittsburgh)

   Modified version of traditional DIFDIF with for diffusive diffusivity problems
   Modified by Nicolas Moreno @ nmoreno@bcamath.org
------------------------------------------------------------------------- */

#include "pair_difdif.h"

#include <cmath>
#include "atom.h"
#include "comm.h"
#include "update.h"
#include "force.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "random_mars.h"
#include "memory.h"
#include "error.h"


using namespace LAMMPS_NS;

#define EPSILON 1.0e-10

/* ---------------------------------------------------------------------- */

PairDIFDIF::PairDIFDIF(LAMMPS *lmp) : Pair(lmp)
{
  writedata = 1;
  random = nullptr;
}

/* ---------------------------------------------------------------------- */

PairDIFDIF::~PairDIFDIF()
{
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);
    memory->destroy(cut);
    memory->destroy(a0);
    memory->destroy(gamma);
    memory->destroy(sigma);
  }

  if (random) delete random;
}

/* ---------------------------------------------------------------------- */

void PairDIFDIF::compute(int eflag, int vflag)
{
  int i,j,ii,jj,inum,jnum,itype,jtype;
  double xtmp,ytmp,ztmp,delx,dely,delz,evdwl,fpair;
  double vxtmp,vytmp,vztmp,delvx,delvy,delvz;
  double rsq,r,rinv,dot,wd,randnum,factor_DIFDIF;
  int *ilist,*jlist,*numneigh,**firstneigh;


  //New variables for modified friction coefficient and corresponding noise.
  double sigvar, gammaij;
  

  evdwl = 0.0;
  ev_init(eflag,vflag);

  double *gammai = atom->gammai;
  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  double *special_lj = force->special_lj;
  int newton_pair = force->newton_pair;
  double dtinvsqrt = 1.0/sqrt(update->dt);

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

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

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      factor_DIFDIF = special_lj[sbmask(j)];
      j &= NEIGHMASK;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;
      jtype = type[j];

      if (rsq < cutsq[itype][jtype]) {
        r = sqrt(rsq);
        if (r < EPSILON) continue;     // r can be 0.0 in DIFDIF systems
        rinv = 1.0/r;
        delvx = vxtmp - v[j][0];
        delvy = vytmp - v[j][1];
        delvz = vztmp - v[j][2];
        dot = delx*delvx + dely*delvy + delz*delvz;
        wd = 1.0 - r/cut[itype][jtype];
        randnum = random->gaussian();

        // conservative force = a0 * wd
        // drag force = -gamma * wd^2 * (delx dot delv) / r
        // random force = sigma * wd * rnd * dtinvsqrt;

        fpair = a0[itype][jtype]*wd;
        
        //old
        //fpair -= gamma[itype][jtype]*wd*wd*dot*rinv;
        //fpair += sigma[itype][jtype]*wd*randnum*dtinvsqrt;
        
        //new
        gammaij=0.5*(gammai[i]+gammai[j]);  //Mean friction between particles     

        sigvar = sqrt(2.0*force->boltz*temperature*gammaij); //sigma[itype][jtype]
            

        fpair -= gammaij*wd*wd*dot*rinv;
        fpair += sigvar*wd*randnum*dtinvsqrt;
        //endnew
        
        fpair *= factor_DIFDIF*rinv;
        
        f[i][0] += delx*fpair;
        f[i][1] += dely*fpair;
        f[i][2] += delz*fpair;
        if (newton_pair || j < nlocal) {
          f[j][0] -= delx*fpair;
          f[j][1] -= dely*fpair;
          f[j][2] -= delz*fpair;
        }

        if (eflag) {
          // unshifted eng of conservative term:
          // evdwl = -a0[itype][jtype]*r * (1.0-0.5*r/cut[itype][jtype]);
          // eng shifted to 0.0 at cutoff
          evdwl = 0.5*a0[itype][jtype]*cut[itype][jtype] * wd*wd;
          evdwl *= factor_DIFDIF;
        }

        if (evflag) ev_tally(i,j,nlocal,newton_pair,
                             evdwl,0.0,fpair,delx,dely,delz);
      }
    }
  }

  if (vflag_fdotr) virial_fdotr_compute();
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairDIFDIF::allocate()
{
  int i,j;
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag,n+1,n+1,"pair:setflag");
  for (i = 1; i <= n; i++)
    for (j = i; j <= n; j++)
      setflag[i][j] = 0;

  memory->create(cutsq,n+1,n+1,"pair:cutsq");

  memory->create(cut,n+1,n+1,"pair:cut");
  memory->create(a0,n+1,n+1,"pair:a0");
  memory->create(gamma,n+1,n+1,"pair:gamma");
  memory->create(sigma,n+1,n+1,"pair:sigma");
  for (i = 0; i <= atom->ntypes; i++)
    for (j = 0; j <= atom->ntypes; j++)
      sigma[i][j] = gamma[i][j] = 0.0;
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairDIFDIF::settings(int narg, char **arg)
{
  if (narg != 6) error->all(FLERR,"Illegal pair_style command");

  temperature = utils::numeric(FLERR,arg[0],false,lmp);
  cut_global = utils::numeric(FLERR,arg[1],false,lmp);
  seed = utils::inumeric(FLERR,arg[2],false,lmp);
  
  disType = utils::inumeric(FLERR,arg[3],false,lmp);
  
  ///reading the mean or expected value.
  double gmean  = utils::numeric(FLERR,arg[4],false,lmp);
  double sGamma  = utils::numeric(FLERR,arg[5],false,lmp);
  //Modified for distType variable. So far only 5 types if more added then increases the conditional for distype
  if (disType>7) error->all(FLERR, "Illegal number of distribution type");  


  // initialize Marsaglia RNG with processor-unique seed
  if (seed <= 0) error->all(FLERR,"Illegal pair_style command");
  delete random;
  random = new RanMars(lmp,seed + comm->me);

  double *gammai = atom->gammai;
  ///Setting a distributed gamma for each particle in the procs.
  for (int ii = 0; ii<atom->nlocal;ii++){
    if (disType==0) gammai[ii] = gmean;
    if (disType==1) gammai[ii] = random->Halfgaussian(gmean,sGamma);
    if (disType==2) gammai[ii] = random->lognormal(gmean, sGamma);
    if (disType==3) gammai[ii] = random->onesideLevy(1.0, gmean);
    if (disType==4) gammai[ii] = random->maxboltz(gmean);
    if (disType==5) gammai[ii] = random->gengamma(gmean/1.3, 1.3, sGamma);

  }

// reset cutoffs that have been explicitly set
  if (allocated) {
    int i,j;
    for (i = 1; i <= atom->ntypes; i++)
      for (j = i; j <= atom->ntypes; j++)
        if (setflag[i][j]) cut[i][j] = cut_global;
  }
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairDIFDIF::coeff(int narg, char **arg)
{
  if (narg < 4 || narg > 5)
    error->all(FLERR,"Incorrect args for pair coefficients");
  if (!allocated) allocate();

  int ilo,ihi,jlo,jhi;
  utils::bounds(FLERR,arg[0],1,atom->ntypes,ilo,ihi,error);
  utils::bounds(FLERR,arg[1],1,atom->ntypes,jlo,jhi,error);

  double a0_one = utils::numeric(FLERR,arg[2],false,lmp);
  double gamma_one = utils::numeric(FLERR,arg[3],false,lmp);

  double cut_one = cut_global;
  if (narg == 5) cut_one = utils::numeric(FLERR,arg[4],false,lmp);

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo,i); j <= jhi; j++) {
      a0[i][j] = a0_one;

      //gamma[i][j] = gamma_one;
      cut[i][j] = cut_one;
      setflag[i][j] = 1;
      count++;
    }

  }
  
  if (count == 0) error->all(FLERR,"Incorrect args for pair coefficients");
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairDIFDIF::init_style()
{
  if (comm->ghost_velocity == 0)
    error->all(FLERR,"Pair DIFDIF requires ghost atoms store velocity");

  // if newton off, forces between atoms ij will be double computed
  // using different random numbers

  if (force->newton_pair == 0 && comm->me == 0) error->warning(FLERR,
      "Pair DIFDIF needs newton pair on for momentum conservation");

  neighbor->request(this,instance_me);
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairDIFDIF::init_one(int i, int j)
{
  if (setflag[i][j] == 0) error->all(FLERR,"All pair coeffs are not set");

  sigma[i][j] = sqrt(2.0*force->boltz*temperature*gamma[i][j]);

  cut[j][i] = cut[i][j];
  a0[j][i] = a0[i][j];
  gamma[j][i] = gamma[i][j];
  sigma[j][i] = sigma[i][j];

  return cut[i][j];
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairDIFDIF::write_restart(FILE *fp)
{
  write_restart_settings(fp);

  int i,j;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      fwrite(&setflag[i][j],sizeof(int),1,fp);
      if (setflag[i][j]) {
        fwrite(&a0[i][j],sizeof(double),1,fp);
        fwrite(&gamma[i][j],sizeof(double),1,fp);
        fwrite(&cut[i][j],sizeof(double),1,fp);
      }
    }
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairDIFDIF::read_restart(FILE *fp)
{
  read_restart_settings(fp);

  allocate();

  int i,j;
  int me = comm->me;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      if (me == 0) utils::sfread(FLERR,&setflag[i][j],sizeof(int),1,fp,nullptr,error);
      MPI_Bcast(&setflag[i][j],1,MPI_INT,0,world);
      if (setflag[i][j]) {
        if (me == 0) {
          utils::sfread(FLERR,&a0[i][j],sizeof(double),1,fp,nullptr,error);
          utils::sfread(FLERR,&gamma[i][j],sizeof(double),1,fp,nullptr,error);
          utils::sfread(FLERR,&cut[i][j],sizeof(double),1,fp,nullptr,error);
        }
        MPI_Bcast(&a0[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&gamma[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&cut[i][j],1,MPI_DOUBLE,0,world);
      }
    }
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairDIFDIF::write_restart_settings(FILE *fp)
{
  fwrite(&temperature,sizeof(double),1,fp);
  fwrite(&cut_global,sizeof(double),1,fp);
  fwrite(&seed,sizeof(int),1,fp);
  fwrite(&mix_flag,sizeof(int),1,fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairDIFDIF::read_restart_settings(FILE *fp)
{
  if (comm->me == 0) {
    utils::sfread(FLERR,&temperature,sizeof(double),1,fp,nullptr,error);
    utils::sfread(FLERR,&cut_global,sizeof(double),1,fp,nullptr,error);
    utils::sfread(FLERR,&seed,sizeof(int),1,fp,nullptr,error);
    utils::sfread(FLERR,&mix_flag,sizeof(int),1,fp,nullptr,error);
  }
  MPI_Bcast(&temperature,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&cut_global,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&seed,1,MPI_INT,0,world);
  MPI_Bcast(&mix_flag,1,MPI_INT,0,world);

  // initialize Marsaglia RNG with processor-unique seed
  // same seed that pair_style command initially specified

  if (random) delete random;
  random = new RanMars(lmp,seed + comm->me);
}

/* ----------------------------------------------------------------------
   proc 0 writes to data file
------------------------------------------------------------------------- */

void PairDIFDIF::write_data(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++)
    fprintf(fp,"%d %g %g\n",i,a0[i][i],gamma[i][i]);
}

/* ----------------------------------------------------------------------
   proc 0 writes all pairs to data file
------------------------------------------------------------------------- */

void PairDIFDIF::write_data_all(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++)
    for (int j = i; j <= atom->ntypes; j++)
      fprintf(fp,"%d %d %g %g %g\n",i,j,a0[i][j],gamma[i][j],cut[i][j]);
}

/* ---------------------------------------------------------------------- */

double PairDIFDIF::single(int /*i*/, int /*j*/, int itype, int jtype, double rsq,
                       double /*factor_coul*/, double factor_DIFDIF, double &fforce)
{
  double r,rinv,wd,phi;

  r = sqrt(rsq);
  if (r < EPSILON) {
    fforce = 0.0;
    return 0.0;
  }

  rinv = 1.0/r;
  wd = 1.0 - r/cut[itype][jtype];
  fforce = a0[itype][jtype]*wd * factor_DIFDIF*rinv;

  phi = 0.5*a0[itype][jtype]*cut[itype][jtype] * wd*wd;
  return factor_DIFDIF*phi;
}
