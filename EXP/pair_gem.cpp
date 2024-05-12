// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "pair_gem.h"

#include <cmath>
#include <cstring>
#include "atom.h"
#include "comm.h"
#include "force.h"
#include "neigh_list.h"
#include "math_const.h"
#include "memory.h"
#include "error.h"


using namespace LAMMPS_NS;
using namespace MathConst;

/* ---------------------------------------------------------------------- */

PairGem::PairGem(LAMMPS *lmp) : Pair(lmp)
{
  writedata = 1;
}

/* ---------------------------------------------------------------------- */

PairGem::~PairGem()
{
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);

    memory->destroy(parameter_sig);
    memory->destroy(parameter_eps);
    memory->destroy(parameter_n);
    memory->destroy(cut);
  }
}

/* ---------------------------------------------------------------------- */
void PairGem::compute(int eflag, int vflag)
{
  int i,j,ii,jj,inum,jnum,itype,jtype;
  double xtmp,ytmp,ztmp,delx,dely,delz,evdwl,fpair;
  double r,rsq,arg,factor_lj;
  double red_dist,red_dist_pow,intmed_e;
  int *ilist,*jlist,*numneigh,**firstneigh;

  evdwl = 0.0;
  ev_init(eflag,vflag);

  double **x = atom->x;
  double **f = atom->f;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  double *special_lj = force->special_lj;
  int newton_pair = force->newton_pair;

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
    itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      factor_lj = special_lj[sbmask(j)];
      j &= NEIGHMASK;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;
      jtype = type[j];

      if (rsq < cutsq[itype][jtype]) {
        r = sqrt(rsq);

        red_dist = r / parameter_sig[itype][jtype];
        red_dist_pow = pow(red_dist, parameter_n[itype][jtype]);
        intmed_e = parameter_eps[itype][jtype] * exp(-red_dist_pow);

        if (r > 0.0) fpair = intmed_e * parameter_n[itype][jtype] * red_dist_pow / rsq;
        else fpair = 0.0;
  
        f[i][0] += delx*fpair;
        f[i][1] += dely*fpair;
        f[i][2] += delz*fpair;
        if (newton_pair || j < nlocal) {
          f[j][0] -= delx*fpair;
          f[j][1] -= dely*fpair;
          f[j][2] -= delz*fpair;
        }

        if (eflag)
          evdwl = factor_lj * intmed_e;

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

void PairGem::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag,n+1,n+1,"pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  memory->create(cutsq,n+1,n+1,"pair:cutsq");

  memory->create(parameter_eps,n+1,n+1,"pair:parameter_eps");
  memory->create(parameter_sig,n+1,n+1,"pair:parameter_sig");
  memory->create(parameter_n,n+1,n+1,"pair:parameter_n");
  memory->create(cut,n+1,n+1,"pair:cut");
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairGem::settings(int narg, char **arg)
{
  if (narg != 1) error->all(FLERR,"Illegal pair_style command");

  cut_global = utils::numeric(FLERR,arg[0],false,lmp);

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

void PairGem::coeff(int narg, char **arg)
{
  if (narg < 5 || narg > 6)
    error->all(FLERR,"Incorrect args for pair coefficients");
  if (!allocated) allocate();

  int ilo,ihi,jlo,jhi;
  utils::bounds(FLERR,arg[0],1,atom->ntypes,ilo,ihi,error);
  utils::bounds(FLERR,arg[1],1,atom->ntypes,jlo,jhi,error);

  double parameter_eps_one = utils::numeric(FLERR,arg[2],false,lmp);
  double parameter_sig_one = utils::numeric(FLERR,arg[3],false,lmp);
  double parameter_n_one = utils::numeric(FLERR,arg[4],false,lmp);

  double cut_one = cut_global;
  if (narg == 6) cut_one = utils::numeric(FLERR,arg[5],false,lmp);

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo,i); j <= jhi; j++) {
      parameter_eps[i][j] = parameter_eps_one;
      parameter_sig[i][j] = parameter_sig_one;
      parameter_n[i][j] = parameter_n_one;
      cut[i][j] = cut_one;
      setflag[i][j] = 1;
      count++;
    }
  }

  if (count == 0) error->all(FLERR,"Incorrect args for pair coefficients");
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairGem::init_one(int i, int j)
{
  // always mix prefactors geometrically

  if (setflag[i][j] == 0) {
    parameter_eps[i][j] = sqrt(parameter_eps[i][i]*parameter_eps[j][j]);
    parameter_sig[i][j] = sqrt(parameter_sig[i][i]*parameter_sig[j][j]);
    parameter_n[i][j] = sqrt(parameter_n[i][i]*parameter_n[j][j]);
    cut[i][j] = mix_distance(cut[i][i],cut[j][j]);
  }

  parameter_eps[j][i] = parameter_eps[i][j];
  parameter_sig[j][i] = parameter_sig[i][j];
  parameter_n[j][i] = parameter_n[i][j];
  cut[j][i] = cut[i][j];

  return cut[i][j];
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairGem::write_restart(FILE *fp)
{
  write_restart_settings(fp);

  int i,j;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      fwrite(&setflag[i][j],sizeof(int),1,fp);
      if (setflag[i][j]) {
        fwrite(&parameter_eps[i][j],sizeof(double),1,fp);
        fwrite(&parameter_sig[i][j],sizeof(double),1,fp);
        fwrite(&parameter_n[i][j],sizeof(double),1,fp);
        fwrite(&cut[i][j],sizeof(double),1,fp);
      }
    }
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairGem::read_restart(FILE *fp)
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
          utils::sfread(FLERR,&parameter_eps[i][j],sizeof(double),1,fp,nullptr,error);
          utils::sfread(FLERR,&parameter_sig[i][j],sizeof(double),1,fp,nullptr,error);
          utils::sfread(FLERR,&parameter_n[i][j],sizeof(double),1,fp,nullptr,error);
          utils::sfread(FLERR,&cut[i][j],sizeof(double),1,fp,nullptr,error);
        }
        MPI_Bcast(&parameter_eps[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&parameter_sig[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&parameter_n[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&cut[i][j],1,MPI_DOUBLE,0,world);
      }
    }
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairGem::write_restart_settings(FILE *fp)
{
  fwrite(&cut_global,sizeof(double),1,fp);
  fwrite(&mix_flag,sizeof(int),1,fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairGem::read_restart_settings(FILE *fp)
{
  if (comm->me == 0) {
    utils::sfread(FLERR,&cut_global,sizeof(double),1,fp,nullptr,error);
    utils::sfread(FLERR,&mix_flag,sizeof(int),1,fp,nullptr,error);
  }
  MPI_Bcast(&cut_global,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&mix_flag,1,MPI_INT,0,world);
}

/* ----------------------------------------------------------------------
   proc 0 writes to data file
------------------------------------------------------------------------- */

void PairGem::write_data(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++)
    fprintf(fp,"%d %g %g %g %g\n",i,parameter_eps[i][i],parameter_sig[i][i],parameter_n[i][i],cut[i][i]);
}

/* ----------------------------------------------------------------------
   proc 0 writes all pairs to data file
------------------------------------------------------------------------- */

void PairGem::write_data_all(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++)
    for (int j = i; j <= atom->ntypes; j++)
      fprintf(fp,"%d %d %g %g %g %g\n",i,j,parameter_eps[i][j],parameter_sig[i][j],parameter_n[i][j],cut[i][j]);
}

/* ---------------------------------------------------------------------- */
double PairGem::single(int /*i*/, int /*j*/, int itype, int jtype, double rsq,
                        double /*factor_coul*/, double factor_lj,
                        double &fforce)
{
  double r,philj,red_dist,red_dist_pow;

  r = sqrt(rsq);

  red_dist = r / parameter_sig[itype][jtype];
  red_dist_pow = pow(red_dist, parameter_n[itype][jtype]);

  philj = parameter_eps[itype][jtype] * exp(-red_dist_pow);
  fforce = philj * red_dist_pow / rsq * parameter_n[itype][jtype];

  return factor_lj*philj;
}

/* ---------------------------------------------------------------------- */

void *PairGem::extract(const char *str, int &dim)
{
  dim = 2;
  if (strcmp(str,"eps") == 0) return (void *) parameter_eps;
  if (strcmp(str,"sig") == 0) return (void *) parameter_sig;
  if (strcmp(str,"n") == 0) return (void *) parameter_n;
  return nullptr;
}
