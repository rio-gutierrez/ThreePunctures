/* ThreePunctures:  File  "Equations.c"*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include "cctk_Parameters.h"
#include "TP_utilities.h"
#include "ThreePunctures.h"

/* U.d0[ivar]   = U[ivar];  (ivar = 0..nvar-1) */
/* U.d1[ivar]   = U[ivar]_x;  */
/* U.d2[ivar]   = U[ivar]_y;  */
/* U.d3[ivar]   = U[ivar]_z;  */
/* U.d11[ivar]  = U[ivar]_xx; */
/* U.d12[ivar]  = U[ivar]_xy; */
/* U.d13[ivar]  = U[ivar]_xz;*/
/* U.d22[ivar]  = U[ivar]_yy;*/
/* U.d23[ivar]  = U[ivar]_yz;*/
/* U.d33[ivar]  = U[ivar]_zz;*/

CCTK_REAL
BY_KKofxyz (CCTK_REAL x, CCTK_REAL y, CCTK_REAL z)
{
  DECLARE_CCTK_PARAMETERS;
  int i, j;
  CCTK_REAL r_plus, r2_plus, r3_plus,
            r_minus, r2_minus, r3_minus,
            r_third, r2_third, r3_third,
            np_Pp, nm_Pm, nt_Pt,
            Aij, AijAij,
            n_plus[3], n_minus[3], n_third[3],
            np_Sp[3],  nm_Sm[3],   nt_St[3];

  r2_plus  = (x - par_b) * (x - par_b) + y * y + z * z;
  r2_minus = (x + par_b) * (x + par_b) + y * y + z * z;
  r2_third = (x - .5 * par_b) * (x - .5 * par_b) + y * y + z * z;   // change this to whatever is more convenient (*)
  r_plus   = sqrt (r2_plus);
  r_minus  = sqrt (r2_minus);
  r_third  = sqrt (r2_third);
  r3_plus  = r_plus * r2_plus;
  r3_minus = r_minus * r2_minus;
  r3_third = r_third * r2_third;

  n_plus[0]  = (x - par_b) / r_plus;
  n_minus[0] = (x + par_b) / r_minus;
  n_third[0] = (x - .5 * par_b) / r_third;    // change this according to (*)
  n_plus[1]  = y / r_plus;
  n_minus[1] = y / r_minus;
  n_third[1] = y / r_third;
  n_plus[2]  = z / r_plus;
  n_minus[2] = z / r_minus;
  n_third[2] = z / r_third;

  /* dot product: np_Pp = (n_+).(P_+); nm_Pm = (n_-).(P_-); nt_Pt = (n_third).(P_third)  */
  np_Pp = 0;
  nm_Pm = 0;
  nt_Pt = 0;

  for (i = 0; i < 3; i++)
  {
    np_Pp += n_plus[i]  * par_P_plus[i];
    nm_Pm += n_minus[i] * par_P_minus[i];
    nt_Pt += n_third[i] * par_P_third[i];
  }
  /* cross product: np_Sp[i] = [(n_+) x (S_+)]_i; nm_Sm[i] = [(n_-) x (S_-)]_i; nt_St[i] = [(n_third) x (S_third)]_i */
  np_Sp[0] = n_plus[1] * par_S_plus[2] - n_plus[2] * par_S_plus[1];
  np_Sp[1] = n_plus[2] * par_S_plus[0] - n_plus[0] * par_S_plus[2];
  np_Sp[2] = n_plus[0] * par_S_plus[1] - n_plus[1] * par_S_plus[0];
  nm_Sm[0] = n_minus[1] * par_S_minus[2] - n_minus[2] * par_S_minus[1];
  nm_Sm[1] = n_minus[2] * par_S_minus[0] - n_minus[0] * par_S_minus[2];
  nm_Sm[2] = n_minus[0] * par_S_minus[1] - n_minus[1] * par_S_minus[0];
  nt_St[0] = n_third[1] * par_S_third[2] - n_third[2] * par_S_third[1];
  nt_St[1] = n_third[2] * par_S_third[0] - n_third[0] * par_S_third[2];
  nt_St[2] = n_third[0] * par_S_third[1] - n_third[1] * par_S_third[0];

  AijAij = 0;

  for (i = 0; i < 3; i++)
  {
    for (j = 0; j < 3; j++)
    {				/* Bowen-York-Curvature :*/
      Aij =
        + 1.5 * (par_P_plus[i] * n_plus[j] + par_P_plus[j] * n_plus[i]
                      + np_Pp * n_plus[i] * n_plus[j]) / r2_plus
        + 1.5 * (par_P_minus[i] * n_minus[j] + par_P_minus[j] * n_minus[i]
                      + nm_Pm * n_minus[i] * n_minus[j]) / r2_minus
        + 1.5 * (par_P_third[i] * n_third[j] + par_P_third[j] * n_third[i]
                      + nt_Pt * n_third[i] * n_third[j]) / r2_third
        - 3. * (np_Sp[i] * n_plus[j]  + np_Sp[j] * n_plus[i])  / r3_plus
        - 3. * (nm_Sm[i] * n_minus[j] + nm_Sm[j] * n_minus[i]) / r3_minus
        - 3. * (nt_St[i] * n_third[j] + nt_St[j] * n_third[i]) / r3_third;

      if (i == j)
	Aij -= +1.5 * (np_Pp / r2_plus + nm_Pm / r2_minus + nt_Pt / r2_third);
      AijAij += Aij * Aij;
    }
  }

  return AijAij;
}


void
BY_Aijofxyz (CCTK_REAL x, CCTK_REAL y, CCTK_REAL z, CCTK_REAL Aij[3][3])
{
  DECLARE_CCTK_PARAMETERS;
  int i, j;

  CCTK_REAL r_plus, r2_plus, r3_plus,
            r_minus, r2_minus, r3_minus,
            r_third, r2_third, r3_third,
            np_Pp, nm_Pm, nt_Pt,
            n_plus[3], n_minus[3], n_third[3],
            np_Sp[3], nm_Sm[3], nt_St[3];


  r2_plus  = (x - par_b) * (x - par_b) + y * y + z * z;
  r2_minus = (x + par_b) * (x + par_b) + y * y + z * z;
  r2_third = (x - .5 * par_b) * (x - .5 * par_b) + y * y + z * z;   // change this according to (*)

  r2_plus  = sqrt (pow (r2_plus, 2)  + pow (TP_epsilon, 4));
  r2_minus = sqrt (pow (r2_minus, 2) + pow (TP_epsilon, 4));
  r2_third = sqrt (pow (r2_third, 2) + pow (TP_epsilon, 4));

  if (r2_plus < pow(TP_Tiny,2))
    r2_plus = pow(TP_Tiny,2);

  if (r2_minus < pow(TP_Tiny,2))
    r2_minus = pow(TP_Tiny,2);

  if (r2_third < pow(TP_Tiny,2))
    r2_third = pow(TP_Tiny,2);

  r_plus   = sqrt (r2_plus);
  r_minus  = sqrt (r2_minus);
  r_third  = sqrt (r2_third);
  r3_plus  = r_plus * r2_plus;
  r3_minus = r_minus * r2_minus;
  r3_third = r_third * r2_third;

  n_plus[0]  = (x - par_b) / r_plus;
  n_minus[0] = (x + par_b) / r_minus;
  n_third[0] = (x - .5 * par_b) / r_third;   // change this according to (*)
  n_plus[1]  = y / r_plus;
  n_minus[1] = y / r_minus;
  n_third[1] = y / r_third;
  n_plus[2]  = z / r_plus;
  n_minus[2] = z / r_minus;
  n_third[2] = z / r_third;

  /* dot product: np_Pp = (n_+).(P_+); nm_Pm = (n_-).(P_-); nt_Pt = (n_third).(P_third) */
  np_Pp = 0;
  nm_Pm = 0;
  nt_Pt = 0;

  for (i = 0; i < 3; i++)
  {
    np_Pp += n_plus[i]  * par_P_plus[i];
    nm_Pm += n_minus[i] * par_P_minus[i];
    nt_Pt += n_third[i] * par_P_third[i];
  }

  /* cross product: np_Sp[i] = [(n_+) x (S_+)]_i; nm_Sm[i] = [(n_-) x (S_-)]_i; nt_St[i] = [(n_third) x (S_third)]_i */
  np_Sp[0] = n_plus[1] * par_S_plus[2] - n_plus[2] * par_S_plus[1];
  np_Sp[1] = n_plus[2] * par_S_plus[0] - n_plus[0] * par_S_plus[2];
  np_Sp[2] = n_plus[0] * par_S_plus[1] - n_plus[1] * par_S_plus[0];
  nm_Sm[0] = n_minus[1] * par_S_minus[2] - n_minus[2] * par_S_minus[1];
  nm_Sm[1] = n_minus[2] * par_S_minus[0] - n_minus[0] * par_S_minus[2];
  nm_Sm[2] = n_minus[0] * par_S_minus[1] - n_minus[1] * par_S_minus[0];
  nt_St[0] = n_third[1] * par_S_third[2] - n_third[2] * par_S_third[1];
  nt_St[1] = n_third[2] * par_S_third[0] - n_third[0] * par_S_third[2];
  nt_St[2] = n_third[0] * par_S_third[1] - n_third[1] * par_S_third[0];

  for (i = 0; i < 3; i++)
  {
    for (j = 0; j < 3; j++)
    {				/* Bowen-York-Curvature :*/
      Aij[i][j] =
        + 1.5 * (par_P_plus[i] * n_plus[j] + par_P_plus[j] * n_plus[i]
		          + np_Pp * n_plus[i] * n_plus[j]) / r2_plus
        + 1.5 * (par_P_minus[i] * n_minus[j] + par_P_minus[j] * n_minus[i]
              + nm_Pm * n_minus[i] * n_minus[j]) / r2_minus
        + 1.5 * (par_P_third[i] * n_third[j] + par_P_third[j] * n_third[i]
              + nt_Pt * n_third[i] * n_third[j]) / r2_third
        - 3. * (np_Sp[i] * n_plus[j] + np_Sp[j] * n_plus[i]) / r3_plus
        - 3. * (nm_Sm[i] * n_minus[j] + nm_Sm[j] * n_minus[i]) / r3_minus
        - 3. * (nt_St[i] * n_third[j] + nt_St[j] * n_third[i]) / r3_third;

      if (i == j)
	      Aij[i][j] -= +1.5 * (np_Pp / r2_plus + nm_Pm / r2_minus + nt_Pt / r2_third);
    }
  }
}



/*-----------------------------------------------------------*/
/********           Nonlinear Equations                ***********/
/*-----------------------------------------------------------*/
void
NonLinEquations (CCTK_REAL rho_adm,
     CCTK_REAL A, CCTK_REAL B, CCTK_REAL X, CCTK_REAL R,
		 CCTK_REAL x, CCTK_REAL r, CCTK_REAL phi,
		 CCTK_REAL y, CCTK_REAL z, derivs U, CCTK_REAL *values)
{
  DECLARE_CCTK_PARAMETERS;
  CCTK_REAL r_plus, r_minus, r_third, psi, psi2, psi4, psi7;

  r_plus  = sqrt ((x - par_b) * (x - par_b) + y * y + z * z);
  r_minus = sqrt ((x + par_b) * (x + par_b) + y * y + z * z);
  r_third = sqrt ((x - .5 * par_b) * (x - .5 * par_b) + y * y + z * z);    // change this according to (*)

  psi =
    1. + .5 * par_m_plus / r_plus + .5 * par_m_minus / r_minus + .5 * par_m_third / r_third + U.d0[0];
  psi2 = psi * psi;
  psi4 = psi2 * psi2;
  psi7 = psi * psi2 * psi4;

  values[0] =
    U.d11[0] + U.d22[0] + U.d33[0] + 0.125 * BY_KKofxyz (x, y, z) / psi7 +
    2.0 * Pi / psi2/psi * rho_adm;

}



/*-----------------------------------------------------------*/
/********               Linear Equations                ***********/
/*-----------------------------------------------------------*/
void
LinEquations (CCTK_REAL A, CCTK_REAL B, CCTK_REAL X, CCTK_REAL R,
	      CCTK_REAL x, CCTK_REAL r, CCTK_REAL phi,
	      CCTK_REAL y, CCTK_REAL z, derivs dU, derivs U, CCTK_REAL *values)
{
  DECLARE_CCTK_PARAMETERS;
  CCTK_REAL r_plus, r_minus, r_third, psi, psi2, psi4, psi8;

  r_plus  = sqrt ((x - par_b) * (x - par_b) + y * y + z * z);
  r_minus = sqrt ((x + par_b) * (x + par_b) + y * y + z * z);
  r_third = sqrt ((x - .5 * par_b) * (x - .5 * par_b) + y * y + z * z);    // change this according to (*)

  psi =
    1. + .5 * par_m_plus / r_plus + .5 * par_m_minus / r_minus + .5 * par_m_third / r_third + U.d0[0];
  psi2 = psi * psi;
  psi4 = psi2 * psi2;
  psi8 = psi4 * psi4;

  values[0] = dU.d11[0] + dU.d22[0] + dU.d33[0]
    - 0.875 * BY_KKofxyz (x, y, z) / psi8 * dU.d0[0];
}

/*-----------------------------------------------------------*/