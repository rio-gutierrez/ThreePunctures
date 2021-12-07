
#include <iostream>
#include <fstream>
#include <iomanip>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

using namespace std;

extern "C"
void ThreePunctures_Metadata (CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  if (CCTK_MyProc(cctkGH) == 0)
  {
    ofstream o;
    o.open(string(string(out_dir) + "/ThreePunctures.bbh").c_str());

    o << setprecision(19);

    o << "\
# ==================================\n\
# Numerical Relativity Metadata file\n\
# ==================================\n\
#\n\
# This file contains information about the simulation provided by the\n\
# ThreePunctures thorn.  The format is described in the NR Data Format Document\n\
# http://arxiv.org/abs/0709.0093 [draft SVN r707].\n\
" << endl;

    o << "[metadata]" << endl;
    o << "initial-ADM-energy               = " << *E << endl;
    o << "initial-ADM-angular-momentumx    = " << *J1 << endl;
    o << "initial-ADM-angular-momentumy    = " << *J2 << endl;
    o << "initial-ADM-angular-momentumz    = " << *J3 << endl;
    o << "initial-separation of p- and p+  = " << par_b * 2. << endl;
    o << "initial-data-type                = Bowen-York" << endl;
    o << "initial-data-bibtex-keys         = Bowen:1980yu Brandt:1997tf Ansorg:2004ds" << endl;

    o << "initial-bh-position1x            = " << par_b + center_offset[0] << endl;
    o << "initial-bh-position1y            = " << center_offset[1] << endl;
    o << "initial-bh-position1z            = " << center_offset[2] << endl;

    o << "initial-bh-position2x            = " << -par_b + center_offset[0] << endl;
    o << "initial-bh-position2y            = " << center_offset[1] << endl;
    o << "initial-bh-position2z            = " << center_offset[2] << endl;

    o << "initial-bh-position3x            = " << .5 * par_b + center_offset[0] << endl;  // change this according to (*)
    o << "initial-bh-position3y            = " << center_offset[1] << endl;
    o << "initial-bh-position3z            = " << center_offset[2] << endl;

    o << "initial-bh-momentum1x            = " << par_P_plus[0] << endl;
    o << "initial-bh-momentum1y            = " << par_P_plus[1] << endl;
    o << "initial-bh-momentum1z            = " << par_P_plus[2] << endl;

    o << "initial-bh-momentum2x            = " << par_P_minus[0] << endl;
    o << "initial-bh-momentum2y            = " << par_P_minus[1] << endl;
    o << "initial-bh-momentum2z            = " << par_P_minus[2] << endl;

    o << "initial-bh-momentum3x            = " << par_P_third[0] << endl;
    o << "initial-bh-momentum3y            = " << par_P_third[1] << endl;
    o << "initial-bh-momentum3z            = " << par_P_third[2] << endl;

    o << "initial-bh-spin1x                = " << par_S_plus[0] << endl;
    o << "initial-bh-spin1y                = " << par_S_plus[1] << endl;
    o << "initial-bh-spin1z                = " << par_S_plus[2] << endl;

    o << "initial-bh-spin2x                = " << par_S_minus[0] << endl;
    o << "initial-bh-spin2y                = " << par_S_minus[1] << endl;
    o << "initial-bh-spin2z                = " << par_S_minus[2] << endl;

    o << "initial-bh-spin3x                = " << par_S_third[0] << endl;
    o << "initial-bh-spin3y                = " << par_S_third[1] << endl;
    o << "initial-bh-spin3z                = " << par_S_third[2] << endl;

    o << "initial-bh-puncture-adm-mass1    = " << *mp_adm << endl;
    o << "initial-bh-puncture-adm-mass2    = " << *mm_adm << endl;
    o << "initial-bh-puncture-adm-mass3    = " << *mt_adm << endl;

    o << "initial-bh-puncture-bare-mass1   = " << *mp << endl;
    o << "initial-bh-puncture-bare-mass2   = " << *mm << endl;
    o << "initial-bh-puncture-bare-mass3   = " << *mt << endl;
  }
}