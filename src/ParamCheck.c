/* ThreePunctures:  File  "ThreePunctures.c"*/

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "TP_utilities.h"
#include "ThreePunctures.h"

/* -------------------------------------------------------------------*/
void
ThreePunctures_ParamCheck (CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  if (use_sources)
  {
    if (! CCTK_IsFunctionAliased ("Set_Rho_ADM"))
      CCTK_WARN (0, "Matter sources have been enabled, "
                 "but there is no aliased function for matter sources.");
  }
}
