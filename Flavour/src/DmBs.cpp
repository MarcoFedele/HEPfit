/*
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "DmBs.h"
#include "StandardModel.h"

double  DmBs::computeThValue()
{
    return(2. * SM.getCBs() * AmpBs(FULLNLO).abs());
}

double  RmBs::computeThValue()
{
    return RBs(FULLNLO).abs() - 1.;
}

double  RmBsm1::computeThValue()
{
    return RBsm1(FULLNLO).real();
}
