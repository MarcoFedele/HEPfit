/*
 * Copyright (C) 2018 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "MyObs.h"
#include "StandardModel.h"

MyObs::MyObs(const StandardModel& SM_i)
: ThObservable(SM_i)
{
};

double MyObs::computeThValue()
{
    return 0.;
}
