/* 
 * Copyright (C) 2018 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef MYOBS_H
#define MYOBS_H

class StandardModel;
#include "ThObservable.h"
#include "QCD.h"
#include "OrderScheme.h"

class MyObs : public ThObservable {
public:   
    /**
     * constructor
     * @param Flavour
     */
    MyObs(const StandardModel& SM_i);
    
    /**
     * 
     */
    double computeThValue();
    
    
protected:
    
private:

};

#endif /* MYOBS_H */

