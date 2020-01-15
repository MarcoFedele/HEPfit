/* 
 * Copyright (C) 2018 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef LOOPMEDIATORSDMMATCHING_H
#define LOOPMEDIATORSDMMATCHING_H

#include "gslpp.h"
#include "StandardModelMatching.h"

class LoopMediatorsDM;

/**
 * @class LoopMediatorsDM
 * @ingroup FlavourWilsonCoefficient
 * @brief A class for the matching in the NPSMEFTd6.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details  
 */
class LoopMediatorsDMMatching : public StandardModelMatching {
public:
    LoopMediatorsDMMatching(const LoopMediatorsDM & LoopMediatorsDM_i);
    
    virtual ~LoopMediatorsDMMatching();
    
    /**
     *
     * @brief Updates to new FlavourWilsonCoefficient parameter sets.
     * @return
     */
    
    void updateLoopMediatorsDMParameters();
    std::vector<WilsonCoefficient>& CMBMll(QCD::lepton lepton);

private:
    const LoopMediatorsDM & myLoopMediatorsDM;
    
    gslpp::complex C9NPmu;
    gslpp::complex C10NPmu;
    
    WilsonCoefficient mcBMll;
    std::vector<WilsonCoefficient> vmcBMll;
};

#endif /* LOOPMEDIATORSDMMATCHING_H */

