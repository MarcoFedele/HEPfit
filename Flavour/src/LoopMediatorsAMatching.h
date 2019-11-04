/* 
 * Copyright (C) 2018 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef LOOPMEDIATORSAMATCHING_H
#define LOOPMEDIATORSAMATCHING_H

#include "gslpp.h"
#include "StandardModelMatching.h"

class LoopMediatorsA;

/**
 * @class LoopMediatorsA
 * @ingroup FlavourWilsonCoefficient
 * @brief A class for the matching in the NPSMEFTd6.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details  
 */
class LoopMediatorsAMatching : public StandardModelMatching {
public:
    LoopMediatorsAMatching(const LoopMediatorsA & LoopMediatorsA_i);
    
    virtual ~LoopMediatorsAMatching();
    
    /**
     *
     * @brief Updates to new FlavourWilsonCoefficient parameter sets.
     * @return
     */
    
    void updateLoopMediatorsAParameters();
    
    std::vector<WilsonCoefficient>& CMbsg();
    
    std::vector<WilsonCoefficient>& CMprimebsg();
    
    std::vector<WilsonCoefficient>& CMBMll(QCD::lepton lepton);
    
    std::vector<WilsonCoefficient>& CMprimeBMll(QCD::lepton lepton);
    
    std::vector<WilsonCoefficient>& CMdbs2();

private:
    const LoopMediatorsA & myLoopMediatorsA;
    
    gslpp::complex C7NP;
    gslpp::complex C7pNP;
    gslpp::complex C8NP;
    gslpp::complex C8pNP;
    
    gslpp::complex C9NPmu;
    gslpp::complex C9pNPmu;
    gslpp::complex C10NPmu;
    gslpp::complex C10pNPmu;
    gslpp::complex CSNPmu;
    gslpp::complex CSpNPmu;
    gslpp::complex CPNPmu;
    gslpp::complex CPpNPmu;
    
    gslpp::complex C1NP;
    gslpp::complex C2NP;
    gslpp::complex C3NP;
    gslpp::complex C4NP;
    gslpp::complex C5NP;
    gslpp::complex C1pNP;
    gslpp::complex C2pNP;
    gslpp::complex C3pNP;
    
    double WCscale;
    
    WilsonCoefficient mcbsg, mcprimebsg;
    WilsonCoefficient mcBMll, mcprimeBMll;
    WilsonCoefficient mcdbs2;
    std::vector<WilsonCoefficient> vmcbsg, vmcprimebsg;
    std::vector<WilsonCoefficient> vmcBMll, vmcprimeBMll;
    std::vector<WilsonCoefficient> vmcdbs2;
};

#endif /* LOOPMEDIATORSAMATCHING_H */

