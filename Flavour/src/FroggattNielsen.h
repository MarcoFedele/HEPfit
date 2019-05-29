/* 
 * Copyright (C) 2018 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef FROGGATTNIELSEN_H
#define FROGGATTNIELSEN_H

#include "StandardModel.h"
#include "gslpp.h"
#include "ThObservable.h"

/**
 * @class FroggattNielsen
 * @brief Model for NP contributions to Axions.
 */
class FroggattNielsen: public StandardModel {
public:

    static const int NFroggattNielsenvars = 1;

    static const std::string FroggattNielsenvars[NFroggattNielsenvars];
    
    /**
     * @brief FlavourWilsonCoefficient constructor
     */
    FroggattNielsen();
    
    /**
     * @brief FlavourWilsonCoefficient destructor
     */
    ~FroggattNielsen();
    
    virtual bool InitializeModel();
    
    virtual bool Init(const std::map<std::string, double>& DPars);
    
    virtual bool PreUpdate();
    
    virtual bool Update(const std::map<std::string, double>& DPars);
    
    virtual bool PostUpdate();
    
    virtual bool CheckParameters(const std::map<std::string, double>& DPars);
    
    virtual bool setFlag(const std::string name, const bool value);
    
    /*virtual LoopMediatorsMatching& getMatching() const
    {
        return LoopMediatorsM.getObj();
    }*/
    
    /**
     *
     * @return \f$ m_a $\f
     */
    double getma() const {
        return ma;
    }
    
protected: 
    
    virtual void setParameter(const std::string, const double&);
    //mutable Matching<LoopMediatorsMatching,LoopMediators> LoopMediatorsM;
    
    

private:
    
    double ma;
    
      
};



///////////////////////////////////////////////////////////////////////////
// Observables




class test : public ThObservable {
public:
    /**
     * @brief Constructor.
     */
   test(const StandardModel& SM_i);
     
   ~test();

   /**
     * @brief Two positivity conditions of the Higgs potential.
     * @return
     */
    double computeThValue();
    const FroggattNielsen *myFroggattNielsen;
};

#endif /* FROGGATTNIELSEN_H */

