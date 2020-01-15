/*
 * Copyright (C) 2018 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef LOOPMEDIATORSDM_H
#define LOOPMEDIATORSDM_H

#include "StandardModel.h"
#include "gslpp.h"
#include "LoopMediatorsDMMatching.h"
#include "ThObservable.h"

/**
 * @class LoopMediatorsDM
 * @brief Model for NP contributions to flavour.
 */
class LoopMediatorsDM: public StandardModel {
public:

    static const int NLoopMediatorsDMvars = 7;

    static const std::string LoopMediatorsDMvars[NLoopMediatorsDMvars];

    /**
     * @brief FlavourWilsonCoefficient constructor
     */
    LoopMediatorsDM();

    /**
     * @brief FlavourWilsonCoefficient destructor
     */
    ~LoopMediatorsDM();

    virtual bool InitializeModel();

    virtual bool Init(const std::map<std::string, double>& DPars);

    virtual bool PreUpdate();

    virtual bool Update(const std::map<std::string, double>& DPars);

    virtual bool PostUpdate();

    virtual bool CheckParameters(const std::map<std::string, double>& DPars);

    virtual bool setFlag(const std::string name, const bool value);

    double F9(double x, double y);
    double G9(double x, double y);

    virtual LoopMediatorsDMMatching& getMatching() const
    {
        return LoopMediatorsDMM.getObj();
    }

    /**
     *
     * @return \f$ C_9 $\f
     */
    double getC9() const {
        return C9;
    }

    /**
     *
     * @return \f$ C_10 $\f
     */
    double getC10() const {
        return C10;
    }

protected:

    virtual void setParameter(const std::string, const double&);
    mutable Matching<LoopMediatorsDMMatching,LoopMediatorsDM> LoopMediatorsDMM;

private:

    double C9;
    double C10;

    double ysybgD;
    double gmuV;
    double gmuA;
    double QB;
    double mB_NP;
    double mchi_NP;
    double mV_NP;
    
};



class C9_LMD : public ThObservable {
public:
    /**
     * @brief Constructor.
     */
   C9_LMD(const StandardModel& SM_i);

   ~C9_LMD();

   /**
     * @brief Two positivity conditions of the Higgs potential.
     * @return
     */
    double computeThValue();
    const LoopMediatorsDM * myLM;
};



class C10_LMD : public ThObservable {
public:
    /**
     * @brief Constructor.
     */
   C10_LMD(const StandardModel& SM_i);

   ~C10_LMD();

   /**
     * @brief Two positivity conditions of the Higgs potential.
     * @return
     */
    double computeThValue();
    const LoopMediatorsDM * myLM;
};

#endif /* LOOPMEDIATORSDM_H */
