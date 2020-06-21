/* 
 * Copyright (C) 2018 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef AXIONS_H
#define AXIONS_H

#include "StandardModel.h"
#include "gslpp.h"
#include "ThObservable.h"

/**
 * @class Axions
 * @brief Model for NP contributions to Axions.
 */
class Axions: public StandardModel {
public:

    static const int NAxionsvars = 15;

    static const std::string Axionsvars[NAxionsvars];
    
    /**
     * @brief FlavourWilsonCoefficient constructor
     */
    Axions();
    
    /**
     * @brief FlavourWilsonCoefficient destructor
     */
    ~Axions();
    
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
     * @return \f$ m_a \cos^2\beta $\f
     */
    double getgae() const {
        return gae;
    }
    
    /**
     *
     * @return \f$ g_{a\gamma} $\f
     */
    double getgag() const {
        return gag;
    }
    
    /**
     *
     * @return a_G117B15A
     */
    double geta_G117B15A() const {
        return a_G117B15A;
    }
    
    /**
     *
     * @return b_G117B15A
     */
    double getb_G117B15A() const {
        return b_G117B15A;
    }
    
    /**
     *
     * @return a_R548
     */
    double geta_R548() const {
        return a_R548;
    }
    
    /**
     *
     * @return b_R548
     */
    double getb_R548() const {
        return b_R548;
    }
    
    /**
     *
     * @return a_PG1351489
     */
    double geta_PG1351489() const {
        return a_PG1351489;
    }
    
    /**
     *
     * @return b_PG1351489
     */
    double getb_PG1351489() const {
        return b_PG1351489;
    }
    
    /**
     *
     * @return a_L113
     */
    double geta_L113() const {
        return a_L113;
    }
    
    /**
     *
     * @return b_L113
     */
    double getb_L113() const {
        return b_L113;
    }
    
    /**
     *
     * @return a_L192
     */
    double geta_L192() const {
        return a_L192;
    }
    
    /**
     *
     * @return b_L192
     */
    double getb_L192() const {
        return b_L192;
    }
    
    /**
     *
     * @return a_TRGB
     */
    double geta_TRGB() const {
        return a_TRGB;
    }
    
    /**
     *
     * @return b_TRGB
     */
    double getb_TRGB() const {
        return b_TRGB;
    }
    
    /**
     *
     * @return Y_HBR
     */
    double getY_HBR() const {
        return Y_HBR;
    }
    
protected: 
    
    virtual void setParameter(const std::string, const double&);
    //mutable Matching<LoopMediatorsMatching,LoopMediators> LoopMediatorsM;

private:
    
    double gae, gag;
    
    double a_G117B15A, b_G117B15A;
    double a_R548, b_R548;
    double a_PG1351489, b_PG1351489;
    double a_L113, b_L113;
    double a_L192, b_L192;
    
    double a_TRGB, b_TRGB;
    
    double Y_HBR;
    
      
};



///////////////////////////////////////////////////////////////////////////
// Observables




class loggag : public ThObservable {
public:
    /**
     * @brief Constructor.
     */
   loggag(const StandardModel& SM_i);
     
   ~loggag();

   /**
     * @brief Two positivity conditions of the Higgs potential.
     * @return
     */
    double computeThValue();
    const Axions *myAxions;
};

class loggae : public ThObservable {
public:
    /**
     * @brief Constructor.
     */
   loggae(const StandardModel& SM_i);
     
   ~loggae();

   /**
     * @brief Two positivity conditions of the Higgs potential.
     * @return
     */
    double computeThValue();
    const Axions *myAxions;
};

class mac2 : public ThObservable {
public:
    /**
     * @brief Constructor.
     */
   mac2(const StandardModel& SM_i);
     
   ~mac2();

   /**
     * @brief Two positivity conditions of the Higgs potential.
     * @return
     */
    double computeThValue();
    const Axions *myAxions;
};

class G117B15A : public ThObservable {
public:
    /**
     * @brief Constructor.
     */
   G117B15A(const StandardModel& SM_i);
     
   ~G117B15A();

   /**
     * @brief Two positivity conditions of the Higgs potential.
     * @return
     */
    double computeThValue();
    const Axions *myAxions;
};

class R548 : public ThObservable {
public:
    /**
     * @brief Constructor.
     */
   R548(const StandardModel& SM_i);
     
   ~R548();

   /**
     * @brief Two positivity conditions of the Higgs potential.
     * @return
     */
    double computeThValue();
    const Axions *myAxions;
};

class PG1351489 : public ThObservable {
public:
    /**
     * @brief Constructor.
     */
   PG1351489(const StandardModel& SM_i);
     
   ~PG1351489();

   /**
     * @brief Two positivity conditions of the Higgs potential.
     * @return
     */
    double computeThValue();
    const Axions *myAxions;
};

class L113 : public ThObservable {
public:
    /**
     * @brief Constructor.
     */
   L113(const StandardModel& SM_i);
     
   ~L113();

   /**
     * @brief Two positivity conditions of the Higgs potential.
     * @return
     */
    double computeThValue();
    const Axions *myAxions;
};

class L192 : public ThObservable {
public:
    /**
     * @brief Constructor.
     */
   L192(const StandardModel& SM_i);
     
   ~L192();

   /**
     * @brief Two positivity conditions of the Higgs potential.
     * @return
     */
    double computeThValue();
    const Axions *myAxions;
};

class TRGB : public ThObservable {
public:
    /**
     * @brief Constructor.
     */
   TRGB(const StandardModel& SM_i);
     
   ~TRGB();

   /**
     * @brief Two positivity conditions of the Higgs potential.
     * @return
     */
    double computeThValue();
    const Axions *myAxions;
};

class HBR : public ThObservable {
public:
    /**
     * @brief Constructor.
     */
   HBR(const StandardModel& SM_i);
     
   ~HBR();

   /**
     * @brief Two positivity conditions of the Higgs potential.
     * @return
     */
    double computeThValue();
    const Axions *myAxions;
};

class Xenon : public ThObservable {
public:
    /**
     * @brief Constructor.
     */
   Xenon(const StandardModel& SM_i);
     
   ~Xenon();

   /**
     * @brief Two positivity conditions of the Higgs potential.
     * @return
     */
    double computeThValue();
    const Axions *myAxions;
};

#endif /* AXIONS_H */

