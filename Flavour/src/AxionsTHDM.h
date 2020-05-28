/*
 * Copyright (C) 2018 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef AXIONSTHDM_H
#define AXIONSTHDM_H

#include "StandardModel.h"
#include "gslpp.h"
#include "ThObservable.h"

/**
 * @class Axions
 * @brief Model for NP contributions to Axions.
 */
class AxionsTHDM: public StandardModel {
public:

    static const int NAxionsvars = 23;

    static const std::string Axionsvars[NAxionsvars];

    /**
     * @brief FlavourWilsonCoefficient constructor
     */
    AxionsTHDM();

    /**
     * @brief FlavourWilsonCoefficient destructor
     */
    ~AxionsTHDM();

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

    /**
     *
     * @return \f$ \tan\beta $\f
     */
    double gettanb() const {
        return tanb;
    }

    /**
     *
     * @return \f$ model type $\f
     */
    double getmodel() const {
        return model;
    }

    /**
     *
     * @return \f$ e_L $\f
     */
    double geteps_L() const {
        return eps_L;
    }

    /**
     *
     * @return \f$ E/N $\f
     */
    double getEoN() const {
        return EoN;
    }

    /**
     *
     * @return \f$ E/N $\f
     */
    double getcgamma() const {
        return cgamma;
    }

    /**
     *
     * @return \f$ \Chi_3 $\f
     */
    double getChi3() const {
        return Chi3;
    }

    /**
     *
     * @return \f$ parameter in C_p and C_n $\f
     */
    double getC_pn_err_0() const {
        return C_pn_err_0;
    }

    /**
     *
     * @return \f$ parameter in C_p and C_n $\f
     */
    double getC_pn_err_1() const {
        return C_pn_err_1;
    }

    /**
     *
     * @return \f$ parameter in C_p and C_n $\f
     */
    double getC_pn_err_2() const {
        return C_pn_err_2;
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

    /**
     *
     * @return \f$ g_{a\gamma} $\f
     */
    double gag() const;

    /**
     *
     * @return \f$ C_{ae} $\f
     */
    double Cae() const;

    /**
     *
     * @return \f$ g_{ae} $\f
     */
    double gae() const;

    /**
     *
     * @return \f$ g_{ap} $\f
     */
    double gap() const;

    /**
     *
     * @return \f$ g_{an} $\f
     */
    double gan() const;

protected:

    virtual void setParameter(const std::string, const double&);
    //mutable Matching<LoopMediatorsMatching,LoopMediators> LoopMediatorsM;



private:

    double ma, tanb, sinb, cosb;

    double model, eps_L, EoN, cgamma,Chi3;
            
    double C_pn_err_0, C_pn_err_1, C_pn_err_2;

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




class gagTHDM : public ThObservable {
public:
    /**
     * @brief Constructor.
     */
   gagTHDM(const StandardModel& SM_i);

   ~gagTHDM();

   /**
     * @brief Two positivity conditions of the Higgs potential.
     * @return
     */
    double computeThValue();
    const AxionsTHDM *myAxions;
};

class gaeTHDM : public ThObservable {
public:
    /**
     * @brief Constructor.
     */
   gaeTHDM(const StandardModel& SM_i);

   ~gaeTHDM();

   /**
     * @brief Two positivity conditions of the Higgs potential.
     * @return
     */
    double computeThValue();
    const AxionsTHDM *myAxions;
};

class logtbTHDM : public ThObservable {
public:
    /**
     * @brief Constructor.
     */
   logtbTHDM(const StandardModel& SM_i);

   ~logtbTHDM();

   /**
     * @brief Two positivity conditions of the Higgs potential.
     * @return
     */
    double computeThValue();
    const AxionsTHDM *myAxions;
};

class logmaTHDM : public ThObservable {
public:
    /**
     * @brief Constructor.
     */
   logmaTHDM(const StandardModel& SM_i);

   ~logmaTHDM();

   /**
     * @brief Two positivity conditions of the Higgs potential.
     * @return
     */
    double computeThValue();
    const AxionsTHDM *myAxions;
};

class loggagTHDM : public ThObservable {
public:
    /**
     * @brief Constructor.
     */
   loggagTHDM(const StandardModel& SM_i);

   ~loggagTHDM();

   /**
     * @brief Two positivity conditions of the Higgs potential.
     * @return
     */
    double computeThValue();
    const AxionsTHDM *myAxions;
};

class loggaeTHDM : public ThObservable {
public:
    /**
     * @brief Constructor.
     */
   loggaeTHDM(const StandardModel& SM_i);

   ~loggaeTHDM();

   /**
     * @brief Two positivity conditions of the Higgs potential.
     * @return
     */
    double computeThValue();
    const AxionsTHDM *myAxions;
};

class mac2THDM : public ThObservable {
public:
    /**
     * @brief Constructor.
     */
   mac2THDM(const StandardModel& SM_i);

   ~mac2THDM();

   /**
     * @brief Two positivity conditions of the Higgs potential.
     * @return
     */
    double computeThValue();
    const AxionsTHDM *myAxions;
};

class G117B15ATHDM : public ThObservable {
public:
    /**
     * @brief Constructor.
     */
   G117B15ATHDM(const StandardModel& SM_i);

   ~G117B15ATHDM();

   /**
     * @brief Two positivity conditions of the Higgs potential.
     * @return
     */
    double computeThValue();
    const AxionsTHDM *myAxions;
};

class R548THDM : public ThObservable {
public:
    /**
     * @brief Constructor.
     */
   R548THDM(const StandardModel& SM_i);

   ~R548THDM();

   /**
     * @brief Two positivity conditions of the Higgs potential.
     * @return
     */
    double computeThValue();
    const AxionsTHDM *myAxions;
};

class PG1351489THDM : public ThObservable {
public:
    /**
     * @brief Constructor.
     */
   PG1351489THDM(const StandardModel& SM_i);

   ~PG1351489THDM();

   /**
     * @brief Two positivity conditions of the Higgs potential.
     * @return
     */
    double computeThValue();
    const AxionsTHDM *myAxions;
};

class L113THDM : public ThObservable {
public:
    /**
     * @brief Constructor.
     */
   L113THDM(const StandardModel& SM_i);

   ~L113THDM();

   /**
     * @brief Two positivity conditions of the Higgs potential.
     * @return
     */
    double computeThValue();
    const AxionsTHDM *myAxions;
};

class L192THDM : public ThObservable {
public:
    /**
     * @brief Constructor.
     */
   L192THDM(const StandardModel& SM_i);

   ~L192THDM();

   /**
     * @brief Two positivity conditions of the Higgs potential.
     * @return
     */
    double computeThValue();
    const AxionsTHDM *myAxions;
};

class TRGBTHDM : public ThObservable {
public:
    /**
     * @brief Constructor.
     */
   TRGBTHDM(const StandardModel& SM_i);

   ~TRGBTHDM();

   /**
     * @brief Two positivity conditions of the Higgs potential.
     * @return
     */
    double computeThValue();
    const AxionsTHDM *myAxions;
};

class HBRTHDM : public ThObservable {
public:
    /**
     * @brief Constructor.
     */
   HBRTHDM(const StandardModel& SM_i);

   ~HBRTHDM();

   /**
     * @brief Two positivity conditions of the Higgs potential.
     * @return
     */
    double computeThValue();
    const AxionsTHDM *myAxions;
};

class GaNTHDM : public ThObservable {
public:
    /**
     * @brief Constructor.
     */
   GaNTHDM(const StandardModel& SM_i);

   ~GaNTHDM();

   /**
     * @brief Two positivity conditions of the Higgs potential.
     * @return
     */
    double computeThValue();
    const AxionsTHDM *myAxions;
};

class ganTHDM : public ThObservable {
public:
    /**
     * @brief Constructor.
     */
   ganTHDM(const StandardModel& SM_i);

   ~ganTHDM();

   /**
     * @brief Two positivity conditions of the Higgs potential.
     * @return
     */
    double computeThValue();
    const AxionsTHDM *myAxions;
};

class CaeTHDM : public ThObservable {
public:
    /**
     * @brief Constructor.
     */
   CaeTHDM(const StandardModel& SM_i);

   ~CaeTHDM();

   /**
     * @brief Two positivity conditions of the Higgs potential.
     * @return
     */
    double computeThValue();
    const AxionsTHDM *myAxions;
};

#endif /* AXIONS_H */
