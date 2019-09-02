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

    static const int NFroggattNielsenvars = 38;

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
     * @return \f$ \epsilon $\f
     */
    double geteps() const {
        return eps;
    }
    
    /**
     *
     * @return \f$ m_u $\f
     */
    double getmup() const {
        return mup;
    }
    
    /**
     *
     * @return \f$ m_c $\f
     */
    double getmcharm() const {
        return mcharm;
    }
    
    /**
     *
     * @return \f$ m_t $\f
     */
    double getmtop() const {
        return mtop;
    }
    
    /**
     *
     * @return \f$ m_d $\f
     */
    double getmdown() const {
        return mdown;
    }
    
    /**
     *
     * @return \f$ m_s $\f
     */
    double getmstrange() const {
        return mstrange;
    }
    
    /**
     *
     * @return \f$ m_b $\f
     */
    double getmbottom() const {
        return mbottom;
    }
    
    /**
     *
     * @return \f$ s_12^{CKM} $\f
     */
    double gets12CKM() const {
        return s12CKM;
    }
    
    /**
     *
     * @return \f$ c_12^{CKM} $\f
     */
    double getc12CKM() const {
        return c12CKM;
    }
    
    /**
     *
     * @return \f$ s_13^{CKM} $\f
     */
    double gets13CKM() const {
        return s13CKM;
    }
    
    /**
     *
     * @return \f$ c_13^{CKM} $\f
     */
    double getc13CKM() const {
        return c13CKM;
    }
    
    /**
     *
     * @return \f$ s_23^{CKM} $\f
     */
    double gets23CKM() const {
        return s23CKM;
    }
    
    /**
     *
     * @return \f$ c_23^{CKM} $\f
     */
    double getc23CKM() const {
        return c23CKM;
    }
    
    /**
     *
     * @return \f$ \delta_{CKM} $\f
     */
    double getdeltaCKM() const {
        return deltaCKM;
    }
    
    /**
     *
     * @return \f$ Q_{L1} $\f
     */
    double getQL1() const {
        return QL1;
    }
    
    /**
     *
     * @return \f$ Q_{L2} $\f
     */
    double getQL2() const {
        return QL2;
    }
    
    /**
     *
     * @return \f$ Q_{L3} $\f
     */
    double getQL3() const {
        return QL3;
    }
    
    /**
     *
     * @return \f$ Q_{Ru1} $\f
     */
    double getQRu1() const {
        return QRu1;
    }
    
    /**
     *
     * @return \f$ Q_{Ru2} $\f
     */
    double getQRu2() const {
        return QRu2;
    }
    
    /**
     *
     * @return \f$ Q_{Ru3} $\f
     */
    double getQRu3() const {
        return QRu3;
    }
    
    /**
     *
     * @return \f$ Q_{Rd1} $\f
     */
    double getQRd1() const {
        return QRd1;
    }
    
    /**
     *
     * @return \f$ Q_{Rd2} $\f
     */
    double getQRd2() const {
        return QRd2;
    }
    
    /**
     *
     * @return \f$ Q_{Rd3} $\f
     */
    double getQRd3() const {
        return QRd3;
    }
    
    /**
     *
     * @return \f$ th_{L1} $\f
     */
    double getthL1() const {
        return thL1;
    }
    
    /**
     *
     * @return \f$ th_{L2} $\f
     */
    double getthL2() const {
        return thL2;
    }
    
    /**
     *
     * @return \f$ th_{L3} $\f
     */
    double getthL3() const {
        return thL3;
    }
    
    /**
     *
     * @return \f$ phi_{L1} $\f
     */
    double getphiL1() const {
        return phiL1;
    }
    
    /**
     *
     * @return \f$ phi_{L2} $\f
     */
    double getphiL2() const {
        return phiL2;
    }
    
    /**
     *
     * @return \f$ phi_{L3} $\f
     */
    double getphiL3() const {
        return phiL3;
    }
    
    /**
     *
     * @return \f$ th_{Ru1} $\f
     */
    double getthRu1() const {
        return thRu1;
    }
    
    /**
     *
     * @return \f$ th_{Ru2} $\f
     */
    double getthRu2() const {
        return thRu2;
    }
    
    /**
     *
     * @return \f$ th_{Ru3} $\f
     */
    double getthRu3() const {
        return thRu3;
    }
    
    /**
     *
     * @return \f$ phi_{Ru1} $\f
     */
    double getphiRu1() const {
        return phiRu1;
    }
    
    /**
     *
     * @return \f$ phi_{Ru2} $\f
     */
    double getphiRu2() const {
        return phiRu2;
    }
    
    /**
     *
     * @return \f$ phi_{Ru3} $\f
     */
    double getphiRu3() const {
        return phiRu3;
    }
    
    /**
     *
     * @return \f$ th_{Rd1} $\f
     */
    double getthRd1() const {
        return thRd1;
    }
    
    /**
     *
     * @return \f$ th_{Rd2} $\f
     */
    double getthRd2() const {
        return thRd2;
    }
    
    /**
     *
     * @return \f$ th_{Rd3} $\f
     */
    double getthRd3() const {
        return thRd3;
    }
    
    /**
     *
     * @return \f$ phi_{Rd1} $\f
     */
    double getphiRd1() const {
        return phiRd1;
    }
    
    /**
     *
     * @return \f$ phi_{Rd2} $\f
     */
    double getphiRd2() const {
        return phiRd2;
    }
    
    /**
     *
     * @return \f$ phi_{Rd3} $\f
     */
    double getphiRd3() const {
        return phiRd3;
    }
    
    
    ///////////////////
    
    
    /**
     *
     * @return \f$ |CKM| $\f
     */
    gslpp::matrix<gslpp::complex> CKM() const;
    
    /**
     *
     * @return \f$ \lambda_u $\f
     */
    gslpp::matrix<gslpp::complex> lambda_u() const;
    
    /**
     *
     * @return \f$ \lambda_d $\f
     */
    gslpp::matrix<gslpp::complex> lambda_d() const;
    
    /**
     *
     * @return \f$ V $\f
     */
    gslpp::matrix<gslpp::complex> create_V(double th1, double th2, double th3, double phi1, double phi2, double phi3) const;
    
    /**
     *
     * @return \f$ y_u $\f
     */
    gslpp::matrix<gslpp::complex> y_u() const;
    
    /**
     *
     * @return \f$ y_d $\f
     */
    gslpp::matrix<gslpp::complex> y_d() const;
    
    /**
     *
     * @return \f$ F_u $\f
     */
    gslpp::matrix<double> F_u() const;
    
    /**
     *
     * @return \f$ F_d $\f
     */
    gslpp::matrix<double> F_d() const;
    
    /**
     *
     * @return \f$ \delta $\f
     */
    double delta() const;
    
protected: 
    
    virtual void setParameter(const std::string, const double&);
    //mutable Matching<LoopMediatorsMatching,LoopMediators> LoopMediatorsM;
    
    

private:
    
    double eps, mup, mcharm, mtop, mdown, mstrange, mbottom;
    double s12CKM, c12CKM, s13CKM, c13CKM, s23CKM, c23CKM, deltaCKM, JCKM;
    
    double QL1, QL2, QL3;
    double QRu1, QRu2, QRu3;
    double QRd1, QRd2, QRd3;
    
    double thL1, thL2, thL3, phiL1, phiL2, phiL3;
    double thRu1, thRu2, thRu3, phiRu1, phiRu2, phiRu3;
    double thRd1, thRd2, thRd3, phiRd1, phiRd2, phiRd3;
      
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


class Ru11 : public ThObservable {
public:
    /**
     * @brief Constructor.
     */
   Ru11(const StandardModel& SM_i);
     
   ~Ru11();

   /**
     * @brief Two positivity conditions of the Higgs potential.
     * @return
     */
    double computeThValue();
    const FroggattNielsen *myFroggattNielsen;
};


class Ru12 : public ThObservable {
public:
    /**
     * @brief Constructor.
     */
   Ru12(const StandardModel& SM_i);
     
   ~Ru12();

   /**
     * @brief Two positivity conditions of the Higgs potential.
     * @return
     */
    double computeThValue();
    const FroggattNielsen *myFroggattNielsen;
};


class Ru13 : public ThObservable {
public:
    /**
     * @brief Constructor.
     */
   Ru13(const StandardModel& SM_i);
     
   ~Ru13();

   /**
     * @brief Two positivity conditions of the Higgs potential.
     * @return
     */
    double computeThValue();
    const FroggattNielsen *myFroggattNielsen;
};



class Ru21 : public ThObservable {
public:
    /**
     * @brief Constructor.
     */
   Ru21(const StandardModel& SM_i);
     
   ~Ru21();

   /**
     * @brief Two positivity conditions of the Higgs potential.
     * @return
     */
    double computeThValue();
    const FroggattNielsen *myFroggattNielsen;
};


class Ru22 : public ThObservable {
public:
    /**
     * @brief Constructor.
     */
   Ru22(const StandardModel& SM_i);
     
   ~Ru22();

   /**
     * @brief Two positivity conditions of the Higgs potential.
     * @return
     */
    double computeThValue();
    const FroggattNielsen *myFroggattNielsen;
};


class Ru23 : public ThObservable {
public:
    /**
     * @brief Constructor.
     */
   Ru23(const StandardModel& SM_i);
     
   ~Ru23();

   /**
     * @brief Two positivity conditions of the Higgs potential.
     * @return
     */
    double computeThValue();
    const FroggattNielsen *myFroggattNielsen;
};



class Ru31 : public ThObservable {
public:
    /**
     * @brief Constructor.
     */
   Ru31(const StandardModel& SM_i);
     
   ~Ru31();

   /**
     * @brief Two positivity conditions of the Higgs potential.
     * @return
     */
    double computeThValue();
    const FroggattNielsen *myFroggattNielsen;
};


class Ru32 : public ThObservable {
public:
    /**
     * @brief Constructor.
     */
   Ru32(const StandardModel& SM_i);
     
   ~Ru32();

   /**
     * @brief Two positivity conditions of the Higgs potential.
     * @return
     */
    double computeThValue();
    const FroggattNielsen *myFroggattNielsen;
};


class Ru33 : public ThObservable {
public:
    /**
     * @brief Constructor.
     */
   Ru33(const StandardModel& SM_i);
     
   ~Ru33();

   /**
     * @brief Two positivity conditions of the Higgs potential.
     * @return
     */
    double computeThValue();
    const FroggattNielsen *myFroggattNielsen;
};


class Rd11 : public ThObservable {
public:
    /**
     * @brief Constructor.
     */
   Rd11(const StandardModel& SM_i);
     
   ~Rd11();

   /**
     * @brief Two positivity conditions of the Higgs potential.
     * @return
     */
    double computeThValue();
    const FroggattNielsen *myFroggattNielsen;
};


class Rd12 : public ThObservable {
public:
    /**
     * @brief Constructor.
     */
   Rd12(const StandardModel& SM_i);
     
   ~Rd12();

   /**
     * @brief Two positivity conditions of the Higgs potential.
     * @return
     */
    double computeThValue();
    const FroggattNielsen *myFroggattNielsen;
};


class Rd13 : public ThObservable {
public:
    /**
     * @brief Constructor.
     */
   Rd13(const StandardModel& SM_i);
     
   ~Rd13();

   /**
     * @brief Two positivity conditions of the Higgs potential.
     * @return
     */
    double computeThValue();
    const FroggattNielsen *myFroggattNielsen;
};



class Rd21 : public ThObservable {
public:
    /**
     * @brief Constructor.
     */
   Rd21(const StandardModel& SM_i);
     
   ~Rd21();

   /**
     * @brief Two positivity conditions of the Higgs potential.
     * @return
     */
    double computeThValue();
    const FroggattNielsen *myFroggattNielsen;
};


class Rd22 : public ThObservable {
public:
    /**
     * @brief Constructor.
     */
   Rd22(const StandardModel& SM_i);
     
   ~Rd22();

   /**
     * @brief Two positivity conditions of the Higgs potential.
     * @return
     */
    double computeThValue();
    const FroggattNielsen *myFroggattNielsen;
};


class Rd23 : public ThObservable {
public:
    /**
     * @brief Constructor.
     */
   Rd23(const StandardModel& SM_i);
     
   ~Rd23();

   /**
     * @brief Two positivity conditions of the Higgs potential.
     * @return
     */
    double computeThValue();
    const FroggattNielsen *myFroggattNielsen;
};



class Rd31 : public ThObservable {
public:
    /**
     * @brief Constructor.
     */
   Rd31(const StandardModel& SM_i);
     
   ~Rd31();

   /**
     * @brief Two positivity conditions of the Higgs potential.
     * @return
     */
    double computeThValue();
    const FroggattNielsen *myFroggattNielsen;
};


class Rd32 : public ThObservable {
public:
    /**
     * @brief Constructor.
     */
   Rd32(const StandardModel& SM_i);
     
   ~Rd32();

   /**
     * @brief Two positivity conditions of the Higgs potential.
     * @return
     */
    double computeThValue();
    const FroggattNielsen *myFroggattNielsen;
};


class Rd33 : public ThObservable {
public:
    /**
     * @brief Constructor.
     */
   Rd33(const StandardModel& SM_i);
     
   ~Rd33();

   /**
     * @brief Two positivity conditions of the Higgs potential.
     * @return
     */
    double computeThValue();
    const FroggattNielsen *myFroggattNielsen;
};


class FN_delta : public ThObservable {
public:
    /**
     * @brief Constructor.
     */
   FN_delta(const StandardModel& SM_i);
     
   ~FN_delta();

   /**
     * @brief Two positivity conditions of the Higgs potential.
     * @return
     */
    double computeThValue();
    const FroggattNielsen *myFroggattNielsen;
};

#endif /* FROGGATTNIELSEN_H */

