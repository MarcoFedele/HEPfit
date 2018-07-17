/* 
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef WILSONTEMPLATE_H
#define	WILSONTEMPLATE_H

#include "OrderScheme.h"
#include "Expanded.h"
#include <sstream>
#include <stdexcept>

/**
 * @class WilsonTemplate
 * @ingroup StandardModel
 * @brief A template class for the Wilson coefficients. 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details 
 */

// T is the type (double or complex) of the gslpp::vector of Wilson coefficients
template <class T> class WilsonTemplate {
public:
    WilsonTemplate(unsigned int size_i, schemes scheme_i, orders_qcd order_qcd_i, 
                   orders_qed order_qed_i = NOQED)
    {
        gslpp::vector<T> el(size_i, 0.);
        
        size = size_i;
        scheme = scheme_i;
        order_qcd = order_qcd_i;
        order_qed = order_qed_i; 
        mu = -1.;

        if(order_qed_i == NOQED)
        {
            std::vector<gslpp::vector<T> > obj;    
            if(order_qcd < FULLNLO)
                for(uint i = 0; i <= order_qcd; i++)
                    obj.push_back(el);
            else
                throw std::runtime_error("WilsonTemplate::WilsonTemplate(): order_qcd out of range");
           
            coeff = Expanded<gslpp::vector<T> > (obj);
        } else {
            std::vector<std::vector<gslpp::vector<T> > > obj; 
            std::vector<gslpp::vector<T> > tmp;
            if(order_qcd < FULLNLO && order_qed < FULLQED1)
                for(uint i = 0; i <= order_qcd; i++)
                {
                    for(uint j = 0; j <= order_qed; j++) 
                        tmp.push_back(el);
                    obj.push_back(tmp);
                    tmp.clear();
                }
            else
                throw std::runtime_error("WilsonTemplate::WilsonTemplate(): order_qcd and/or order_qed out of range");

            coeff = Expanded<gslpp::vector<T> >(obj);
        }
    };
        

    orders_qcd getOrder_QCD() const 
    {
        return order_qcd;
    }
    
    orders_qed getOrder_QED() const
    {
        return order_qed;
    }

    double getMu() const
    {
        return mu;
    }

    virtual void resetCoefficient()
    {
        gslpp::vector<T> zz(size, 0.);

        if(order_qed == NOQED)
            for(uint i = LO; i <= order_qcd; i++)
              coeff.setOrd(i, zz);
        else
            for(uint i = LO; i <= order_qcd; i++)
                for(uint j = NOQED; j <= order_qed; j++)
              coeff.setOrd(i, j, zz);            
     }

    virtual void setMu(double mu)
    {
        this->mu = mu;
        resetCoefficient();
    }
    
    schemes getScheme() const 
    {
        return scheme;
    }

    void setScheme(schemes scheme)
    {
        this->scheme = scheme;
    }

    unsigned int getSize() const
    {
        return size;
    }

    Expanded<T> getCoeffElement(uint i) const
    {
        Expanded<T> ret;

        if (i >= size) {
            std::stringstream out;
            out << i;
            throw std::runtime_error("WilsonTemplate::getCoeff(): requested element " + out.str() +
                    " not present in the object");
        }
        if(coeff.getN2() == 0)
        {
            std::vector<T> obj;
            for(uint j = 0; j < coeff.getN1(); j++)
                obj.push_back(coeff.getOrd(j)(i));
            ret = Expanded<T>(obj);
            return(ret);
        } else {
            std::vector<std::vector<T> > obj(coeff.getN1());
            for(uint j = 0; j < coeff.getN1(); j++)
                for(uint k = 0; k < coeff.getN2(); k++)                
                    obj[j].push_back(coeff.getOrd(j, k)(i));
            return(Expanded<T>(obj));            
        }
    };

    gslpp::vector<T> getCoeff(orders_qcd order_qcd_i, orders_qed order_qed_i = NOQED) const
    {
        if (order_qcd_i > order_qcd || order_qed_i > order_qed) {
            std::stringstream out;
            out << order_qcd_i << " and " << order_qed_i;
            throw std::runtime_error("WilsonTemplate::getCoeff(): requested order " + out.str() +
                    " not present in the object");
        }
        if(order_qed_i == NOQED)
        {
            if(coeff.getN2() != 0)
                throw std::runtime_error("WilsonTemplate::getCoeff(): wrong arguments");
            return coeff.getOrd(order_qcd_i);
        }
        else
            return coeff.getOrd(order_qcd_i, order_qed_i);
    };
    
    Expanded<gslpp::vector<T> > getCoeff() const
    {
        return coeff;
    }
    
    void setCoeff(const gslpp::vector<T> & v, orders_qcd order_qcd_i, orders_qed order_qed_i = NOQED)
    {
        if (order_qcd_i > order_qcd || order_qed_i > order_qed) {
            std::stringstream out;
            out << order_qcd_i << " and " << order_qed_i;
            throw std::runtime_error("WilsonTemplate::setElem(): order " + out.str() +
                    " not implemented ");
        }
        if (v.size() != size)
            throw std::runtime_error("WilsonTemplate::setElem(): wrong size");

        if(order_qed_i == NOQED)
            coeff.setOrd(order_qcd_i, v);
        else
            coeff.setOrd(order_qcd_i, order_qed_i, v);
    };

    void setCoeff(unsigned int i, T z, orders_qcd order_qcd_i, orders_qed order_qed_i = NOQED) {
        if (i >= size) {
            std::stringstream out;
            out << i;
            throw std::runtime_error("WilsonTemplate::setCoeff(): coefficient index "
                    + out.str() + " out of range");
        }
        if (order_qcd_i > order_qcd || order_qed_i > order_qed) {
            std::stringstream out;
            out << order_qcd_i << " and " << order_qed_i;
            throw std::runtime_error("WilsonTemplate::setCoeff(): order " + out.str() +
                    " not implemented ");
        }
        if (order_qed_i == NOQED) {
            gslpp::vector<gslpp::complex> tmp = coeff.getOrd(order_qcd_i);
            tmp.assign(i, z);
            coeff.setOrd(order_qcd_i, tmp);
        } else {
            gslpp::vector<gslpp::complex> tmp = coeff.getOrd(order_qcd_i, order_qed_i);
            tmp.assign(i, z);
            coeff.setOrd(order_qcd_i, order_qed_i, tmp);

        }
    }

protected:
    Expanded<gslpp::vector<T> > coeff;
    unsigned int size;
    double mu;
    schemes scheme;
    orders_qcd order_qcd;
    orders_qed order_qed;

};

#endif	/* WILSONTEMPLATE_H */

