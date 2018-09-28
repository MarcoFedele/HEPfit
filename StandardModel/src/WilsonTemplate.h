/* 
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef WILSONTEMPLATE_H
#define	WILSONTEMPLATE_H

#include "OrderScheme.h"
#include <sstream>
#include <stdexcept>
#define MAXORDER_QED QED2
#define MAXORDER NNLO
/**
 * @class WilsonTemplate
 * @ingroup StandardModel
 * @brief A template class for the Wilson coefficients. 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
template <class T> class WilsonTemplate {
public:
    WilsonTemplate(unsigned int dim, schemes scheme_i, orders_qcd order_qcd_i , 
                   orders_qed order_qed_i = QED0)
    {
        size = dim;
        scheme = scheme_i;
        order_qcd = order_qcd_i;
        order_qed = order_qed_i; 
        mu = -1.;
        
        for (int i = LO; i <= MAXORDER; i++) {
            if (i <= order_qcd) elem[i] = new T(size, 0.);
            else elem[i] = NULL;
        }
        
        elem[orders_qed(QED0)] = NULL;
        //for (int i = LO_QED; i <= NLO_QED; i++){
        for (int i = QED0; i <= MAXORDER_QED; i++) {
            if (i <= order_qed) elem[i] = new T(size, 0.);
            else elem[i] = NULL;
        }
    }
    
    WilsonTemplate<T>(const WilsonTemplate<T>& orig) 
    {
        size = orig.size;
        scheme = orig.scheme;
        order_qcd = orig.order_qcd;
        order_qed = orig.order_qed;
        mu = orig.mu;
        for (int i = LO; i <= MAXORDER_QED; i++)
            if (orig.elem[i]!= NULL) elem[i] = new T(*(orig.elem[i]));
            else elem[i] = NULL;
    }
    
    virtual ~WilsonTemplate()
    {
        for (int i = LO; i <= MAXORDER_QED; i++)
            if (elem[i] != NULL) delete elem[i];
    }

    orders_qcd getOrder() const 
    {
        return order_qcd;
    }
    
    orders_qed getOrder_qed() const
    {
        return order_qed;
    }

    double getMu() const
    {
        return mu;
    }

    virtual void resetCoefficient()
    {
        for(int i = LO; i <= order_qcd; i++){
            *(elem[i]) = 0.;
        }
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

    unsigned int  getSize() const
    {
        return size;
    }

protected:
    T* elem[MAXORDER_QED+1];
    unsigned int size;
    double mu;
    schemes scheme;
    orders_qcd order_qcd;
    orders_qed order_qed;

    T * Elem(orders_qcd order) const
    {
        if (order > this->order_qcd) {
            std::stringstream out;
            out << order;
            throw std::runtime_error("WilsonTemplate::getElem(): requested order " + out.str() +
                    " not present in the object");
        }
        return elem[order];
    }
    
    T * Elem(orders_qed order_qed) const
    {
        if ((order_qed > this->order_qed)) {
            std::stringstream out;
            out << order_qed;
            throw std::runtime_error("WilsonTemplate::getElem(): requested order_qed " + out.str() +
                    "not present in the object");
        }
        return elem[order_qed];
    }

    void setElem(const T & v, orders_qcd order_qcd_i)
    {
        if (order_qcd_i > order_qcd) {
            std::stringstream out;
            out << order_qcd_i;
            throw std::runtime_error("MatchingCondition::setElem(): order " + out.str() +
                    " not implemented ");
        }
        *elem[order_qcd_i] = v;
    }
    
    void setElem(const T & v, orders_qed order_qed_i)
    {
        if (order_qed_i > order_qed) {
            std::stringstream out;
            out << order_qed_i;
            throw std::runtime_error("MatchingCondition::setElem(): order " + out.str() +
                    " not implemented ");
        }
        *elem[order_qed_i] = v;
    }
};

#endif	/* WILSONTEMPLATE_H */

