/* 
 * File:   EvolDF2.h
 * Author: marco
 *
 * Created on May 20, 2011, 3:55 PMModel
 */

#ifndef EVOLDF2_H
#define	EVOLDF2_H

#include <RGEvolutor.h>
#include <StandardModel.h>
#include <sstream>
#include <gsl/gsl_sf_dilog.h>

class EvolDF2 : public RGEvolutor {
public:
    /**
     * 
     * @brief constructor
     * @param dim 
     * @param scheme
     * @param order
     * @param model 
     */
    EvolDF2(unsigned int dim, schemes scheme, orders order,
            const StandardModel& model);
    
    /**
     * 
     * @brief destructor
     */
    virtual ~EvolDF2();
    
    /**
     * 
     * @param order
     * @param nf number of active flavours
     * @return Anomalous dimension for DeltaF=2 processes
     */
    matrix<double> AnomalousDimension(orders order, unsigned int nf) const;
    
    /**
     * 
     * @param mu low energy scale
     * @param M matching scale
     * @param order
     * @param scheme
     * @return the Wilson coefficients evolved from the scale M to the scale mu
     */
    matrix<double>& Df2Evol(double mu, double M, orders order, 
            schemes scheme = NDR);
    
    /**
     * 
     * @brief Buras et al, hep-ph/9512380
     * @param mu
     * @return the NLO corrective factor for the charm-charm contribution to the kaon oscillations
     */
    double etacc(double mu) const;
    
    /**
     * 
     * @brief Buras et al, hep-ph/9512380
     * @param mu
     * @return the NLO corrective factor for the charm-top contribution to the kaon oscillations
     */
    double etact(double mu) const;
    
    /**
     * 
     * @brief Buras et al, hep-ph/9512380
     * @param mu
     * @return the NLO corrective factor for the top-top contribution to the kaon oscillations
     */
    double etatt(double mu) const;    
    
private:
    double S1tt() const;
    void Df2Evol(double mu, double M, double nf, schemes scheme);
    double a[5], b[5][5][5], c[3][5][5][5], d[3][5][5][5];
    const StandardModel& model;
};

#endif	/* EVOLDF2_H */

