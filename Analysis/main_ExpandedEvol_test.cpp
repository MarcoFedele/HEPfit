/*
 * Copyright (C) 2014 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

/**
 *  To check the two codes one needs to:
 *  1. put the same sw2 and add the missinig matching (to check penguins) in StandardModelMatching
 *  2. replace the Rest function with the large Mt-expression (Rest is computed for the GF^2 normalization)
 *  3. Use ComputeCoeffsmumuStandardNorm()
 */

#include <iostream>
#include <ComputeObservables.h>
#include "HeffDF1.h"

int main(void) {

    /* Define the model configuration file.                        */
    /* Here it is passed as the first argument to the executable.  */
    /* The model configuration file provides the default values of */
    /* the mandatory model parameters.                             */
    std::string ModelConf = "StandardModel.conf";

    /* Define a map for the parameters to be varied. */
    std::map<std::string, double> DPars;

    /* Create objects of the classes ModelFactory and ThObsFactory */
    ModelFactory ModelF;
    ThObsFactory ThObsF;

    /* register user-defined model named ModelName defined in class ModelClass using the following syntax: */
    /* ModelF.addModelToFactory(ModelName, boost::factory<ModelClass*>() ) */

    /* register user-defined ThObservable named ThObsName defined in class ThObsClass using the following syntax: */
    /* ThObsF.addObsToFactory(ThObsName, boost::factory<ThObsClass*>() )*/

    /* Create an object of the class ComputeObservables. */
    ComputeObservables CO(ModelF, ThObsF, ModelConf);
    StandardModel& mySM = *CO.getModel();

    double mub = 5.;
    double muW = 120.;
    
    double myMW = mySM.Mw();

    HeffDF1 Heff("CPMLQB", mySM, NNLO);

    std::cout << "%SUITE_STARTING% Evolutor" << std::endl;
    std::cout << "%SUITE_STARTED% *****" << std::endl;
    Expanded<gslpp::vector<gslpp::complex> > allcoeff(Heff.ComputeCoeff(mub));
    Expanded<gslpp::vector<gslpp::complex> > total(allcoeff);

    for(double mu=8;mu<9;mu+=0.005) {
      allcoeff = Heff.ComputeCoeff(mu);
      total = total +allcoeff;
    }
    std::cout << total << std::endl;
      

//    std::cout << std::endl << "00:" << std::endl;
//    std::cout << Heff.LowScaleCoeff(00) <<  std::endl;
//
//    std::cout << std::endl << "10:" << std::endl;
//    std::cout << Heff.LowScaleCoeff(10) <<  std::endl;
//
//    std::cout << std::endl << "20:" << std::endl;
//    std::cout << Heff.LowScaleCoeff(20) <<  std::endl;
//      
//    std::cout << std::endl << "01:" << std::endl;
//    std::cout << Heff.LowScaleCoeff(01) <<  std::endl;
//
//    std::cout << std::endl << "11:" << std::endl;
//    std::cout << Heff.LowScaleCoeff(11) <<  std::endl;
//
//    std::cout << std::endl << "21:" << std::endl;
//    std::cout << Heff.LowScaleCoeff(21) <<  std::endl;
//      
//    std::cout << std::endl << "02:" << std::endl;
//    std::cout << Heff.LowScaleCoeff(02) <<  std::endl;
//
//    std::cout << std::endl << "12:" << std::endl;
//    std::cout << Heff.LowScaleCoeff(12) <<  std::endl;
//
//    std::cout << std::endl << "22:" << std::endl;
//    std::cout << Heff.LowScaleCoeff(22) <<  std::endl;
      
}