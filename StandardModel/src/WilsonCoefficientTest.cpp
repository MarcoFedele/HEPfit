#include <cstdio>
#include "WilsonCoefficient.h"

int main(void)
{
        std::cout << "QCD" << std::endl;

	WilsonCoefficient coeff(2, NDR, NNLO);
	
        coeff.setMu(5.);
	coeff.setCoeff(0, gslpp::complex(4,2), LO);
	coeff.setCoeff(0, gslpp::complex(1,5), NLO);
	coeff.setCoeff(0, gslpp::complex(3,7), NNLO);
	coeff.setCoeff(1, gslpp::complex(1,8), LO);
	coeff.setCoeff(1, gslpp::complex(5,3), NLO);
	coeff.setCoeff(1, gslpp::complex(6,1), NNLO);

        std::cout << coeff.getCoeffElement(0) << std::endl;
        std::cout << coeff.getCoeff(LO) << std::endl;
        std::cout << coeff.getCoeff() << std::endl;

        std::cout << "QCD+QED" << std::endl;

	WilsonCoefficient coeffe(2, NDR, NLO, QED1);

        coeffe.setMu(5.);
        coeffe.setCoeff(0, gslpp::complex(4,2), LO, QED0);
        coeffe.setCoeff(0, gslpp::complex(1,5), NLO, QED0);
        coeffe.setCoeff(1, gslpp::complex(1,8), LO, QED0);
        coeffe.setCoeff(1, gslpp::complex(5,3), NLO, QED0);
        coeffe.setCoeff(0, gslpp::complex(41,22), LO, QED1);
        coeffe.setCoeff(0, gslpp::complex(11,52), NLO, QED1);
        coeffe.setCoeff(1, gslpp::complex(11,82), LO, QED1);
        coeffe.setCoeff(1, gslpp::complex(51,32), NLO, QED1);


        std::cout << coeffe.getCoeffElement(0) << std::endl;
        std::cout << coeffe.getCoeff(LO, QED1) << std::endl;
        std::cout << coeffe.getCoeff() << std::endl;

	return 0;

}
