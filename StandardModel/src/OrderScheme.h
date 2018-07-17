/* 
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef ORDERSCHEME_H
#define	ORDERSCHEME_H

//#define MAXORDER FULLNNNLO
//#define MAXORDER_QED FULLNLO_QED

/**
 * @enum schemes
 * @ingroup StandardModel
 * @brief An enum type for regularization schemes.
 */
enum schemes
{
    NDR = 0, /**< Naive dimensional regularization (NDR) scheme */
    HV, /**< 't Hooft-Veltman (HV) scheme */
    LRI /**< Regularization-Independent (RI) renormalization schemes with the Landau gauge */
};

/**
 * @enum orders
 * @ingroup StandardModel
 * @brief An enum type for orders in %QCD.
 */
enum orders_qcd
{
    LO = 0, /**< Leading order */
    NLO, /**< Next-to-leading order */
    NNLO, /**< Next-to-next-to-leading order */
    NNNLO, /**< Next-to-next-to-next-to-leading order */
    NOQCD,
    FULLNLO, /**< Full NLO = LO + NLO */
    FULLNNLO, /**< Full NNLO = LO + NLO + NNLO */
    FULLNNNLO /**< Full NNLO = LO + NLO + NNLO + NNNLO */        
};

/**
 * @enum orders_qed
 * @ingroup StandardModel
 * @brief An enum type for orders in electroweak.
 */
enum orders_qed // WARNING: don't change the ordering, it matters in HeffDF1
{
    NOQED = -1,
    QED0, /* Leading order e/s */
    QED1, /* */
    QED2, /**< Next-to-leading order e */
    FULLQED1, /* all terms up to QED1 included */
    FULLQED2 /* all terms up to QED2 included */
};

#endif	/* ORDERSCHEME_H */