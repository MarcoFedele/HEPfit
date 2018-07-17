/* 
 * Copyright (C) 2018 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "Expanded.h"

template <> Expanded<complex> Expanded<complex>::conjugate() const {
    uint ord1 = getN1();
    uint ord2 = getN2();
    if (ord2 == 0) {
        std::vector<complex> res(ord1);
        for (uint i = 0; i < ord1; i++)
            res[i] = getOrd(i).conjugate();
        return Expanded<complex>(res);
    } else {
        std::vector<std::vector<complex> > res(ord1, std::vector<complex>(ord2));
        for (uint i = 0; i < ord1; i++)
            for (uint j = 0; j < ord2; j++)
                res[i][j] = getOrd(i, j).conjugate();
        return Expanded<complex>(res);
    }
}

template <> Expanded<matrix<double> > Expanded<matrix<double> >::transpose() const {
        if (ord2 == 0) {
            std::vector<matrix<double> > res(ord1,matrix<double>(data1[0].size_i(),data1[0].size_j(),0.));
            for (uint i = 0; i < ord1; i++)
                    res[i] = data1[i].transpose();
            return Expanded<matrix<double>>(res);
        } else {
            std::vector<std::vector<matrix<double> > > res(ord1, std::vector<matrix<double> >(ord2, matrix<double>(data2[0][0].size_i(),data2[0][0].size_j(),0.)));
            for (uint i1 = 0; i1 < ord1; i1++)
                for (uint i2 = 0; i2 < ord2; i2++)
                            res[i1][i2] = data2[i1][i2].transpose();
            return Expanded<matrix<double>>(res);
        }
    }

template <> Expanded<matrix<complex> > Expanded<matrix<complex> >::transpose() const {
        if (ord2 == 0) {
            std::vector<matrix<complex> > res(ord1,matrix<complex>(data1[0].size_i(),data1[0].size_j(),0.));
            for (uint i = 0; i < ord1; i++)
                    res[i] = data1[i].transpose();
            return Expanded<matrix<complex>>(res);
        } else {
            std::vector<std::vector<matrix<complex> > > res(ord1, std::vector<matrix<complex> >(ord2, matrix<complex>(data2[0][0].size_i(),data2[0][0].size_j(),0.)));
            for (uint i1 = 0; i1 < ord1; i1++)
                for (uint i2 = 0; i2 < ord2; i2++)
                            res[i1][i2] = data2[i1][i2].transpose();
            return Expanded<matrix<complex>>(res);
        }
    }

template <> Expanded<double> Expanded<complex>::real() const {
    uint ord1 = getN1();
    uint ord2 = getN2();
    if (ord2 == 0) {
        std::vector<double> res(ord1, 0.);
        for (uint i = 0; i < ord1; i++)
            res[i] = getOrd(i).real();
        return Expanded<double>(res);
    } else {
        std::vector<std::vector<double> > res(ord1, std::vector<double>(ord2, 0.));
        for (uint i = 0; i < ord1; i++)
            for (uint j = 0; j < ord2; j++)
                res[i][j] = getOrd(i, j).real();
        return Expanded<double>(res);
    }
}

template <> Expanded<double> Expanded<complex>::abs2() const {
    return ((*this) * ((*this).conjugate())).real();
}

template <>
matrix<double> Expanded<matrix<double> >::Zero() const {
    uint id, jd;
    if(ord2 == 0) { 
        id=data1[0].size_i();
        jd=data1[0].size_j();
    } else {
        id=data2[0][0].size_i();
        jd=data2[0][0].size_j();        
    }
    return matrix<double>(id,jd,0.);
}

template <>
matrix<complex> Expanded<matrix<complex> >::Zero() const {
    uint id, jd;
    if(ord2 == 0) { 
        id=data1[0].size_i();
        jd=data1[0].size_j();
    } else {
        id=data2[0][0].size_i();
        jd=data2[0][0].size_j();        
    }
    return matrix<complex>(id,jd,0.);
}

template<>
vector<double> Expanded<vector<double> >::Zero() const {
    uint id;
    if(ord2 == 0) { 
        id=data1[0].size();
    } else {
        id=data2[0][0].size();
    }
    return vector<double>(id,0.);
}

template<>
vector<complex> Expanded<vector<complex> >::Zero() const {
    uint id;
    if(ord2 == 0) { 
        id=data1[0].size();
    } else {
        id=data2[0][0].size();
    }
    return vector<complex>(id,0.);
}


/***************** Begin UnExpanded * Expanded *****************/

// scalar * X -> commute
    template <class T>
    Expanded<T> operator*(const double& ue, const Expanded<T>& ex) {
        return (ex * ue);
    }

    template <class T>
    typename std::enable_if<
      isComp(T), Expanded<T> >::type
     operator*(const complex& ue, const Expanded<T>& ex) {
        return (ex * ue);
    }

    Expanded<complex> operator*(const complex& ue, const Expanded<double>& ex) {
        return (ex * ue);
    }
    Expanded<vector<complex> > operator*(const complex& ue, const Expanded<vector<double> >& ex) {
        return (ex * ue);
    }
    Expanded<matrix<complex> > operator*(const complex& ue, const Expanded<matrix<double> >& ex) {
        return (ex * ue);
    }
    
// X * scalar -> commute

    Expanded<vector<double> > operator*(const vector<double>& ue, const Expanded<double>& ex) {
        return (ex * ue);
    }
    Expanded<vector<complex> > operator*(const vector<double>& ue, const Expanded<complex>& ex) {
        return (ex * ue);
    }
    Expanded<matrix<double> > operator*(const matrix<double>& ue, const Expanded<double>& ex) {
        return (ex * ue);
    }
    Expanded<matrix<complex> > operator*(const matrix<double>& ue, const Expanded<complex>& ex) {
        return (ex * ue);
    }

    template <class T>
    typename std::enable_if<
      isS(T), Expanded<vector<complex> > >::type
    operator*(const vector<complex>& ue, const Expanded<T>& ex) {
        return (ex * ue);
    }
    template <class T>
    typename std::enable_if<
      isS(T), Expanded<matrix<complex> > >::type
    operator*(const matrix<complex>& ue, const Expanded<T>& ex) {
        return (ex * ue);
    }

    // both vector -> commute
    
    Expanded<double> operator*(const vector<double>& ue, const Expanded<vector<double> >& ex) {
        return (ex * ue);
    }
    template <class T, class Q>
    typename std::enable_if<
     isV(T) && isV(Q) && !(isVec(T,double) && isVec(Q,double)), Expanded<complex> >::type
     operator*(const T& ue, const Expanded<Q>& ex) {
        return (ex * ue);
    }

//  vector * matrix

    Expanded<vector<double> > operator*(const vector<double>& ue, const Expanded<matrix<double> >& ex) {
        return (ex.transpose() * ue);
    }

    template <class T, class Q>
    typename std::enable_if<
     isV(T) && isM(Q) && !(isVec(T,double) && isMat(Q,double)), Expanded<vector<complex> > >::type
     operator*(const T& ue, const Expanded<Q>& ex) {
        return (ex.transpose() * ue);
    }

//  matrix * vector

    Expanded<vector<double> > operator*(const matrix<double>& ue, const Expanded<vector<double> >& ex) {
        return (ex * ue.transpose());
    }

    template <class T, class Q>
    typename std::enable_if<
     isM(T) && isV(Q) && !(isMat(T,double) && isVec(Q,double)), Expanded<vector<complex> > >::type
     operator*(const T& ue, const Expanded<Q>& ex) {
        return (ex * ue.transpose());
    }

//  matrix * matrix

    Expanded<matrix<double> > operator*(const matrix<double>& ue, const Expanded<matrix<double> >& ex) {
        return ((ex.transpose() * ue.transpose()).transpose());
    }

    template <class T, class Q>
    typename std::enable_if<
     isM(T) && isM(Q) && !(isMat(T,double) && isMat(Q,double)), Expanded<matrix<complex> > >::type
     operator*(const T& ue, const Expanded<Q>& ex) {
        return ((ex.transpose() * ue.transpose()).transpose());
    }

/***************** End UnExpanded * Expanded *****************/    
