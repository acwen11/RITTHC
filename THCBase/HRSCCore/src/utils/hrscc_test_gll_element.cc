#include <cassert>
#include <cmath>
#include <iostream>

#include <boost/preprocessor/arithmetic/add.hpp>
#include <boost/preprocessor/repetition/repeat_from_to.hpp>

#include <cctk.h>

#include <utils.hh>

#include <hrscc_gll_element.hh>

#define MIN_ELEMENT_ORDER 0
#define EPSILON_D 1e-10
#define EPSILON_I 1e-14
#define EPSILON_T 1e-12

using namespace hrscc;

template<int N>
CCTK_REAL function(CCTK_REAL x) {
    CCTK_REAL out = 0;
    for(int i = 0; i < N; ++i) {
        out += std::pow(x, i);
    }
    return out;
}

template<int N>
CCTK_REAL diff_function(CCTK_REAL x) {
    CCTK_REAL out = 0;
    for(int i = 1; i < N; ++i) {
        out += i*std::pow(x, i-1);
    }
    return out;
}

template<int N>
CCTK_REAL int_function(CCTK_REAL x) {
    CCTK_REAL out = 0;
    for(int i = 0; i < N; ++i) {
        out += 1.0/static_cast<CCTK_REAL>(i+1) * std::pow(x, i+1);
    }
    return out;
}

// We need to do this "by-hand" since we are working in 1D here
template<int N>
void DLT(CCTK_REAL const * func, CCTK_REAL * dlt_func) {
    for(int i = 0; i < N+1; ++i) {
        dlt_func[i] = 0;
        for(int j = 0; j < N+1; ++j) {
            dlt_func[i] += GLLElement<N>::dltop[i][j]*func[j];
        }
    }
}

// We need to do this "by-hand" since we are working in 1D here
template<int N>
void IDLT(CCTK_REAL const * func, CCTK_REAL * idlt_func) {
    for(int i = 0; i < N+1; ++i) {
        idlt_func[i] = 0;
        for(int j = 0; j < N+1; ++j) {
            idlt_func[i] += GLLElement<N>::idltop[i][j]*func[j];
        }
    }
}

int main(void) {
    CCTK_REAL L2error;
    CCTK_REAL Linferror;
    CCTK_REAL interror;
    CCTK_REAL const shape[3] = {2, 2, 2};
    int const stride[3] = {1, 0, 0};

#define TEST_ELEMENT(z, N, unused)                                             \
    {                                                                          \
        CCTK_REAL cfunc[N+1];                                                  \
        CCTK_REAL dfunc[N+1];                                                  \
        CCTK_REAL tfunc[N+1];                                                  \
        CCTK_REAL itfunc[N+1];                                                 \
        CCTK_REAL error[N+1];                                                  \
                                                                               \
        for(int i = 0; i < N+1; ++i) {                                         \
            cfunc[i] = function<N>(GLLElement<N>::node[i]);                    \
        }                                                                      \
                                                                               \
        GLLElement<N>::diff<policy::x>(shape, &cfunc[0], stride,               \
                &dfunc[0], stride);                                            \
                                                                               \
        for(int i = 0; i < N+1; ++i) {                                         \
            error[i] = utils::pow<2>(dfunc[i] -                                \
                    diff_function<N>(GLLElement<N>::node[i]));                 \
        }                                                                      \
                                                                               \
        L2error = std::sqrt(GLLElement<N>::integrate(shape, error, stride));   \
        Linferror = error[0];                                                  \
        for(int i = 1; i < N+1; ++i) {                                         \
            Linferror = Linferror < error[i] ? error[i] : Linferror;           \
        }                                                                      \
        Linferror = std::sqrt(Linferror);                                      \
                                                                               \
        assert(L2error < EPSILON_D);                                           \
        assert(Linferror < EPSILON_D);                                         \
                                                                               \
        CCTK_REAL const inte = int_function<N>(1.0) - int_function<N>(-1.0);   \
        CCTK_REAL intn = 0.0;                                                  \
        /* We need to do this "by-hand" since we are working in 1D here */     \
        for(int i = 0; i < N+1; ++i) {                                         \
            intn += cfunc[i]*GLLElement<N>::weight[i];                         \
        }                                                                      \
                                                                               \
        interror = std::abs(inte - intn);                                      \
                                                                               \
        assert(interror < EPSILON_I);                                          \
                                                                               \
        DLT<N>(cfunc, tfunc);                                                  \
        IDLT<N>(tfunc, itfunc);                                                \
                                                                               \
        for(int i = 0; i < N+1; ++i) {                                         \
            error[i] = utils::pow<2>(itfunc[i] - cfunc[i]);                    \
        }                                                                      \
        L2error = std::sqrt(GLLElement<N>::integrate(shape, error, stride));   \
                                                                               \
        Linferror = error[0];                                                  \
        for(int i = 1; i < N+1; ++i) {                                         \
            Linferror = Linferror < error[i] ? error[i] : Linferror;           \
        }                                                                      \
        Linferror = std::sqrt(Linferror);                                      \
                                                                               \
        std::cout << "N = " << N;                                              \
        std::cout << "\tL2error = " << L2error;                                \
        std::cout << "\tLinferror = " << Linferror;                            \
        std::cout << std::endl;                                                \
                                                                               \
        assert(L2error < EPSILON_T);                                           \
        assert(Linferror < EPSILON_T);                                         \
    }

    BOOST_PP_REPEAT_FROM_TO(MIN_ELEMENT_ORDER,
            BOOST_PP_ADD(HRSCC_GLL_ELEMENT_MAX_ORDER, 1), TEST_ELEMENT, ~);

#undef TEST_ELEMENT
}
