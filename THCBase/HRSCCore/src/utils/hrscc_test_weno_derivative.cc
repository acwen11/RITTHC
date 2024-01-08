#include <algorithm>
#include <cassert>
#include <fstream>
#include <iostream>
#include <string>

#include <cctk.h>

#include <hrscc_config_par.hh>
#include <hrscc_eno_stencil.hh>
#include <hrscc_fd_reconstruction.hh>
#include <hrscc_typedefs.hh>
#include <hrscc_weno_limiter.hh>
#include <hrscc_weno_stencil.hh>
#include <hrscc_weno_weights.hh>
#include <hrscc_weno_reconstruction.hh>

///////////////////////////////////////////////////////////////////////////////
// Configuration
///////////////////////////////////////////////////////////////////////////////
// Type of WENO reconstruction
#define config_eno_width            3
#define config_weno_type            hrscc::policy::standard
#define config_weno_optim           hrscc::policy::order
#define config_weno_limit           hrscc::policy::dummy

// Parameters for the WENO reconstruction
#define config_weno_alpha           10.0
#define config_weno_epsilon         1.0e-5

// Grid
#define config_x0                   0.0
#define config_x1                   1.0
#define config_n0                   100
#define config_ngrids               5
#define config_grid_factor          2

// Test function
#define config_test_function        x*utils::pow<2>(std::sin(4*M_PI*x))
#define config_diff_test_function   std::sin(4*M_PI*x)*(std::sin(4*M_PI*x) +   \
                                        8*M_PI*x*std::cos(4*M_PI*x));

// Define to output the data onto files
//#define config_output_full_data
// Define to save the Linf errors
//#define config_output_error
// Verbose
#define config_verbose

// Reference value for the convergence order
#define config_check_order
#define config_avg_order_reference  5.15311
#define config_toll                 1e-3

///////////////////////////////////////////////////////////////////////////////

#define my_eno_stencil              hrscc::ENOStencil<config_eno_width>
#define my_weno_weights             hrscc::WENOWeights<config_eno_width>
#define my_weno_stencil             hrscc::WENOStencil<config_eno_width,       \
                                        config_weno_type, config_weno_optim>
#define my_weno_limiter             hrscc::WENOLimiter<config_eno_width,       \
                                        my_weno_stencil::width,                \
                                        config_weno_limit>
#define my_weno_reconstruction      hrscc::WENOReconstruction<my_eno_stencil,  \
                                        my_weno_limiter, my_weno_stencil,      \
                                        my_weno_weights >

#define nghostzones                 my_weno_reconstruction::width
///////////////////////////////////////////////////////////////////////////////

typedef CCTK_REAL (* function)(CCTK_REAL);

class test_function {
    public:
        CCTK_REAL operator()(CCTK_REAL x) const {
            return config_test_function;
        }
};

class diff_test_function {
    public:
        CCTK_REAL operator()(CCTK_REAL x) const {
            return config_diff_test_function;
        }
};

CCTK_REAL * make_grid(CCTK_REAL a, CCTK_REAL b, int npoints,
        bool interfaces = false) {
    CCTK_REAL Delta = (b - a) / (npoints - 1);
    int size = npoints + 2 - (2*nghostzones - 1)*interfaces;
    CCTK_REAL * out = new CCTK_REAL[size];

    out[0] = a - Delta + (nghostzones - 0.5)*interfaces*Delta;

    for(int i = 1; i < size; ++i) {
        out[i] = out[i-1] + Delta;
    }

    return &out[1];
}

CCTK_REAL * free_grid(CCTK_REAL * grid) {
    delete[] &grid[-1];
    return NULL;
}

int main(void) {
    hrscc::config::param::weno_alpha = config_weno_alpha;
    hrscc::config::param::weno_eps   = config_weno_epsilon;

    std::string base_name = "weno_derivative";
    std::ofstream file;

    test_function phi;
    diff_test_function diff_phi;

    CCTK_REAL * xp   = NULL;
    CCTK_REAL * fp   = NULL;
    CCTK_REAL * dfp  = NULL;
    CCTK_REAL * rdfp = NULL;

    cGH * cctkGH = new cGH;
    cctkGH->cctk_ash = new int[3];
    cctkGH->cctk_delta_space = new CCTK_REAL[3];
    cctkGH->cctk_levfac = new int[3];
    cctkGH->cctk_lsh = new int[3];

    cctkGH->cctk_levfac[0] = 1;
    cctkGH->cctk_levfac[1] = 1;
    cctkGH->cctk_levfac[2] = 1;

    int grid_n[config_ngrids];
    CCTK_REAL error[config_ngrids];

    hrscc::FDReconstruction<my_weno_reconstruction > fdrecon;

    int n = config_n0;
    for(int g = 0; g < config_ngrids; ++g) {
        grid_n[g] = n;

        cctkGH->cctk_ash[0] = 1;
        cctkGH->cctk_ash[1] = 1;
        cctkGH->cctk_ash[2] = n;
        cctkGH->cctk_lsh[0] = 1;
        cctkGH->cctk_lsh[1] = 1;
        cctkGH->cctk_lsh[2] = n;

        cctkGH->cctk_delta_space[0] = (config_x1 - config_x0) / n;
        cctkGH->cctk_delta_space[1] = (config_x1 - config_x0) / n;
        cctkGH->cctk_delta_space[2] = (config_x1 - config_x0) / n;

        xp   = new CCTK_REAL[n];
        fp   = new CCTK_REAL[n];
        dfp  = new CCTK_REAL[n];
        rdfp = new CCTK_REAL[n];

        for(int i = 0; i < n; ++i) {
            xp[i]  = config_x0 + i * cctkGH->cctk_delta_space[0];
            fp[i]  = phi(xp[i]);
            dfp[i] = diff_phi(xp[i]);
        }
        fdrecon.diff<hrscc::policy::z, hrscc::policy::minus>(cctkGH, fp, rdfp);

        error[g] = 0;
        for(int i = nghostzones; i < n - nghostzones; ++i) {
            error[g] = std::max(error[g], std::abs(dfp[i] - rdfp[i]));
        }

#ifdef config_output_full_data
        std::ostringstream sn;
        sn << n;
        std::string mystring(sn.str());
        file.open((base_name + "." + mystring + ".dat").c_str());

        file << "# 1:x 2:f(x) 3:df(x) 4:rdf(x)" << std::endl;
        for(int i = nghostzones; i < n - nghostzones; ++i) {
            file << xp[i] << " " << fp[i] << " " << dfp[i] << " " << rdfp[i];
            file << std::endl;
        }
        file.close();
#endif

        delete[] xp;
        delete[] fp;
        delete[] dfp;
        delete[] rdfp;

        n = static_cast<int>(static_cast<CCTK_REAL>(n) * config_grid_factor);
    }

    CCTK_REAL order[config_ngrids];
    CCTK_REAL avg_order;
    order[0] = 0;
    avg_order = 0;
    for(int g = 1; g < config_ngrids; ++g) {
        order[g] = - std::log(error[g]/error[g-1]) /
            std::log(static_cast<CCTK_REAL>(grid_n[g]) /
                    static_cast<CCTK_REAL>(grid_n[g-1]));
        avg_order += order[g] / (config_ngrids-1);
    }

#ifdef config_verbose
    std::cout << "n points\t\tLinfty error\t\torder\n";
#endif

#ifdef config_output_error
    file.open("weno_derivative.error.infty.dat");
#endif
    for(int g = 0; g < config_ngrids; ++g) {
#ifdef config_output_error
        file << grid_n[g] << " " << error[g] << " " << order[g] << std::endl;
#endif
#ifdef config_verbose
        std::cout << grid_n[g] << "\t\t\t" << error[g] << "\t\t" << order[g]
                  << std::endl;
#endif
    }
#ifdef config_output_error
    file.close();
#endif

#ifdef config_check_order
    assert(std::abs(avg_order - config_avg_order_reference) /
           config_avg_order_reference < config_toll);
#endif

#ifdef config_verbose
    std::cout << "\nAverage order = " << avg_order << std::endl;
#endif

    delete[] cctkGH->cctk_lsh;
    delete[] cctkGH->cctk_levfac;
    delete[] cctkGH->cctk_delta_space;
    delete[] cctkGH->cctk_ash;
    delete cctkGH;
};
