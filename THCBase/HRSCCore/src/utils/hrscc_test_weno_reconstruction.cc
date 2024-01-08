#include <cassert>
#include <cmath>
#include <iostream>
#include <fstream>
#include <ostream>
#include <string>

#include <cctk.h>

#include <utils.hh>

#include <hrscc_config_par.hh>
#include <hrscc_eno_stencil.hh>
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
#define config_x0                   0
#define config_x1                   1
#define config_n0                   20
#define config_ngrids               10
#define config_grid_factor          1.2
#define config_nhres                1000

// Number of points used to ``integrate'' the function
#define config_average_npoints      10000
// Test function
#define config_test_function        x*utils::pow<2>(std::sin(4*M_PI*x))

// Define to output the data onto files
//#define config_output_full_data
// Define to save the Linf errors
//#define config_output_error
// Verbose
#define config_verbose


// Reference value for the convergence order
#define config_check_order
#define config_avg_order_reference  4.84094
#define config_toll                 1e-3
///////////////////////////////////////////////////////////////////////////////

#define my_eno_stencil              hrscc::ENOStencil<config_eno_width>
#define my_weno_weights             hrscc::WENOWeights<config_eno_width>
#define my_weno_stencil             hrscc::WENOStencil<config_eno_width,          \
                                        config_weno_type, config_weno_optim>
#define my_weno_limiter             hrscc::WENOLimiter<config_eno_width,          \
                                        my_weno_stencil::width,                   \
                                        config_weno_limit>
#define my_weno_reconstruction      hrscc::WENOReconstruction<my_eno_stencil,     \
                                        my_weno_limiter, my_weno_stencil,         \
                                        my_weno_weights >

#define nghostzones                 my_weno_reconstruction::width
///////////////////////////////////////////////////////////////////////////////

typedef CCTK_REAL (* function)(CCTK_REAL);

template<class F>
CCTK_REAL integrate(F const & f, int N, CCTK_REAL a = -1, CCTK_REAL b = 1) {
    CCTK_REAL xi;
    CCTK_REAL xp;
    CCTK_REAL value = 0;

    for(int i = 1; i < N+1; ++i) {
        xi = static_cast<double>(2*i - 1) / static_cast<double>(2*N) * M_PI;
        xi = std::cos(xi);
        xp = 0.5*(a*(1 - xi) + b*(xi + 1));
        value += M_PI/N * 0.5*(b - a) *  std::sqrt(1 - xi*xi) * f(xp);
    }

    return value;
}

class test_function {
    public:
        CCTK_REAL operator()(CCTK_REAL x) const {
            return config_test_function;
        }
};

CCTK_REAL * make_grid(CCTK_REAL a, CCTK_REAL b, int npoints,
        bool interfaces = false) {
    CCTK_REAL Delta = (b - a) / (npoints - 1);
    int size = npoints + 2 - (2*nghostzones - 1)*interfaces;
    CCTK_REAL * out = new CCTK_REAL[size];

    out[0] = a - Delta + (nghostzones - 0.5)*interfaces*Delta;

    for(int i = 1; i < size; ++i) {
        out[i] = out[0] + Delta*i;
    }

    return &out[1];
}

CCTK_REAL * free_grid(CCTK_REAL * grid) {
    delete[] &grid[-1];
    return NULL;
}

int main(void) {
    int i, j, g;
    CCTK_REAL a, b;

    hrscc::config::param::weno_alpha = config_weno_alpha;
    hrscc::config::param::weno_eps   = config_weno_epsilon;

    test_function phi;

    CCTK_REAL dx_base;
    CCTK_REAL * x_base  = NULL;
    CCTK_REAL * f_base  = NULL;

    CCTK_REAL * x_vert  = NULL;
    CCTK_REAL * f_vert  = NULL;
    CCTK_REAL * fm_vert = NULL;
    CCTK_REAL * fp_vert = NULL;

    int grid_n[config_ngrids];
    CCTK_REAL error[config_ngrids];

    std::string base_name = "weno_reconstruction";
    std::ofstream file;

    my_weno_reconstruction weno;//config_weno_alpha, config_weno_epsilon);

#ifdef config_output_full_data
    CCTK_REAL * x_maxr  = NULL;
    CCTK_REAL * f_maxr  = NULL;

    x_maxr = make_grid(config_x0, config_x1, config_nhres);
    f_maxr = new CCTK_REAL[config_nhres];
    for(i = 0; i < config_nhres; ++i) {
        f_maxr[i] = phi(x_maxr[i]);
    }

    file.open((base_name + ".ref.dat").c_str());
    for(i = 0; i < config_nhres; ++i) {
        file << x_maxr[i] << " " << f_maxr[i] << std::endl;
    }
    file.close();

    x_maxr = free_grid(x_maxr);
    delete[] f_maxr;
#endif

    int n = config_n0;
    int size;
    for(g = 0; g < config_ngrids; ++g) {
        grid_n[g] = n;
        size = n - 2*nghostzones + 3;

        x_base  = make_grid(config_x0, config_x1, n);
        dx_base = x_base[1] - x_base[0];
        f_base  = new CCTK_REAL[n];

        for(i = 0; i < n; ++i) {
            a = 0.5*(x_base[i-1] + x_base[i]);
            b = 0.5*(x_base[i] + x_base[i+1]);
            f_base[i] = integrate(phi, config_average_npoints, a, b) / (b - a);
        }

        x_vert  = make_grid(config_x0, config_x1, n, true);
        f_vert  = new CCTK_REAL[size];
        fm_vert = new CCTK_REAL[size];
        fp_vert = new CCTK_REAL[size];
        for(i = 0; i < size; ++i) {
            f_vert[i] = phi(x_vert[i]);
        }

        error[g] = 0;
        for(i = nghostzones - 1; i < n - nghostzones; ++i) {
            j = i - nghostzones + 1;

            fm_vert[j] = weno.reconstruct<hrscc::policy::minus>(
                    dx_base, &f_base[i]);
            fp_vert[j] = weno.reconstruct<hrscc::policy::plus>(
                    dx_base, &f_base[i]);

            error[g] = std::max(error[g], std::abs(f_vert[j] -
                        0.5*(fm_vert[j] + fp_vert[j])));
        }

#ifdef config_output_full_data
        std::ostringstream sn;
        sn << n;
        std::string mystring(sn.str());
        file.open((base_name + "." + mystring + ".dat").c_str());

        file << "# index 0\n";
        file << "# 1: x 2: < u >\n";
        file << "# index 1\n";
        file << "# 1: x 2: u 3: um 4: up 5: 0.5*(um+up)\n";
        file << "# index 2\n";
        file << "# 1: x 2: um-u 3: up-u 4: 0.5*(um+up)-u\n";

        for(i = 0; i < n; ++i) {
            file << x_base[i] << " " << f_base[i] << std::endl;
        }
        file << std::endl << std::endl;
        for(i = 0; i < size; ++i) {
            file << x_vert[i]  << " " << f_vert[i]  << " ";
            file << fm_vert[i] << " " << fp_vert[i]
                 << 0.5*(fm_vert[i] + fp_vert[i]) << std::endl;
        }
        file << std::endl << std::endl;
        for(i = 0; i < size; ++i) {
            file << x_vert[i] << " " << fm_vert[i] - f_vert[i]
                 << " " << fp_vert[i] - f_vert[i]
                 << 0.5*(fp_vert[i] + fm_vert[i]) - f_vert[i] << std::endl;
        }
        file.close();
#endif

        x_base = free_grid(x_base);
        delete[] f_base;

        x_vert = free_grid(x_vert);
        delete[] f_vert;
        delete[] fm_vert;
        delete[] fp_vert;

        n = static_cast<int>(static_cast<CCTK_REAL>(n) * config_grid_factor);
    }

    CCTK_REAL order[config_ngrids];
    CCTK_REAL avg_order;
    order[0] = 0;
    avg_order = 0;
    for(g = 1; g < config_ngrids; ++g) {
        order[g] = - std::log(error[g]/error[g-1]) /
            std::log(static_cast<CCTK_REAL>(grid_n[g]) /
                    static_cast<CCTK_REAL>(grid_n[g-1]));
        avg_order += order[g] / (config_ngrids-1);
    }

#ifdef config_verbose
    std::cout << "n points\t\tLinfty error\t\torder\n";
#endif

#ifdef config_output_error
    file.open("weno_reconstruction.error.infty.dat");
#endif
    for(g = 0; g < config_ngrids; ++g) {
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

}
