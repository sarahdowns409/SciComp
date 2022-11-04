#include <iostream>
#include "Grid2d.h"
#include "SL_method.h"
#include "Godunov.h"
#include <vector>
#include <cmath>

int main() {
    //------------------Semi-Lagrangian (SL)-------------------//
    SL_method sl;

    // Create grid for problem
    int N = 51;
    int M = N;
    double xmin = -1.;
    double xmax = 1.;
    double ymin = -1.;
    double ymax = 1.;
    Grid2d new_grid = Grid2d(N, M, xmin, xmax, ymin, ymax );
    sl.set_grid(new_grid);

    // Find dt proportional to dx
    double dt =  (new_grid.get_dx())/(10);

    //double dt =  0.1;

    std::cout << "dx: " << new_grid.get_dx() << std::endl;
    std::cout << "dt: " << dt << std::endl;
    std::cout << "dx/dt: " << new_grid.get_dx()/dt << std::endl;

    // allocate space in memory for vals
    std::vector<double> vals;
    vals.assign(N*M, 0.);

    // Initialize the solution
    for (int i = 0; i < N*M; i++){
        vals[i] = sl.find_phi(new_grid.x_from_n(i), new_grid.y_from_n(i));
    }
    std::vector<double> vals_i = vals;


    // Save VTK file for true solution
    //new_grid.print_VTK_format("../trueSol3.vtk");
    //new_grid.print_VTK_format(vals_i, "phi", "../trueSol3.vtk");

    // Compute solution at time 2pi
    std::vector<double> solution;
    solution.assign(N*M, 0.);
    sl.sl_method(vals, 0, 2*M_PI, dt);
    solution = sl.get_sol();

    // Output max error found
    double max_err = sl.error_eval(vals_i, solution);
    std::cout << "Max Error from SL: " << max_err << std::endl;

    // Save VTK file from SL Solutioin
    new_grid.print_VTK_format("../sl_q4.vtk");
    new_grid.print_VTK_format(solution, "phi", "../sl_q4.vtk");


    //------------------GODUNOV-------------------//

    // Initialize godunov method

    Godunov godunov;
    double eps = 0.001 * new_grid.get_dx();
    godunov.set_iter(5000);
    godunov.set_sl(sl);

    godunov.set_grid(new_grid);

    // perturb initial problem
    std::vector<double> new_vals=solution;
    godunov.perturb_vals(new_vals);
    max_err=sl.error_eval(vals_i,new_vals);
    std::cout<<"Perturbation error: "<<max_err <<std::endl;

    godunov.godunov_scheme(new_vals, vals_i, eps);
    new_grid.print_VTK_format("../gdnv_5000.vtk");
    new_grid.print_VTK_format(new_vals, "phi", "../gdnv_5000.vtk");

    max_err = sl.error_eval(vals_i, new_vals);
    std::cout << "Max Error from Godunov: " << max_err << std::endl;

    return 0;

}