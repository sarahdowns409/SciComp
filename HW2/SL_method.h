#ifndef HW2_SL_METHOD_H
#define HW2_SL_METHOD_H

#include "Grid2d.h"
#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>

// Semi-Langrangian Method
class SL_method {
private:
    Grid2d grid;
    std::vector<double> sol;
    std::vector<double> vel_x;
    std::vector<double> vel_y;

    // Compute central finite difference for phi_xx or phi_yy
    double cfd_xx(int i, int j, std::vector <double> &func);
    double cfd_yy(int i, int j, std::vector <double> &func);
    // Compute forward finite difference for phi_xx or phi_yy;
    double ffd_xx(int i, int j, std::vector <double> &func);
    double ffd_yy(int i, int j, std::vector <double> &func);
    // Compute backward finite difference for phi_xx or phi_yy;
    double bfd_xx(int i, int j, std::vector <double> &func);
    double bfd_yy(int i, int j, std::vector <double> &func);
    double finite_difference_xx(int i, int j, std::vector <double> &func);
    double finite_difference_yy(int i, int j, std::vector <double> &func);
    double minmod_xx(int i, int j, std::vector <double> &func);
    double minmod_yy(int i, int j, std::vector <double> &func);

    // Find departure point at given (x, y)
    void find_trajectory(int n, double &x_d, double &y_d, double dt);
    // Return result of quadratic interpolation at (x, y) for func phi
    double quad_interpolation(std::vector<double> &func,double x, double y);

public:
    // Set initial velocity
    void set_velocity(std::vector<double> &vel_x0, std::vector<double> &vel_y0);
    // Set grid
    void set_grid(Grid2d & new_grid){grid = new_grid;}
    // Access solution
    std::vector<double> get_sol(){ return sol; }
    // Compute phi at (x, y)
    double find_phi(double x, double y);
    // Compute solution for semi-lagrangian method over [t0, tf] in T steps
    void sl_method(std::vector<double> &sol_0, double t0, double tf, double dt);
    // Compute error norm for solution versus true solution
    double error_eval(std::vector<double> &tru_sol, std::vector<double> &sol);
    void store_errors(std::vector<double> &sol, double t0, double tf);
};


#endif //HW2_SL_METHOD_H