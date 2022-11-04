#ifndef HW2_GODUNOV_H
#define HW2_GODUNOV_H

#include "Grid2d.h"
#include "SL_method.h"
// Godunov Scheme
class Godunov{
private:
    SL_method sl;
    Grid2d grid;
    int iter;

public:
    // Set sl_method
    void set_sl(SL_method & new_sl){ sl = new_sl; }
    void set_iter(int num_iters){ iter = num_iters; }
    void set_grid(Grid2d & new_grid){ grid = new_grid; }
    void godunov_scheme(std::vector<double> &phi, std::vector<double> &phi_0, double dt);
    void perturb_vals(std::vector<double> &phi);
};

#endif //HW2_GODUNOV_H
