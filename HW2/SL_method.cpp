#include "SL_method.h"

// Returns the departure point
// set initial condition phi(x,y,0)
double SL_method::find_phi(double x, double y) {
    // question 1 and 2
    // return sqrt(pow(x - 0.25, 2) + pow(y,  2)) - 0.2;

    // question 3
    //return pow(x-0.25,2)+pow(y,2)-pow(0.2,2);

    //question 4
    return pow(x,2)+pow(5*y/4- sqrt(abs(x)),2)-0.1;
}

//define the velocity field V
void SL_method::set_velocity(std::vector<double> &vel_x0, std::vector<double> &vel_y0){
    for (int c = 0; c < grid.get_N() * grid.get_M(); c++) {
        vel_x0[c] = -sin(grid.y_from_n(c));
        vel_y0[c] =  cos(grid.x_from_n(c));
    }
}

// finds trajectory
void SL_method::find_trajectory(int n, double &x_d, double &y_d, double dt) {
    // RK2
    double x_0 = grid.x_from_n(n), y_0 = grid.y_from_n(n);
    double x_s = x_0 - 0.5 * dt * vel_x[n];
    double y_s = y_0 - 0.5 * dt * vel_y[n];

    double vx_d = -y_s, vy_d =  x_s;
    x_d = x_0 -       dt * vx_d;
    y_d = y_0 -       dt * vy_d;
}

// central finite difference derivative xx
double SL_method::cfd_xx(int i, int j, std::vector<double> & func) {
    double cfd;
    int n   = grid.n_from_ij(  i    , j);
    int np1 = grid.n_from_ij(i + 1, j);
    int nm1 = grid.n_from_ij(i - 1, j);

    cfd = (func[np1] - 2 * func[n] + func[nm1]) / (pow(grid.get_dx(), 2));
    return cfd;
}

// central finite difference derivative yy
double SL_method::cfd_yy(int i, int j, std::vector<double> &func) {
    double cfd;
    int n   = grid.n_from_ij(i,   j);
    int nm1 = grid.n_from_ij(i, j - 1);
    int np1 = grid.n_from_ij(i, j + 1);

    cfd = (func[np1] -2 * func[n] + func[nm1]) / (pow(grid.get_dy(), 2));
    return cfd;
}

// forward finite difference derivative xx (for edges)
double SL_method::ffd_xx(int i, int j, std::vector<double> &func) {
    double ffd;
    int n   = grid.n_from_ij(  i     , j);
    int np1 = grid.n_from_ij(i + 1, j);
    int np2 = grid.n_from_ij(i + 2, j);
    int np3 = grid.n_from_ij(i + 3, j);

    ffd = (2 * func[n] - 5 * func[np1] + 4 * func[np2] - func[np3]) / pow(grid.get_dx(), 3);
    return ffd;
}

// forward finite difference derivative yy (for edges)
double SL_method::ffd_yy(int i, int j, std::vector<double> &func) {
    double ffd;
    int n   = grid.n_from_ij(i,   j);
    int np1 = grid.n_from_ij(i, j + 1);
    int np2 = grid.n_from_ij(i, j + 2);
    int np3 = grid.n_from_ij(i, j + 3);

    ffd = (2 * func[n] - 5 * func[np1] + 4 * func[np2] - func[np3]) / pow(grid.get_dy(), 3);
    return ffd;
}

// backward finite difference derivative xx (for edges)
double SL_method::bfd_xx(int i, int j, std::vector<double> &func) {
    double bfd;
    int n   = grid.n_from_ij(  i     , j);
    int nm1 = grid.n_from_ij(i - 1, j);
    int nm2 = grid.n_from_ij(i - 2, j);
    int nm3 = grid.n_from_ij(i - 3, j);

    bfd = (2 * func[n] - 5 * func[nm1] + 4 * func[nm2] - func[nm3]) / pow(grid.get_dx(), 3);
    return bfd;
}

// backward finite difference derivative yy (for edges)
double SL_method::bfd_yy(int i, int j, std::vector<double> &func) {
    double bfd;
    int n   = grid.n_from_ij(i,    j);
    int nm1 = grid.n_from_ij(i, j - 1);
    int nm2 = grid.n_from_ij(i, j - 2);
    int nm3 = grid.n_from_ij(i, j - 3);

    bfd = (2 * func[n] - 5 * func[nm1] + 4 * func[nm2] - func[nm3]) / pow(grid.get_dy(), 3);
    return bfd;
}

//compute finite difference derivative for xx
double SL_method::finite_difference_xx(int i, int j, std::vector<double> &func) {
    int N = grid.get_N();

    if ( i == 0 ) return ffd_xx(i, j, func);
    else if ( i == (N - 1) ) return bfd_xx(i, j, func);
    else return cfd_xx(i, j, func);

}

//compute finite difference derivative for yy
double SL_method::finite_difference_yy(int i, int j, std::vector<double> &func) {
    int M = grid.get_M();

    if ( j == 0 ) return ffd_yy(i, j, func);
    else if ( j == (M - 1) ) return bfd_yy(i, j, func);
    else return cfd_yy(i, j, func);
}

//find minmod of xx derivative
double SL_method::minmod_xx(int i, int j, std::vector<double> &func){
    double a = finite_difference_xx(i, j, func);
    double b = finite_difference_xx(i + 1, j, func);
    double c = finite_difference_xx(i, j + 1, func);
    double d = finite_difference_xx(i + 1, j + 1, func);

    int sign_a = std::signbit(a), sign_b = std::signbit(b), sign_c = std::signbit(c), sign_d = std::signbit(d);
    int signs = sign_a + sign_b + sign_c + sign_d;
    std::vector <double> minmods;
    minmods.assign(4, 0.);
    minmods[0] = a;
    minmods[1] = b;
    minmods[2] = c;
    minmods[3] = d;

    if (signs == 0 || signs == 4){ // All positive or all negative
        double minmod = minmods[0];
        for ( int k = 1; k < 4; k++ ){
            minmod = minmod*minmod <= minmod*minmods[k] ? minmod : minmods[k];
        }
        return minmod;
    }
    else return 0;
}

//find minmod of yy derivative
double SL_method::minmod_yy(int i, int j, std::vector<double> & func){
    double a = finite_difference_yy(i, j, func);
    double b = finite_difference_yy(i + 1, j, func);
    double c = finite_difference_yy(i, j + 1, func);
    double d = finite_difference_yy(i + 1, j + 1, func);

    int sign_a = std::signbit(a), sign_b = std::signbit(b), sign_c = std::signbit(c), sign_d = std::signbit(d);
    int signs = sign_a + sign_b + sign_c + sign_d;
    std::vector <double> minmods;
    minmods.assign(4, 0.);
    minmods[0] = a;
    minmods[1] = b;
    minmods[2] = c;
    minmods[3] = d;

    if (signs == 0 || signs == 4){ // All positive or all negative
        double minmod = minmods[0];
        for ( int k = 1; k < 4; k++ ){
            minmod = minmod*minmod <= minmod*minmods[k] ? minmod : minmods[k];
        }
        return minmod;
    }
    else return 0;
}

// Returns the interpolation of the solution at (x,y) (departure point)
double SL_method::quad_interpolation(std::vector<double> &func, double x, double y) {
    //double x_c = x;
    //double y_c = y;

    double dx = grid.get_dx(), dy = grid.get_dy();
    int i, j;

    // If (x,y) are outside domain, set to grid point closest to (x,y) and get i and j
    if (x <= grid.get_xmin() || x >= grid.get_xmax() || y <= grid.get_ymin() || y >= grid.get_ymax() ) {
        // X in grid, but  Y is not:
        if (x > grid.get_xmin() && x < grid.get_xmax()) {
            i = floor( (x - grid.get_xmin()) / dx);
            if (y <= grid.get_ymin()) j = 0;
            else if (y >= grid.get_ymax()) j = grid.get_M() - 2;
        }
        // Y in grid, X is not:
        if (y > grid.get_ymin() && y < grid.get_ymax()) {
            j = floor((y - grid.get_ymin()) / dy);
            if (x <= grid.get_xmin()) i = 0;
            else if (x >= grid.get_xmax() ) i = grid.get_N() - 2;
        }
        // X and Y are both out of grid:
        else{
            if (x <= grid.get_xmin()) i = 0;
            else if (x >= grid.get_xmax()) i = grid.get_N() - 2;
            if (y <= grid.get_ymin()) j = 0;
            else if (y >= grid.get_ymax()) grid.get_M() - 2;
        }
    }
    else {
        i = floor( (x - grid.get_xmin()) / dx);
        j = floor( (y - grid.get_ymin()) / dy);
    }

    double phi;

    double x_i = grid.get_xmin() + i * dx;
    double y_j = grid.get_ymin() + j * dy;
    double x_ip1 = x_i + dx;
    double y_jp1 = y_j + dy;

    // Use quadratic interpolation to get value at x
    // (i.e. think weighted avg)
    // (i, j), (i + 1, j), (i, j + 1), (i + 1, j + 1) are the corners of the cell C
    phi  = func[grid.n_from_ij(i    ,   j  )]  * ( x_ip1 - x   ) * ( y_jp1 - y   ) / (dx*dy) ;
    phi += func[grid.n_from_ij(i+1,   j  )]  * ( x     - x_i ) * ( y_jp1 - y   ) / (dx*dy) ;
    phi += func[grid.n_from_ij(i    , j+1)]  * ( x_ip1 - x   ) * ( y     - y_j ) / (dx*dy) ;
    phi += func[grid.n_from_ij(i+1, j+1)]  * ( x     - x_i ) * ( y     - y_j ) / (dx*dy) ;


    // Finding the minmod of phixx and phiyy
    double minmod_phixx = minmod_xx(i, j, func);
    double minmod_phiyy = minmod_yy(i, j, func);

    phi -= (minmod_phixx * ( x - x_i ) * ( x_ip1 - x )) / 2;
    phi -= (minmod_phiyy * ( y - y_j ) * ( y_jp1 - y )) / 2;

    return phi;
}

// runs SL Method
void SL_method::sl_method(std::vector<double> &sol_0, double t0, double tf, double dt) {
    vel_x.assign(grid.get_N()*grid.get_M(), 0.);
    vel_y.assign(grid.get_N()*grid.get_M(), 0.);
    set_velocity(vel_x, vel_y);

    double x_d = 0, y_d = 0;
    int T = ceil((tf - t0) / (dt));

    // Initialize solution and solve at each timestep
    sol.assign(grid.get_N() * grid.get_M(), 0.);
    for (int t = 0; t < T; t++) {
        for (int n = 0; n < grid.get_N() * grid.get_M(); n++){
            find_trajectory(n, x_d, y_d, dt);
            sol[n] = quad_interpolation(sol_0, x_d, y_d);
        }
        // Update previous solution step
        sol_0 = sol;
    }
}

// calulates max error of found solution vs true
// don't call if true solution is not known!
double SL_method::error_eval(std::vector<double> &tru_sol, std::vector<double> &sol) {
    long length = tru_sol.size();
    double max_error = 0.;
    double error;

    for (int n = 0; n < length; n++) {
        error = sqrt(pow(tru_sol[n] - sol[n], 2));
        if (error > max_error) {
            max_error = error;
        }
        //std::cout << "n"<< n<< " sol[n]" << sol[n] << std::endl;
    }
    return max_error;
}