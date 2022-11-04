#include "Godunov.h"


void Godunov::godunov_scheme(std::vector<double> &phi, std::vector<double> &phi_0, double dt) {
    std::vector<double> S_phi;
    // Creation of S
    S_phi.assign(grid.get_N()*grid.get_M(), 0.);
    for (int n = 0; n < grid.get_N() * grid.get_M(); n++){
        S_phi[n] = phi_0[n] / sqrt(phi_0[n]*phi_0[n] + dt*dt);
    }

    // Initialize current step for loop
    std::vector<double> curr = phi_0;

    int N = grid.get_N(), M = grid.get_M();
    double dx = grid.get_dx(), dy = grid.get_dy();
    double a = 0., b = 0., c = 0., d = 0.;

    // Creation of G
    double G_phi;
    for (int n = 0; n < iter; n++){
        // Set a, b, c, d
        for (int i = 0; i < N; i++){
            for (int j = 0; j < M; j++) {
                // Set derivatives for x
                if (i == 0){
                    b = (curr[grid.n_from_ij(i + 1, j)] - curr[grid.n_from_ij(i, j)]) / dx;
                }
                else if (i == N - 1){
                    a = (curr[grid.n_from_ij(i, j)] - curr[grid.n_from_ij(i - 1, j)]) / dx;
                }
                else{
                    a = (curr[grid.n_from_ij(i, j)] - curr[grid.n_from_ij(i - 1, j)]) / dx;
                    b = (curr[grid.n_from_ij(i + 1, j)] - curr[grid.n_from_ij(i, j)]) / dx;
                }
                // Set derivatives for y
                if (j == 0) {
                    d = (curr[grid.n_from_ij(i, j + 1 )] - curr[grid.n_from_ij(i, j)]) / dy;
                }
                else if (j == M - 1) {
                    c = (curr[grid.n_from_ij(i, j)] - curr[grid.n_from_ij(i, j - 1)]) / dy;
                }
                else{
                    c = (curr[grid.n_from_ij(i, j)] - curr[grid.n_from_ij(i, j - 1)]) / dy;
                    d = (curr[grid.n_from_ij(i, j + 1 )] - curr[grid.n_from_ij(i, j)]) / dy;
                }

                // Find if  a, b, c, d are positive or negative - separate to pos or neg values
                double ap, am, bp, bm, cp, cm, dp, dm;
                //std::cout << "a: " << a << std::endl;
                ap = a >= 0. ? a : 0.;
                am = a <  0. ? a : 0.;
                //std::cout << "ap: " << ap << std::endl;
                //std::cout << "am: " << am << std::endl;
                bp = b >= 0. ? b : 0.;
                bm = b <  0. ? b : 0.;
                cp = c >= 0. ? c : 0.;
                cm = c <  0. ? c : 0.;
                dp = d >= 0. ? d : 0.;
                dm = d <  0. ? d : 0.;

                if (phi_0[grid.n_from_ij(i,j)] > 0){
                    G_phi = sqrt(std::max(ap*ap, bm*bm) + std::max(cp*cp, dm*dm)) - 1;
                }
                else if(phi_0[grid.n_from_ij(i,j)] < 0){
                    G_phi = sqrt(std::max(am*am, bp*bp) + std::max(cm*cm, dp*dp)) - 1;
                }
                else {
                    G_phi = 0;
                }
                // From equation (29) of PDF
                phi[grid.n_from_ij(i, j)] = curr[grid.n_from_ij(i, j)] - dt * S_phi[grid.n_from_ij(i, j)]  * G_phi;
                // End of innermost loop
            }
        }
        // Update for next iteration of outermost loop
        curr = phi;
    }
}

void Godunov::perturb_vals(std::vector<double> &phi) {
    for (int n = 0; n < grid.get_N() * grid.get_M(); n++ ){
        phi[n] = phi[n] + (sin (n))/10;
    }
}