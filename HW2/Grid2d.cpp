#include "Grid2d.h"

Grid2d::Grid2d() {

}

Grid2d::Grid2d(long NN, long MM, double xlo, double xhi, double ylo, double yhi) {
    N = NN;
    M = MM;
    xmin = xlo;
    xmax = xhi;
    ymin = ylo;
    ymax = yhi;

    dx = (xmax - xmin) / (double(N - 1));
    dy = (ymax - ymin) / (double(M - 1));
}

double Grid2d::get_dx() {
    return dx;
}

double Grid2d::get_dy() {
    return dy;
}

double Grid2d::get_xmin() {
    return xmin;
}

double Grid2d::get_xmax() {
    return xmax;
}

double Grid2d::get_ymin() {
    return ymin;
}

double Grid2d::get_ymax() {
    return ymax;
}

long Grid2d::n_from_ij(long i, long j) {
    return N * j + i;
}

long Grid2d::i_from_n(long n) {
    return n % N;
}

long Grid2d::j_from_n(long n) {
    return n / N;
}

double Grid2d::x_from_n(long n) {
    return xmin + double(i_from_n(n)) * dx;
}

double Grid2d::y_from_n(long n) {
    return ymin + double(j_from_n(n)) * dy;
}

void Grid2d::print_VTK_format(std::string output_file) {
    long num_of_leaf_cells;
    long node_of_cell[4];
    double x, y;
    FILE *outFile = fopen(output_file.c_str(),"w");
    fprintf(outFile,"# vtk DataFile Version 2.0 \n");
    fprintf(outFile,"Quadtree Mesh \n");
    fprintf(outFile,"ASCII \n");
    fprintf(outFile,"DATASET UNSTRUCTURED_GRID \n");

    //first output the list of nodes
    fprintf(outFile,"POINTS %d double \n",N*M);
    for (long n=0; n<M*N; n++){
        x = x_from_n(n);
        y = y_from_n(n);
        fprintf(outFile,"%e %e %e\n",x,y,0.0);
    }
    // then output the list of cells. each cell is composed of four nodes
    // num_of_leaves is the total number of cells of the mesh
    num_of_leaf_cells = ( N - 1 ) * ( M - 1 );
    fprintf(outFile,"CELLS %d %d \n",num_of_leaf_cells,5*num_of_leaf_cells);
    for (int i=0; i<N-1; i++)
        for (int j=0; j<M-1; j++)
        {
            fprintf(outFile,"%d %d %d %d %d\n",4,n_from_ij(i,j), n_from_ij(i+1,j), n_from_ij(i+1,j+1), n_from_ij(i,j+1));
        }

    fprintf(outFile,"CELL_TYPES %d \n",num_of_leaf_cells);
    for (long n=0; n<num_of_leaf_cells; n++)    fprintf(outFile,"%d \n",9);
    fprintf(outFile,"POINT_DATA %d \n", int(N * M));
    fclose (outFile);
}

//% this function write the values of the vector F into the vtk file. before
//using it, the .vtk file must have been initialized with all the grid infos
void Grid2d::print_VTK_format( std::vector<double> &F, std::string data_name,
                               std::string file_name )
{
    long num_of_nodes;
    num_of_nodes = N * M;
    FILE *outFile;
    outFile = fopen(file_name.c_str(),"a");
    fprintf(outFile,"SCALARS %s double 1 \n",data_name.c_str());
    fprintf(outFile,"LOOKUP_TABLE default \n");
    for (long n=0; n<num_of_nodes; n++) fprintf(outFile,"%e \n",F[n]);
    fclose (outFile);
}
