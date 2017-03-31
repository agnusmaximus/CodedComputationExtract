#include <iostream>
#include <time.h>
#include <stdlib.h>
#include <cblas.h>
#include "../util.cpp"

using namespace std;

int main(void) {
    //Load data
    double *matrix = load_submatrix(0, N_COLS, N_ROWS);
    double *vec = load_input_vector(N_COLS);
    double *out = generate_empty_vector(N_ROWS);
    long long int t1 = get_time();
    matrix_vector_multiply(matrix, vec, out, N_COLS, N_ROWS);
    long long int t2 = get_time();
    cout << "ELAPSED TIME: " << t2-t1 << endl;
}
