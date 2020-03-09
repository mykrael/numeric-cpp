#include "matrix.h"

// Backwards substitution for upper right triangular matrix
template<int N, class T>
    Matrix<N,1,T> backSub(Matrix<N,N,T> R, Matrix<N,1,T> b){
        Matrix<N,1,T> x;
        for(int i = N-1; i >= 0; i--){
            T sum = T(0);
            for(int j = i+1; j < N; j++){
                sum += R[i][j] * x[j][0];
            }
            x[i][0] = (b[i][0] - sum)/R[i][i];
        }
        return x;
    }

// Forwards substitution for lower left triangular matrix
template<int N, class T>
    Matrix<N,1,T> forSub(Matrix<N,N,T> R, Matrix<N,1,T> b){
        Matrix<N,1,T> x;
        for(int i = 0; i < N; i++){
            T sum = T(0);
            for(int j = 0; j < i; j++){
                sum += R[i][j] * x[j][0];
            }
            x[i][0] = (b[i][0] - sum)/R[i][i];
        }
        return x;
    }

