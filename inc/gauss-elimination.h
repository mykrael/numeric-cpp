#include "matrix.h"

namespace GAUSS{


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

    template<int N, class T>
        Matrix<N, 1, T> solve(Matrix<N, N, T> A, Matrix<N, 1, T> b){
            Matrix<N, 1, T> x;
            // Gauss elimination
            for(int j = 0; j < N-1; j++){
                for(int i = j+1; i < N; i++){
                    T l = A[i][j]/A[j][j];
                    b[i][0] = b[i][0] - l * b[j][0];
                    for(int k = j; k < N; k++){
                        A[i][k] = A[i][k] - l * A[j][k];
                    }
                }
            }
            x = backSub(A, b);

            return x;

        }
}

