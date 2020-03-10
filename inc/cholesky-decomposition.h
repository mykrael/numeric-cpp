#pragma once

#include <cmath>
#include "matrix.h"
#include "math-functions.cpp"
#include "triangular-matrix-substitution.h"

/*
 * Contains all functions and classes using the cholesky decomposition
 */
namespace CHOLESKY{

    /*
     * Stores the matrix L of a Cholesky decomposition
     */
    template<int N, class T>
        struct CholeskyDecomposition{
            Matrix<N, N, T> L;
            CholeskyDecomposition(Matrix<N,N,T> A){
                if(!calculateCholesky(A)){
                    std::cout << "Matrix is not positive definite" << std::endl;
                }
            }
            bool calculateCholesky(Matrix<N, N, T> A){
                Matrix<N, N, T> L;
                for(int i = 0; i < N; i++){
                    /*
                     * Calculate L[i][i]
                     */
                    T sum = T(0);
                    for(int k = 0; k < i; k++){
                        sum += L[i][k] * L[i][k];
                    }
                    sum = A[i][i] - sum;
                    if (sum >= 0){
                        L[i][i] = std::sqrt(sum);
                    }else{
                        return false;
                    }
                    /*
                     * Calculate all entries of L in i-th row)
                     */
                    for(int j = i+1; j < N; j++){
                        sum = T(0);
                        for(int k = 0; k < i; k++){
                            sum += L[j][k] * L[i][k];
                        } 
                        L[j][i] = (A[i][j] - sum)/L[i][i];
                    }
                }
                this->L = L;
                return true;    
            }
            Matrix<N,1,T> solve(Matrix<N,1,T> b){
                Matrix<N,1,T> z = forSub(L, b);
                return backSub(~L, z);

            }

        };



}
