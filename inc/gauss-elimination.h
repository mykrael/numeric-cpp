#pragma once
#include "math-functions.cpp"
#include "matrix.h"
#include "triangular-matrix-substitution.h"

namespace GAUSS{

    enum PIVOT_SEARCH{
        NO_SEARCH,
        COLLUMN_SEARCH
    };

    template<int N, class T>
        bool collumnPivotSearch(Matrix<N, N, T> A, Matrix<N, N, T>& P_j, int j){
            int search_index = j;
            int max_index = search_index;
            T max_entry = T(0);
            bool found_non_zero = max_entry == T(0);
            while(search_index < N){
                search_index++;
                T entry = A[j][search_index];
                if(TYPESPECIFIC::abs(entry > max_entry)){
                    max_index = search_index;
                    max_entry = entry;
                    found_non_zero = true;
                }
            } 
            for(int row = 0; row < N; row++){
                for(int col = 0; col < N; col++){
                    T entry = T(0);
                    if(row == j){
                        if(col == max_index){
                            entry = T(1);
                        }else{
                            entry = T(0);
                        }
                    }else if(row == max_index){
                        if(col == j){
                            entry = T(1);
                        }else{
                            entry = T(0);
                        }
                    }else{
                        if(col == row){
                            entry = T(1);
                        }else{
                            entry = T(0);
                        }
                    }
                    P_j[row][col] = entry;
                }
            }
            return found_non_zero;
        }

    template<int N, class T>
        bool gaussStep(Matrix<N, N, T>& A, Matrix<N, 1, T>& b, Matrix<N, N, T>& L_j, 
                Matrix<N,N,T>& L_j_inv, int j){
            Matrix<N, 1, T> x(b);
            Matrix<N, N, T> I;
            L_j = I;
            L_j_inv = I;
            for(int i = j+1; i < N; i++){
                T l = A[i][j]/A[j][j];
                L_j[i][j] = -l;
                L_j_inv[i][j] = l;
                b[i][0] = b[i][0] - l * b[j][0];
                for(int k = j; k < N; k++){
                    A[i][k] = A[i][k] - l * A[j][k];
                    if(j != k && A[i][k] == 0){
                        return false;
                    }
                }
            }
            return true;
        }

    template<int N, class T>
        struct LRCombination{
           Matrix<N,N,T> L; 
           Matrix<N,N,T> R;
           LRCombination(Matrix<N,N,T> L, Matrix<N,N,T> R) : L(L), R(R){};
           LRCombination(Matrix<N,N,T> A, Matrix<N,1,T> b){
                Matrix<N,N,T> L_j;
                Matrix<N,N,T> L_j_inv;
                Matrix<N,N,T> L;
                Matrix<N,N,T> R(A);
                for(int j = 0; j < N-1; j++){
                    if(!gaussStep(A, b, L_j, L_j_inv, j)){
                        std::cout << "Matrix is not regular or a pivot element is zero." << std::endl; 
                        LRCombination<N,T> error;
                    }
                    R = L_j * R;
                    L = L * L_j_inv;
                }
                this->R = R;
                this->L = L;
            }
           LRCombination(){};
           Matrix<N,1,T> solve(Matrix<N,1,T> b){
               Matrix<N,1,T> z;
               Matrix<N,1,T> x;
               z = forSub(L, b);
               x = backSub(R, z);
               return x;
           }
        };


    template<int N, class T>
        Matrix<N, 1, T> solve(Matrix<N, N, T> A, Matrix<N, 1, T> b, PIVOT_SEARCH search_algo = COLLUMN_SEARCH){
            Matrix<N, N, T> L_j;
            Matrix<N, N, T> L_j_inv;
            Matrix<N, N, T> P_j;
            Matrix<N, 1, T> old_b(b);
            Matrix<N, 1, T> x;
            // Gauss elimination
            for(int j = 0; j < N-1; j++){
                if(search_algo == COLLUMN_SEARCH){
                    if(!collumnPivotSearch(A, P_j, j)){
                        std::cout << "ERROR: The Matrix is not regular!" << std::endl;
                        return old_b;
                    }
                    A = P_j * A;
                    b = P_j * b;
                }
                gaussStep(A, b, L_j, L_j_inv, j);
            }
            x = backSub(A, b);
            return x;
        }


}

