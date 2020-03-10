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
            int search_index = j+1;
            int max_index = search_index;
            T max_entry = T(0);
            bool found_non_zero = !(A[j][j] == T(0));
            while(search_index < N){
                T entry = A[search_index][j];
                if(TYPESPECIFIC::abs(entry) > max_entry){
                    max_index = search_index;
                    max_entry = entry;
                    found_non_zero = true;
                }
                search_index++;
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
                Matrix<N,N,T>& L_j_inv, Matrix<N, N, T>& P_j, int j, PIVOT_SEARCH search_algo = COLLUMN_SEARCH){
            Matrix<N, 1, T> x(b);
            Matrix<N, N, T> I;
            L_j = I;
            L_j_inv = I;
            if(search_algo == COLLUMN_SEARCH){
                if(!collumnPivotSearch(A, P_j, j)){
                    std::cout << "ERROR: The Matrix is not regular in step " << j << std::endl;
                    return false;
                }
                A = P_j * A;
                b = P_j * b;
            }
            for(int i = j+1; i < N; i++){
                T l = A[i][j]/A[j][j];
                L_j[i][j] = -l;
                L_j_inv[i][j] = l;
                b[i][0] = b[i][0] - l * b[j][0];
                for(int k = j; k < N; k++){
                    A[i][k] = A[i][k] - l * A[j][k];
                }
            }
            return true;
        }

    template<int N, class T>
        struct LRCombination{
           Matrix<N,N,T> L; 
           Matrix<N,N,T> R;
           Matrix<N,N,T> P;
           LRCombination(Matrix<N,N,T> L, Matrix<N,N,T> R) : L(L), R(R){};
           LRCombination(Matrix<N,N,T> A, Matrix<N,1,T> b, PIVOT_SEARCH search_algo = COLLUMN_SEARCH){
                Matrix<N,N,T> L_j;
                Matrix<N,N,T> L_j_inv;
                Matrix<N,N,T> L;
                Matrix<N,N,T> P;
                Matrix<N,N,T> P_j;
                Matrix<N,N,T> R(A);
                
                for(int j = 0; j < N-1; j++){
                    if(!gaussStep(A, b, L_j, L_j_inv, P_j, j, search_algo)){
                        return;
                    }
                    L = P_j * L * P_j * L_j_inv;
                    P = P_j * P;
                }
                this->R = A;
                this->L = L;
                this->P = P;
            }
           LRCombination(){};
           Matrix<N,1,T> solve(Matrix<N,1,T> b){
               Matrix<N,1,T> z;
               Matrix<N,1,T> x;
               z = forSub(L, P * b);
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
                if(!gaussStep(A, b, L_j, L_j_inv, P_j, j, search_algo)){
                    return old_b;    
                };
            }
            x = backSub(A, b);
            return x;
        }


}

