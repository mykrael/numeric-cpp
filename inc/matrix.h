#pragma once
#include <iostream>
#include <array>


template<int COLS, class T>
class Row{
public:
    T r[COLS];
    typedef T * iterator;
    typedef const T * const_interator;
    Row(){};
    Row(T a[COLS]){
        for(int col = 0; col < COLS; col++){
           r[col] = a[col];
        }
    }
    T& operator[](int col){
        return r[col];
    }
    iterator begin(){return &r[0];}
    iterator end(){return &r[COLS];}
    friend std::ostream & operator << (std::ostream & out, Row row){
        for(auto e: row){
            out << " " << e;
        }
        return out;
    }
};


template<int ROWS, int COLS, class T>
class Matrix{
private:
    Row<COLS, T> m[ROWS];
public:
    typedef Row<COLS, T>* iterator;
    typedef const Row<COLS, T>* const_iterator;
    
    Matrix(){
        int count = 0;
        for(auto& row: *this){
            for(auto &e: row){
                if(count/COLS == count%COLS){
                    e = T(1);
                }else{
                    e = T(0);
                }
                count++;
            }
        }
    }

    Matrix(const Matrix& mat){
        for(int row = 0; row < ROWS; row++){
            for(int col = 0; col < COLS; col++){
                (*this)[row][col] = mat.m[row].r[col];
            }
        }
    }

    Matrix(T a[ROWS][COLS]){
        for(int row = 0; row < ROWS; row++){
            m[row] = Row<COLS, T>(a[row]);
        }

    }

   template <typename... Args>
      Matrix(Args... args){
           T a[] = {args...};
           int count = 0;
           while(count < ROWS * COLS){
               for(auto e : a){
                   if(count < ROWS * COLS){
                       (*this)[count/COLS][count%COLS] = e;
                       count++;
                   }else{
                       break;
                   }
               }
           }
       }

    Row<COLS, T>& operator[](int row){
        return m[row];
    };
    iterator begin(){return &m[0];}
    iterator end(){return &m[ROWS];}

    Matrix& operator=(Matrix other){
        for(int row = 0; row < ROWS; row++){
            for(int col = 0; col < COLS; col++){
                (*this)[row][col] = other[row][col];
            }
        }
        return (*this);
    }

    // Printing
    void print(){
        int counter = 0;
        for(auto row: m){
            if(counter == 0){
                std::cout << "\n[";
            }else{
                std::cout << " ";
            }
            for(auto e: row){
                counter++;
                std::cout << e;
                if(counter == ROWS * COLS){
                    std::cout << "]";
                }else{
                    std::cout << "\t\t";
                }
            }
            std::cout <<  std::endl;
        }
        std::cout << "\n" <<  std::endl;
    }

   friend std::ostream & operator << (std::ostream & out, Matrix mat){
        out << "[";
        for(auto row: mat){
            out << row << ";";
        }
        out << " ]";
        return out;
    }

    // Inputing
    void input(){
        int count = 0;
        T i;
        std::cin >> i;
        while(!std::cin.fail() && count < ROWS * COLS){
            (*this)[count/COLS][count%COLS] = i;
            count++;
            if(count%COLS != 0){
                std::cout << "\x1b[A";
            }
            for(int tabs = 0; tabs < count%COLS; tabs++){
                std::cout << "\t";
            }

            if(count < ROWS * COLS){
                std::cin >> i;
            }
        }
        std::cout << std::endl;
       

    }

    // Basic Calculation
    Matrix operator+(Matrix other){
    T res[ROWS][COLS];
        for(int row = 0; row < ROWS; row++){
            for(int col = 0; col < COLS; col++){
                res[row][col] = other[row][col] + (*this)[row][col];
            }
        }
        Matrix result(res);
        return result;
    }

    Matrix operator*(T other){
        T res[ROWS][COLS];
        for(int row = 0; row < ROWS; row++){
            for(int col = 0; col < COLS; col++){
                res[row][col] = other * (*this)[row][col];
            }
        }
        Matrix result(res);
        return result;
   }

   template<int OTHER_COLS>
   Matrix<ROWS, OTHER_COLS, T> operator*(Matrix<COLS, OTHER_COLS, T>  other){
       T result[ROWS][OTHER_COLS];
       for(int row = 0; row < ROWS; row++){
           for(int col = 0; col < OTHER_COLS; col++){
               T sum = T(0);
               for(int elem = 0; elem < COLS; elem++){
                   sum += (*this)[row][elem] * other[elem][col];
               }
               result[row][col] = sum;
           }
       }
       Matrix<ROWS, OTHER_COLS, T> m(result);
       return m;
   }

   Matrix operator-(Matrix other){
       Matrix result;
       result = (*this) + T(-1) * other;
       return result;
   }

   friend Matrix operator*(T scalar, Matrix mat){
       Matrix result;      
       result = mat * scalar;
       return result;
   }

   // Transpose
   Matrix<COLS, ROWS, T> operator~(){
       Matrix<COLS, ROWS, T> result;
       for(int row = 0; row < ROWS; row++){
           for(int col = 0; col < COLS; col++){
               result[col][row] = (*this)[row][col];
           }
       }
       return result;
   }

        
};
