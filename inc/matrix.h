#pragma once
#include <iostream>
#include <array>


template<int COLS, class T>
class Row{
private:
    T r[COLS];
public:
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
        for(auto row: *this){
            for(auto &e: row){
                e = T(0);
            }
        }
    }

    Matrix(T a[ROWS][COLS]){
        for(int row = 0; row < ROWS; row++){
            m[row] = Row<COLS, T>(a[row]);
        }

    }

   Row<COLS, T>& operator[](int row){
        return m[row];
    };
    iterator begin(){return &m[0];}
    iterator end(){return &m[ROWS];}
    void print(){
        for(auto row: m){
            for(auto e: row){
                std::cout << e << " ";
            }
            std::cout << std::endl;
        }
    }


    // Basic Calculation
    Matrix operator+(Matrix& other){
        T res[ROWS][COLS];
        for(int row = 0; row < ROWS; row++){
            for(int col = 0; col < COLS; col++){
                res[row][col] = other[row][col] + (*this)[row][col];
            }
        }
        return Matrix(res);
    }
    Matrix operator*(T other){
        T res[ROWS][COLS];
        for(int row = 0; row < ROWS; row++){
            for(int col = 0; col < COLS; col++){
                res[row][col] = other * (*this)[row][col];
            }
        }
        return Matrix(res);
   }

   template<int OTHER_COLS>
   Matrix<ROWS, OTHER_COLS, T> operator*(Matrix<COLS, OTHER_COLS, T>& other){
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

   friend Matrix operator*(T scalar, Matrix& mat){
       Matrix result = mat * scalar;
       return result;
   }
        
   friend std::ostream & operator << (std::ostream & out, Matrix mat){
        out << "[";
        for(auto row: mat){
            out << row << ";";
        }
        out << " ]";
        return out;
    }
};
