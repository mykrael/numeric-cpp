#include <iostream>
using namespace std;
#include "matrix.h"
#include "gauss-elimination.h"

int main(){
    Matrix<3,3,float> A; 
    Matrix<3,1,float> b;
    cout << "Input A: " << endl;
    A.input();
    A.print();
    cout << "Input B: " << endl;
    b.input();
    b.print();
    
    Matrix<3,1,float> x;
    x = GAUSS::solve(A, b);
    cout << "Solved x: " << endl;
    x.print();
}
