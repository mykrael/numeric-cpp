#include <iostream>
using namespace std;
#include "matrix.h"
#include "gauss-elimination.h"

int main(){
    Matrix<3,3,double> A(2.0,2.0,2.0,4.0,5.0,6.0,7.0,8.0,11.0);
    Matrix<3,1,double> x(5.0,7.0,9.0);
    Matrix<3,1,double> b;
    b = A * x;
    x.print();
    A.print();
    b.print();
    cout << "_______________________________" << endl;
    x = GAUSS::solve(A, b);
    x.print();
    
}
