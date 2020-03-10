#include <iostream>
using namespace std;
#include "matrix.h"
#include "gauss-elimination.h"

int main(){
    Matrix<3,3,float> A;
    Matrix<3,1,float> b;
    Matrix<3,1,float> x;
    

    cout << "Enter a regular 3x3 matrix: " << endl;
    A.input();
    cout << "Enter a 3x1 vector as the right side of the equation: " << endl;
    b.input();

    cout << "---------------------------------" << endl;
    x = GAUSS::solve(A, b);
    

    cout << "Solution: " << endl;
    x.print();

    cout << "matrix * solution = right side:" << endl;
    (A*x).print();

    GAUSS::LRCombination<3,float> LRCom(A,b);
    cout << "Result of LR-Decomposition:" << endl;
    cout << "L: " << endl;
    LRCom.L.print();

    cout << "R: " << endl;
    LRCom.R.print();

    cout << "solution: " << endl;
    x = LRCom.solve(b);
    x.print();

    cout << "P*L*R : " << endl;

    (~LRCom.P*LRCom.L*LRCom.R).print();
    

}
