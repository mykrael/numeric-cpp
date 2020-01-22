#include <iostream>
using namespace std;
#include "matrix.h"
#include "gauss-elimination.h"

int main(){
    Matrix<3,3,float> A; 
    Matrix<3,3,float> B;
    cout << "Input A: " << endl;
    A.input();
    A.print();
    cout << "Input B: " << endl;
    B.input();
    B.print();
    
    A = A - B;
    A.print();
}
