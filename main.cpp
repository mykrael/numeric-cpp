#include <iostream>
using namespace std;
#include "matrix.h"
#include "gauss-elimination.h"

int main(){
    Matrix<3,2,float> A; 
    cout << "Input A:" << endl;
    A.input();
    Matrix<3,3,float> result;
    result = A * ~A;
    result.print();
}
