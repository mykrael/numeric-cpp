#include <iostream>
using namespace std;
#include "matrix.h"

int main(){
    Matrix<3,3,int> A;
    Matrix<3,3,int> B;
    Matrix<3,3,int> res;
    res = A + A + 3 * A * B;
    res.print();
}
