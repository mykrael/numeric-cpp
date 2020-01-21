#include <iostream>
using namespace std;
#include "matrix.h"

int main(){
    float a[2][2] = {{1,2},{3,4}};
    Matrix<2,2,float> m(a);
    Matrix<2,2,float> n = m;
    cout << "n = m = " << n << endl;
    m[0][0] = 11;
    cout << "changed first element of m. n should remain the same: " << n << " but m is different: " << m << endl;
    n[0][1] = 12;
    cout << "changed first element of n. m should remain the same: " << m << " but m is different: " << n << endl;
    
    
    float b[2][3] = {{1,2,3},{4,5,6}};
    float c[3][2] = {{1,2},{3,4},{5,6}};
    
    Matrix<2,3,float> K(b);
    Matrix<3,2,float> L(c);

    m = K*L;
    K.print();
    L.print();
    m.print();
}
