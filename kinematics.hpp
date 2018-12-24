#include <cstdio>
#include <cstdlib>
#include <string>
#include <cmath>
#include <Eigen/Core>
#include <iostream>

using namespace std;
using namespace Eigen;

class KM
{
public:
    
    // Member functions
    KM();
    VectorXd forwardKinematics();
    VectorXd inverseKinematics();
    // Member variables
    VectorXd JXdeg = VectorXd(6);
    VectorXd Position =VectorXd(6);

protected:
    
    Vector4d DH_J1;
    Vector4d DH_J2;
    Vector4d DH_J3; 
    Vector4d DH_J4; 
    Vector4d DH_J5; 
    Vector4d DH_J6; 
        
    Matrix4d J1;
    Matrix4d J2;
    Matrix4d J3;
    Matrix4d J4;
    Matrix4d J5;
    Matrix4d J6;

    Matrix4d R01;
    Matrix4d R02;
    Matrix4d R03;
    Matrix4d R04;
    Matrix4d R05;
    Matrix4d R06;
    Matrix4d R0T;

    void matrixPrint(Matrix4d matrix);
    void vectorPrint(VectorXd vec);

private:

    VectorXd ALPHA = VectorXd(6); 
    VectorXd D = VectorXd(6); 
    VectorXd A = VectorXd(6);

    Matrix4d BASE;
    Matrix4d TOOL;    
    Matrix4d DH_deal(Vector4d DH_JX);
    
    double to_rad(double degree);
    double to_degree(double rad);
};
