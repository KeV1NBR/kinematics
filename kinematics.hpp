#include <Eigen/Core>

class KM
{
public:
    
    // Member functions
    KM(Eigen::VectorXd inBASE, Eigen::VectorXd inTOOL, Eigen::VectorXd inALPHA, Eigen::VectorXd inD, Eigen::VectorXd inA);
    Eigen::VectorXd forwardKinematics();
    Eigen::VectorXd inverseKinematics();
    // Member variables
    Eigen::VectorXd JXdeg = Eigen::VectorXd(6);
    Eigen::VectorXd Position = Eigen::VectorXd(6);

protected:
    
    Eigen::Vector4d DH_J1;
    Eigen::Vector4d DH_J2;
    Eigen::Vector4d DH_J3; 
    Eigen::Vector4d DH_J4; 
    Eigen::Vector4d DH_J5; 
    Eigen::Vector4d DH_J6; 
        
    Eigen::Matrix4d J1;
    Eigen::Matrix4d J2;
    Eigen::Matrix4d J3;
    Eigen::Matrix4d J4;
    Eigen::Matrix4d J5;
    Eigen::Matrix4d J6;

    Eigen::Matrix4d R01;
    Eigen::Matrix4d R02;
    Eigen::Matrix4d R03;
    Eigen::Matrix4d R04;
    Eigen::Matrix4d R05;
    Eigen::Matrix4d R06;
    Eigen::Matrix4d R0T;

    void matrixPrint(Eigen::Matrix4d matrix);
    void vectorPrint(Eigen::VectorXd vec);

private:

    Eigen::VectorXd ALPHA = Eigen::VectorXd(6); 
    Eigen::VectorXd D = Eigen::VectorXd(6); 
    Eigen::VectorXd A = Eigen::VectorXd(6);

    Eigen::Matrix4d BASE;
    Eigen::Matrix4d TOOL;    
    Eigen::Matrix4d DH_deal(Eigen::Vector4d DH_JX);
    
    static double to_rad(double degree);
    static double to_degree(double rad);
};
