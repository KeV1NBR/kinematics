#include<cstdio>
#include<cstdlib>
#include<cmath>
#include <iostream>

#include <Eigen/Core>
#include <Eigen/LU>

#include "kinematics.hpp"

using namespace std;
using namespace Eigen;

KM::KM(VectorXd inBASE, VectorXd inTOOL, VectorXd inALPHA, VectorXd inD, VectorXd inA) : JXdeg(6), Position(6), ALPHA(6), D(6), A(6)

{

    this-> ALPHA = inALPHA;
    this-> D     = inD;
    this-> A     = inA;
    
    this-> BASE  = DH_deal(inBASE);
    this->TOOL  << DH_deal(inTOOL);

   // JXdeg    << 0, 0, 0, 0, 0, 0;
   // Position << 0, 0, 0, 0, 0, 0;
    
    DH_J1 << to_rad(JXdeg(0)), to_rad(ALPHA(0)), D(0), A(0);
    DH_J2 << to_rad(JXdeg(1)), to_rad(ALPHA(1)), D(1), A(1);
    DH_J3 << to_rad(JXdeg(2)), to_rad(ALPHA(2)), D(2), A(2);
    DH_J4 << to_rad(JXdeg(3)), to_rad(ALPHA(3)), D(3), A(3);
    DH_J5 << to_rad(JXdeg(4)), to_rad(ALPHA(4)), D(4), A(4);
    DH_J6 << to_rad(JXdeg(5)), to_rad(ALPHA(5)), D(5), A(5);

}

VectorXd KM::forwardKinematics(VectorXd inJXdeg)
{   
    JXdeg = inJXdeg;

    DH_J1(0) =to_rad(JXdeg(0));        
    DH_J2(0) =to_rad(JXdeg(1));        
    DH_J3(0) =to_rad(JXdeg(2)-90);        
    DH_J4(0) =to_rad(JXdeg(3));        
    DH_J5(0) =to_rad(JXdeg(4));        
    DH_J6(0) =to_rad(JXdeg(5)+180);

    J1 = DH_deal(DH_J1);
    J2 = DH_deal(DH_J2);
    J3 = DH_deal(DH_J3);
    J4 = DH_deal(DH_J4);
    J5 = DH_deal(DH_J5);
    J6 = DH_deal(DH_J6);

    matrixPrint(J1);
    matrixPrint(J2);
    matrixPrint(J3);
    matrixPrint(J4);
    matrixPrint(J5);
    matrixPrint(J6);

    R01= BASE * J1;
    R02= R01  * J2;
    R03= R02  * J3;
    R04= R03  * J4;
    R05= R04  * J5;
    R06= R05  * J6;
    R0T= R06  * TOOL;
    
    matrixPrint(R01);
    matrixPrint(R02);
    matrixPrint(R03);
    matrixPrint(R04);
    matrixPrint(R05);
    matrixPrint(R06);
    matrixPrint(R0T);

    double p = atan2(sqrt( (R0T(0,2)*R0T(0,2)) + (R0T(1,2)*R0T(1,2))), -R0T(2,2));
 
             Position  <<  R0T(0,3), R0T(1,3), R0T(2,3), to_degree(atan2(R0T(2,0)/p , R0T(2,1)/p)), to_degree(p), to_degree(atan2(R0T(0,2)/p , R0T(1,2)/p)); 
  
    return Position;
}

VectorXd KM::inverseKinematics(VectorXd inPosition)
{
    Position = inPosition;

    Matrix4d R0T_off;
    Matrix4d TOOL_inv;
    Matrix4d R06_rm;
    Matrix4d R03_trans;
    VectorXd J23_FWD(10);
    VectorXd J23_MID(10);
    Matrix4d R36;

    R0T << cos(to_rad(Position(3)))* cos(to_rad(Position(5)))- cos(to_rad(Position(4)))* sin(to_rad(Position(3)))* sin(to_rad(Position(5))),
           cos(to_rad(Position(5)))* sin(to_rad(Position(3)))+ cos(to_rad(Position(3)))* cos(to_rad(Position(4)))* sin(to_rad(Position(5))),
           sin(to_rad(Position(4)))* sin(to_rad(Position(5))),
           Position(0),

           cos(to_rad(Position(4)))* cos(to_rad(Position(5)))* sin(to_rad(Position(3)))+ cos(to_rad(Position(3)))* sin(to_rad(Position(5))),
           cos(to_rad(Position(3)))* cos(to_rad(Position(4)))* cos(to_rad(Position(5)))- sin(to_rad(Position(3)))* sin(to_rad(Position(5))),
           cos(to_rad(Position(5)))* sin(to_rad(Position(4))),
           Position(1),

           sin(to_rad(Position(3)))* sin(to_rad(Position(4))),
           cos(to_rad(Position(3)))* sin(to_rad(Position(4))),
           -cos(to_rad(Position(4))),
           Position(2),
            
           0, 0, 0, 1;

    R06_rm << -1, 0, 0, 0,
              -sin(to_rad(180))*cos(DH_J6(1)), cos(to_rad(180))*cos(DH_J6(1)), sin(DH_J6(1)), 0,
               sin(to_rad(180))*sin(DH_J6(1)), cos(to_rad(180))*sin(DH_J6(1)), cos(DH_J6(1)), -DH_J6(2),
               0, 0, 0, 1;


    R0T_off = BASE * R0T;
    R0T_off(0,0) = R0T_off(0,0) * -1;
    TOOL_inv= TOOL.inverse();
    R06 = R0T_off * TOOL_inv;
    R05 = R06*R06_rm;

    if(Position(0)>=0 && Position(1)>=0)
        DH_J1(0) = atan((R05(1,3))/(R05(0,3)));
    else if(Position(0)>0 && Position(1)<0)
        DH_J1(0) = atan((R05(1,3))/(R05(0,3)));
    else if(Position(0)<0 && Position(1)<0)
        DH_J1(0) = atan((R05(1,3))/(R05(0,3)))-to_rad(180);
    else if(Position(0)<0 && Position(1)>0)
        DH_J1(0) = atan((R05(1,3))/(R05(0,3)))+to_rad(180);

       J23_FWD(0)= sqrt((abs(R05(1, 3)) * abs(R05(1, 3))) + abs(R05(0, 3)) * abs(R05(0, 3))); 
       J23_FWD(1)= R05(2, 3) - DH_J1(2);
       J23_FWD(2)= J23_FWD(0) - DH_J1(3); 
       J23_FWD(3)= sqrt((J23_FWD(1) * J23_FWD(1)) + (J23_FWD(2) * J23_FWD(2)));
       J23_FWD(4)= atan(J23_FWD(1) / J23_FWD(2)); 
       J23_FWD(5)= acos(((DH_J2(3) * DH_J2(3)) + (J23_FWD(3) * J23_FWD(3)) - (abs(DH_J4(2)) * abs(DH_J4(2)))) /(2 * DH_J2(3) * J23_FWD(3)));
       J23_FWD(6)= to_rad(180) - acos(((abs(DH_J4(2))) * (abs(DH_J4(2))) + (DH_J2(3) * DH_J2(3)) - (J23_FWD(3) * J23_FWD(3))) / 
                          (2 * abs(DH_J4(2)) * (DH_J2(3))));
       J23_FWD(7)= 0;
       J23_FWD(8)= -J23_FWD(4) - J23_FWD(5);
       J23_FWD(9)= J23_FWD(6);


       J23_MID(0)= J23_FWD(0);
       J23_MID(1)= J23_FWD(1);
       J23_MID(2)= -J23_FWD(2);
       J23_MID(3)= sqrt((J23_MID(1)*J23_MID(1))+(J23_MID(2)*J23_MID(2)));
       J23_MID(4)= acos(((DH_J2(3)*DH_J2(3))+(J23_MID(3)*J23_MID(3))-(abs(DH_J4(2))*abs(DH_J4(2))))/(2*DH_J2(3)*J23_MID(3)));
       J23_MID(5)= atan(J23_MID(2)/J23_MID(1));
       J23_MID(6)= to_rad(180)-acos(((abs(DH_J4(2)))*(abs(DH_J4(2)))+(DH_J2(3)*DH_J2(3))-(J23_MID(3)*J23_MID(3)))/(2*abs(DH_J4(2))*(DH_J2(3))));
       J23_MID(7)= to_rad(90)-J23_MID(4)-J23_MID(5);
       J23_MID(8)= to_rad(-180)+J23_MID(7);
       J23_MID(9)= J23_MID(6);
    
    if(J23_FWD(2)>=0)
    {
        DH_J2(0)=J23_FWD(8);
        DH_J3(0)=J23_FWD(9)-to_rad(90);
    }
    else
    {
        DH_J2(0)=J23_MID(8);
        DH_J3(0)=J23_MID(9);
    }

    J1 = DH_deal(DH_J1);
    J2 = DH_deal(DH_J2);
    J3 = DH_deal(DH_J3);

    R01 = BASE * J1;
    R02 = R01  * J2;
    R03 = R02  * J3;
    R03_trans = R03.transpose();
    R36 = R03_trans * R05;

    DH_J5(0) = atan2( sqrt(1-(R36(2,2)*R36(2,2))), R36(2,2));//
    if(JXdeg(4)>=0 && DH_J5(0)>=0)
    {
        DH_J4(0) = atan2(R36(1,2), R36(0,2)); //
        if(R36(2,1)<0)
            DH_J6(0) = atan2(-R36(2,1),R36(2,0))-to_rad(180);//
        else
            DH_J6(0) = atan2(-R36(2,1),R36(2,0))+to_rad(180);//

    }

    else
    {
        DH_J5(0) = atan2(-sqrt(1-(R36(2,2)*R36(2,2))), R36(2,2)); //
        if(DH_J5(0)>=0)
            DH_J4(0) = atan2(R36(1,2), R36(0,2));//
        else
            DH_J4(0) = atan2(-R36(1,2), -R36(0,2));//

        if(R36(2,1)<0)
            DH_J6(0) = atan2(-R36(2,1),R36(2,0))-to_rad(180);//
        else
            DH_J6(0) = atan2(-R36(2,1),R36(2,0))+to_rad(180);//
            
    }

    JXdeg  <<  to_degree(DH_J1(0)),to_degree(DH_J2(0)), to_degree(DH_J3(0)+to_rad(90)), to_degree(DH_J4(0)), to_degree(DH_J5(0)), to_degree(DH_J6(0));
    return JXdeg;
}

double KM::to_rad(double degree)
{
    return degree*M_PI/180;
}

double KM::to_degree(double rad)
{
    return rad*180/M_PI;
}

Matrix4d KM::DH_deal(Vector4d DH_JX)
{
    Matrix4d JX; 
             JX << cos(DH_JX(0)), -sin(DH_JX(0))*cos(DH_JX(1)), sin(DH_JX(0))*sin(DH_JX(1)) , DH_JX(3)*cos(DH_JX(0)),
                   sin(DH_JX(0)), cos(DH_JX(0))*cos(DH_JX(1)) , -cos(DH_JX(0))*sin(DH_JX(1)), DH_JX(3)*sin(DH_JX(0)),
                   0            , sin(DH_JX(1))               , cos(DH_JX(1))               , DH_JX(2)              ,
                   0            , 0                           , 0                           , 1                     ;

    return JX;
}

void KM::matrixPrint(Matrix4d matrix)
{
    int i,j;

    for(i=0; i<4;i++)
    {
        for(j=0;j<4;j++)
        {
            printf("%lf\t",matrix(i,j));
        }
        printf("\n");
    }
    printf("\n");
}

void KM:: vectorPrint(VectorXd vec)
{
    int i;
    for(i=0;i<10;i++)
        printf("%lf\n",vec(i));

    printf("\n\n");

}
