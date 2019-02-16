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

    JXdeg    << 0.0001, -90, 90, 0.0001, 90, 0.0001;
    Position << 286.83, 0.000437, 438.52, 0.0001, 179.9999, 0.0001;
    
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

    //matrixPrint(J1);
    //matrixPrint(J2);
    //matrixPrint(J3);
    //matrixPrint(J4);
    //matrixPrint(J5);
    //matrixPrint(J6);

    R01= BASE * J1;
    R02= R01  * J2;
    R03= R02  * J3;
    R04= R03  * J4;
    R05= R04  * J5;
    R06= R05  * J6;
    R0T= R06  * TOOL;
    
    //matrixPrint(R01);
    //matrixPrint(R02);
    //matrixPrint(R03);
    //matrixPrint(R04);
    //matrixPrint(R05);
    //matrixPrint(R06);
    //matrixPrint(R0T);

    double p = atan2(sqrt( (R0T(0,2)*R0T(0,2)) + (R0T(1,2)*R0T(1,2))), -R0T(2,2));
 
             Position  <<  R0T(0,3), R0T(1,3), R0T(2,3), to_degree(atan2(R0T(2,0)/p , R0T(2,1)/p)), to_degree(p), to_degree(atan2(R0T(0,2)/p , R0T(1,2)/p)); 
  
    return Position;
}

VectorXd KM::inverseKinematics(VectorXd inPosition)
{
    double O4 = inPosition(0);
    double O5 = inPosition(1);
    double O6 = inPosition(2);
    double O7 = inPosition(3);
    double O8 = inPosition(4);
    double O9 = inPosition(5);
    double V8 = 'F';

    Position = inPosition;

 if (this->JXdeg(0) == 0)
    this->JXdeg(0) = .0001;
  if (this->JXdeg(1) == 0)
    this->JXdeg(1) = .0001;
  if (this->JXdeg(2) == 0)
    this->JXdeg(2) = .0001;
  if (this->JXdeg(3) == 0)
    this->JXdeg(3) = .0001;
  if (this->JXdeg(4) == 0)
    this->JXdeg(4) = .0001;
  if (this->JXdeg(5) == 0)
    this->JXdeg(5) = .0001;

  if (O4 == 0)
    O4 = .0001;
  if (O5 == 0)
    O5 = .0001;
  if (O6 == 0)
    O6 = .0001;
  if (O7 == 0)
    O7 = .0001;
  if (O8 == 0)
    O8 = .0001;
  if (O9 == 0)
    O9 = .0001;
  #quadrant
  if (O4>0 && O5>0)
    double V9 = 1;
  else if (O4>0 && O5<0)
    doubleV9 = 2;
  else if (O4<0 && O5<0)
    double V9 = 3;
  else if (O4<0 and O5>0)
    double V9 = 4;
  // DH TABLE
  double D13 = to_rad(ALPHA(0));
  double D14 = to_rad(ALPHA(1));
  double D15 = to_rad(ALPHA(2));
  double D16 = to_rad(ALPHA(3));
  double D17 = to_rad(ALPHA(4));
  double D18 = to_rad(ALPHA(5));
  double E13 = D(0);
  double E14 = D(1);
  double E15 = D(2);
  double E16 = D(3);
  double E17 = D(4);
  double E18 = D(5);
  double F13 = A(0);
  double F14 = A(1);
  double F15 = A(2);
  double F16 = A(3);
  double F17 = A(4);
  double F18 = A(5);
 // WORK FRAME INPUT
  double H13 = 0; 
  double H14 = 0; 
  double H15 = 0; 
  double H16 = 0; 
  double H17 = 0; 
  double H18 = 0; 
  // TOOL FRAME INPUT
  double J13 = 0; 
  double J14 = 0; 
  double J15 = 0; 
  double J16 = 0; 
  double J17 = 0; 
  double J18 = 0; 
  // WORK FRAME TABLE
  double N30 = cos(to_rad(H18))*cos(to_rad(H17));
  double O30 = -sin(to_rad(H18))*cos(to_rad(H16))+cos(to_rad(H18))*sin(to_rad(H17))*sin(to_rad(H16));
  double P30 = sin(to_rad(H18))*sin(to_rad(H16))+cos(to_rad(H18))*sin(to_rad(H17))*cos(to_rad(H16));
  double Q30 = H13;
  double N31 = sin(to_rad(H18))*cos(to_rad(H17));
  double O31 = cos(to_rad(H18))*cos(to_rad(H16))+sin(to_rad(H18))*sin(to_rad(H17))*sin(to_rad(H16));
  double P31 = -cos(to_rad(H18))*sin(to_rad(H16))+sin(to_rad(H18))*sin(to_rad(H17))*cos(to_rad(H16));
  double Q31 = H14;
  double N32 = -sin(to_rad(H18));
  double O32 = cos(to_rad(H17))*sin(to_rad(H16));
  double P32 = cos(to_rad(H17))*cos(to_rad(H16));
  double Q32 = H15;
  double N33 = 0;
  double O33 = 0;
  double P33 = 0;
  double Q33 = 1;
  // R 0-T
  double X30 = cos(to_rad(O7))*cos(to_rad(O9))-cos(to_rad(O8))*sin(to_rad(O7))*sin(to_rad(O9));
  double Y30 = cos(to_rad(O9))*sin(to_rad(O7))+cos(to_rad(O7))*cos(to_rad(O8))*sin(to_rad(O9));
  double Z30 = sin(to_rad(O8))*sin(to_rad(O9));
  double AA30 = O4;
  double X31 = cos(to_rad(O8))*cos(to_rad(O9))*sin(to_rad(O7))+cos(to_rad(O7))*sin(to_rad(O9));
  double Y31 = cos(to_rad(O7))*cos(to_rad(O8))*cos(to_rad(O9))-sin(to_rad(O7))*sin(to_rad(O9));
  double Z31 = cos(to_rad(O9))*sin(to_rad(O8));
  double AA31 = O5;
  double X32 = sin(to_rad(O7))*sin(to_rad(O8));
  double Y32 = cos(to_rad(O7))*sin(to_rad(O8));
  double Z32 = -cos(to_rad(O8));
  double AA32 = O6;
  double X33 = 0;
  double Y33 = 0;
  double Z33 = 0;
  double AA33 = 1;
  // R 0-T   offset by work frame
  double X36 = ((N30*X30)+(O30*X31)+(P30*X32)+(Q30*X33))*-1;
  double Y36 = (N30*Y30)+(O30*Y31)+(P30*Y32)+(Q30*Y33);
  double Z36 = (N30*Z30)+(O30*Z31)+(P30*Z32)+(Q30*Z33);
  double AA36 = (N30*AA30)+(O30*AA31)+(P30*AA32)+(Q30*AA33);
  double X37 = (N31*X30)+(O31*X31)+(P31*X32)+(Q31*X33);
  double Y37 = (N31*Y30)+(O31*Y31)+(P31*Y32)+(Q31*Y33);
  double Z37 = (N31*Z30)+(O31*Z31)+(P31*Z32)+(Q31*Z33);
  double AA37 = (N31*AA30)+(O31*AA31)+(P31*AA32)+(Q31*AA33);
  double X38 = (N32*X30)+(O32*X31)+(P32*X32)+(Q32*X33);
  double Y38 = (N32*Y30)+(O32*Y31)+(P32*Y32)+(Q32*Y33);
  double Z38 = (N32*Z30)+(O32*Z31)+(P32*Z32)+(Q32*Z33);
  double AA38 = (N32*AA30)+(O32*AA31)+(P32*AA32)+(Q32*AA33);
  double X39 = (N33*X30)+(O33*X31)+(P33*X32)+(Q33*X33);
  double Y39 = (N33*Y30)+(O33*Y31)+(P33*Y32)+(Q33*Y33);
  double Z39 = (N33*Z30)+(O33*Z31)+(P33*Z32)+(Q33*Z33);
  double AA39 = (N33*AA30)+(O33*AA31)+(P33*AA32)+(Q33*AA33);
  // TOOL FRAME
  double X42 = cos(to_rad(J18))*cos(to_rad(J17));
  double Y42 = -sin(to_rad(J18))*cos(to_rad(J16))+cos(to_rad(J18))*sin(to_rad(J17))*sin(to_rad(J16));
  double Z42 = sin(to_rad(J18))*sin(to_rad(J16))+cos(to_rad(J18))*sin(to_rad(J17))*cos(to_rad(J16));
  double AA42 = (J13);
  double X43 = sin(to_rad(J18))*cos(to_rad(J17));
  double Y43 = cos(to_rad(J18))*cos(to_rad(J16))+sin(to_rad(J18))*sin(to_rad(J17))*sin(to_rad(J16));
  double Z43 = -cos(to_rad(J18))*sin(to_rad(J16))+sin(to_rad(J18))*sin(to_rad(J17))*cos(to_rad(J16));
  double AA43 = (J14);
  double X44 = -sin(to_rad(J18));
  double Y44 = cos(to_rad(J17))*sin(to_rad(J16));
  double Z44 = cos(to_rad(J17))*cos(to_rad(J16));
  double AA44 = (J15);
  double X45 = 0;
  double Y45 = 0;
  double Z45 = 0;
  double AA45 = 1;
  // INVERT TOOL FRAME
  double X48 = X42;
  double Y48 = X43;
  double Z48 = X44;
  double AA48 = (X48*AA42)+(Y48*AA43)+(Z48*AA44);
  double X49 = Y42;
  double Y49 = Y43;
  double Z49 = Y44;
  double AA49 = (X49*AA42)+(Y49*AA43)+(Z49*AA44);
  double X50 = Z42;
  double Y50 = Z43;
  double Z50 = Z44;
  double AA50 = (X50*AA42)+(Y50*AA43)+(Z50*AA44);
  double X51 = 0;
  double Y51 = 0;
  double Z51 = 0;
  double AA51 = 1;
  // R 0-6
  double X54 =(X36*X48)+(Y36*X49)+(Z36*X50)+(AA36*X51);
  double Y54 =(X36*Y48)+(Y36*Y49)+(Z36*Y50)+(AA36*Y51);
  double Z54 =(X36*Z48)+(Y36*Z49)+(Z36*Z50)+(AA36*Z51);
  double AA54 =(X36*AA48)+(Y36*AA49)+(Z36*AA50)+(AA36*AA51);
  double X55 =(X37*X48)+(Y37*X49)+(Z37*X50)+(AA37*X51);
  double Y55 =(X37*Y48)+(Y37*Y49)+(Z37*Y50)+(AA37*Y51);
  double Z55 =(X37*Z48)+(Y37*Z49)+(Z37*Z50)+(AA37*Z51);
  double AA55 =(X37*AA48)+(Y37*AA49)+(Z37*AA50)+(AA37*AA51);
  double X56 =(X38*X48)+(Y38*X49)+(Z38*X50)+(AA38*X51);
  double Y56 =(X38*Y48)+(Y38*Y49)+(Z38*Y50)+(AA38*Y51);
  double Z56 =(X38*Z48)+(Y38*Z49)+(Z38*Z50)+(AA38*Z51);
  double AA56 =(X38*AA48)+(Y38*AA49)+(Z38*AA50)+(AA38*AA51);
  double X57 =(X39*X48)+(Y39*X49)+(Z39*X50)+(AA39*X51);
  double Y57 =(X39*Y48)+(Y39*Y49)+(Z39*Y50)+(AA39*Y51);
  double Z57 =(X39*Z48)+(Y39*Z49)+(Z39*Z50)+(AA39*Z51);
  double AA57 =(X39*AA48)+(Y39*AA49)+(Z39*AA50)+(AA39*AA51);
  // REMOVE R 0-6
  double X60 =cos(to_rad(180));
  double Y60 =sin(to_rad(180));
  double Z60 = 0;
  double AA60 = 0;
  double X61 =-sin(to_rad(180))*cos(D18);
  double Y61 =cos(to_rad(180))*cos(D18);
  double Z61 =sin(D18);
  double AA61 = 0;
  double X62 =sin(to_rad(180))*sin(D18);
  double Y62 =-cos(to_rad(180))*sin(D18);
  double Z62 =cos(D18);
  double AA62 = -E18;
  double X63 = 0;
  double Y63 = 0;
  double Z63 = 0;
  double AA63 = 1;
  // R 0-5 (center spherica wrist)
  double X66 =(X54*X60)+(Y54*X61)+(Z54*X62)+(AA54*X63);
  double Y66 =(X54*Y60)+(Y54*Y61)+(Z54*Y62)+(AA54*Y63);
  double Z66 =(X54*Z60)+(Y54*Z61)+(Z54*Z62)+(AA54*Z63);
  double AA66 =(X54*AA60)+(Y54*AA61)+(Z54*AA62)+(AA54*AA63);
  double X67 =(X55*X60)+(Y55*X61)+(Z55*X62)+(AA55*X63);
  double Y67 =(X55*Y60)+(Y55*Y61)+(Z55*Y62)+(AA55*Y63);
  double Z67 =(X55*Z60)+(Y55*Z61)+(Z55*Z62)+(AA55*Z63);
  double AA67 =(X55*AA60)+(Y55*AA61)+(Z55*AA62)+(AA55*AA63);
  double X68 =(X56*X60)+(Y56*X61)+(Z56*X62)+(AA56*X63);
  double Y68 =(X56*Y60)+(Y56*Y61)+(Z56*Y62)+(AA56*Y63);
  double Z68 =(X56*Z60)+(Y56*Z61)+(Z56*Z62)+(AA56*Z63);
  double AA68 =(X56*AA60)+(Y56*AA61)+(Z56*AA62)+(AA56*AA63);
  double X69 =(X57*X60)+(Y57*X61)+(Z57*X62)+(AA57*X63);
  double Y69 =(X57*Y60)+(Y57*Y61)+(Z57*Y62)+(AA57*Y63);
  double Z69 =(X57*Z60)+(Y57*Z61)+(Z57*Z62)+(AA57*Z63);
  double AA69 =(X57*AA60)+(Y57*AA61)+(Z57*AA62)+(AA57*AA63);
  // CALCULATE J1 ANGLE
  double O13 = atan((AA67)/(AA66));
  if (V9 == 1)
    double P13 = to_degree(O13);
  if (V9 == 2)
    double P13 = to_degree(O13);
  if (V9 == 3)
    double P13 = -180 + to_degree(O13);
  if (V9 == 4)
    double P13 = 180 + to_degree(O13);
  // CALCULATE J2 ANGLE	FWD
 double O18 = sqrt(((abs(AA67))*(abs(AA67)))+((abs(AA66))(abs(AA66))));
 double O19 = AA68-E13;
 double O20 = O18-F13;
 double O21 = sqrt((O19*O19)+(O20*O20));
 double O22 = to_degree(atan(O19/O20));
 double O23 = to_degree(acos(((F14*F14)+(O21*O21)-(abs(E16)*abs(E16)))/(2*F14*O21)));
 double O24 = 180-to_degree(acos(((abs(E16)*abs(E16))+(F14*F14)-(O21*O21))/(2*abs(E16)*F14)));
 double O26 = -(O22+O23);
 double O27 = O24;
 // CALCULATE J2 ANGLE	MID
 double P20 = -O20;
 double P21 = sqrt((O19*O19)+(P20*P20));
 double P22 = to_degree(acos(((F14*F14)+(P21*P21)-(abs(E16)*abs(E16)))/(2*F14*P21)));
 double P23 = to_degree(atan(P20/O19));
 double P24 = 180-to_degree(acos(((abs(E16)*abs(E16))+(F14*F14)-(P21*P21))/(2*abs(E16)*F14)));
 double P25 = 90-(P22+P23);
 double P26 = -180+P25;
 double P27 = P24;
  // J1,J2,J3
  double Q4 = P13;
  if (O20<0)
    double Q5 = P26;
  else
    double Q5 = O26;
  if (O20<0)
    double Q6 = P27;
  else
    double Q6 = O27;
  // J1
 double N36 =cos(to_rad(Q4));
 double O36 =-sin(to_rad(Q4))*cos(D13);
 double P36 =sin(to_rad(Q4))*sin(D13);
 double Q36 =F13*cos(to_rad(Q4));
 double N37 =sin(to_rad(Q4));
 double O37 =cos(to_rad(Q4))*cos(D13);
 double P37 =-cos(to_rad(Q4))*sin(D13);
 double Q37 =F13*sin(to_rad(Q4));
 double N38 = 0;
 double O38 =sin(D13);
 double P38 =cos(D13);
 double Q38 =E13;
 double N39 = 0;
 double O39 = 0;
 double P39 = 0;
 double Q39 = 1;
 // J2
 double N42 =cos(to_rad(Q5));
 double O42 =-sin(to_rad(Q5))*cos(D14);
 double P42 =sin(to_rad(Q5))*sin(D14);
 double Q42 =F14*cos(to_rad(Q5));
 double N43 =sin(to_rad(Q5));
 double O43 =cos(to_rad(Q5))*cos(D14);
 double P43 =-cos(to_rad(Q5))*sin(D14);
 double Q43 =F14*sin(to_rad(Q5));
 double N44 = 0;
 double O44 =sin(D14);
 double P44 =cos(D14);
 double Q44 =E14;
 double N45 = 0;
 double O45 = 0;
 double P45 = 0;
 double Q45 = 1;
 // J3
 double N48 =cos(to_rad((Q6)-90));
 double O48 =-sin(to_rad((Q6)-90))*cos(D15);
 double P48 =sin(to_rad((Q6)-90))*sin(D15);
 double Q48 =F15*cos(to_rad((Q6)-90));
 double N49 =sin(to_rad((Q6)-90));
 double O49 =cos(to_rad((Q6)-90))*cos(D15);
 double P49 =-cos(to_rad((Q6)-90))*sin(D15);
 double Q49 =F15*sin(to_rad((Q6)-90));
 double N50 =0;
 double O50 =sin(D15);
 double P50 =cos(D15);
 double Q50 =E15;
 double N51 =0;
 double O51 =0;
 double P51 =0;
 double Q51 =0;
 // R 0-1
 double S33 =(N30*N36)+(O30*N37)+(P30*N38)+(Q30*N39);
 double T33 =(N30*O36)+(O30*O37)+(P30*O38)+(Q30*O39);
 double U33 =(N30*P36)+(O30*P37)+(P30*P38)+(Q30*P39);
 double V33 =(N30*Q36)+(O30*Q37)+(P30*Q38)+(Q30*Q39);
 double S34 =(N31*N36)+(O31*N37)+(P31*N38)+(Q31*N39);
 double T34 =(N31*O36)+(O31*O37)+(P31*O38)+(Q31*O39);
 double U34 =(N31*P36)+(O31*P37)+(P31*P38)+(Q31*P39);
 double V34 =(N31*Q36)+(O31*Q37)+(P31*Q38)+(Q31*Q39);
 double S35 =(N32*N36)+(O32*N37)+(P32*N38)+(Q32*N39);
 double T35 =(N32*O36)+(O32*O37)+(P32*O38)+(Q32*O39);
 double U35 =(N32*P36)+(O32*P37)+(P32*P38)+(Q32*P39);
 double V35 =(N32*Q36)+(O32*Q37)+(P32*Q38)+(Q32*Q39);
 double S36 =(N33*N36)+(O33*N37)+(P33*N38)+(Q33*N39);
 double T36 =(N33*O36)+(O33*O37)+(P33*O38)+(Q33*O39);
 double U36 =(N33*P36)+(O33*P37)+(P33*P38)+(Q33*P39);
 double V36 =(N33*Q36)+(O33*Q37)+(P33*Q38)+(Q33*Q39);
 // R 0-2
 double S39 =(S33*N42)+(T33*N43)+(U33*N44)+(V33*N45);
 double T39 =(S33*O42)+(T33*O43)+(U33*O44)+(V33*O45);
 double U39 =(S33*P42)+(T33*P43)+(U33*P44)+(V33*P45);
 double V39 =(S33*Q42)+(T33*Q43)+(U33*Q44)+(V33*Q45);
 double S40 =(S34*N42)+(T34*N43)+(U34*N44)+(V34*N45);
 double T40 =(S34*O42)+(T34*O43)+(U34*O44)+(V34*O45);
 double U40 =(S34*P42)+(T34*P43)+(U34*P44)+(V34*P45);
 double V40 =(S34*Q42)+(T34*Q43)+(U34*Q44)+(V34*Q45);
 double S41 =(S35*N42)+(T35*N43)+(U35*N44)+(V35*N45);
 double T41 =(S35*O42)+(T35*O43)+(U35*O44)+(V35*O45);
 double U41 =(S35*P42)+(T35*P43)+(U35*P44)+(V35*P45);
 double V41 =(S35*Q42)+(T35*Q43)+(U35*Q44)+(V35*Q45);
 double S42 =(S36*N42)+(T36*N43)+(U36*N44)+(V36*N45);
 double T42 =(S36*O42)+(T36*O43)+(U36*O44)+(V36*O45);
 double U42 =(S36*P42)+(T36*P43)+(U36*P44)+(V36*P45);
 double V42 =(S36*Q42)+(T36*Q43)+(U36*Q44)+(V36*Q45);
 // R 0-3
 double S45 =(S39*N48)+(T39*N49)+(U39*N50)+(V39*N51);
 double T45 =(S39*O48)+(T39*O49)+(U39*O50)+(V39*O51);
 double U45 =(S39*P48)+(T39*P49)+(U39*P50)+(V39*P51);
 double V45 =(S39*Q48)+(T39*Q49)+(U39*Q50)+(V39*Q51);
 double S46 =(S40*N48)+(T40*N49)+(U40*N50)+(V40*N51);
 double T46 =(S40*O48)+(T40*O49)+(U40*O50)+(V40*O51);
 double U46 =(S40*P48)+(T40*P49)+(U40*P50)+(V40*P51);
 double V46 =(S40*Q48)+(T40*Q49)+(U40*Q50)+(V40*Q51);
 double S47 =(S41*N48)+(T41*N49)+(U41*N50)+(V41*N51);
 double T47 =(S41*O48)+(T41*O49)+(U41*O50)+(V41*O51);
 double U47 =(S41*P48)+(T41*P49)+(U41*P50)+(V41*P51);
 double V47 =(S41*Q48)+(T41*Q49)+(U41*Q50)+(V41*Q51);
 double S48 =(S42*N48)+(T42*N49)+(U42*N50)+(V42*N51);
 double T48 =(S42*O48)+(T42*O49)+(U42*O50)+(V42*O51);
 double U48 =(S42*P48)+(T42*P49)+(U42*P50)+(V42*P51);
 double V48 =(S42*Q48)+(T42*Q49)+(U42*Q50)+(V42*Q51);
 // R 0-3 transposed
 double S51 =S45;
 double T51 =S46;
 double U51 =S47;
 double S52 =T45;
 double T52 =T46;
 double U52 =T47;
 double S53 =U45;
 double T53 =U46;
 double U53 =U47;
 // R 3-6 (spherical wrist  orietation)
 double X72 =(S51*X66)+(T51*X67)+(U51*X68);
 double Y72 =(S51*Y66)+(T51*Y67)+(U51*Y68);
 double Z72 =(S51*Z66)+(T51*Z67)+(U51*Z68);
 double X73 =(S52*X66)+(T52*X67)+(U52*X68);
 double Y73 =(S52*Y66)+(T52*Y67)+(U52*Y68);
 double Z73 =(S52*Z66)+(T52*Z67)+(U52*Z68);
 double X74 =(S53*X66)+(T53*X67)+(U53*X68);
 double Y74 =(S53*Y66)+(T53*Y67)+(U53*Y68);
 double Z74 =(S53*Z66)+(T53*Z67)+(U53*Z68);
 // WRIST ORENTATION
 double R7 = to_degree(atan2(Z73,Z72));
 double R8 = to_degree(atan2(+sqrt(1-(Z74*Z74)),Z74));
  if (Y74 < 0)
    double R9 = to_degree(atan2(-Y74,X74))-180;
  else
    double R9 = to_degree(atan2(-Y74,X74))+180;
  double S7 = to_degree(atan2(-Z73,-Z72));
  double S8 = to_degree(atan2(-sqrt(1-(Z74*Z74)),Z74));
  if (Y74 < 0)
    S9 = to_degree(atan2(Y74,-X74))+180;
  else
    double S9 = to_degree(atan2(Y74,-X74))-180;
  if (V8 == "F")
    double Q8 = R8;
  else
    double Q8 = S8;
  if(Q8>0)
    double Q7 = R7;
  else
    double Q7 = S7;
  if(Q8<0)
    double Q9 = S9;
  else
    double Q9 = R9;

    
    JXdeg  <<  Q4, Q5, Q6, Q7, Q8, Q9;
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
