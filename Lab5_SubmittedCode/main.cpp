/*
ENGO431 Lab 3
Claire Mah
Vincent Cung
Jad Shehadeh
Zian Zahid*/
#include <iostream>
#include <iomanip>
#include "LSA.h"

using namespace std;

int main(){

//run program using test data from lecture slides
LSA test;
test.c = 152.15;

cout << fixed << setprecision(10) << "---------------BEGINNING SPACE INTERSECTION FOR TEST DATA-------------" << endl;
test.SpaceInt_data("in_int_test_obs.txt", "in_int_test_EOP.txt", "in_int_test_x0.txt", "in_int_test_c_size_std_Tol.txt");
test.SI_Run_LSA();

MatToText(test.SI_All_X, "out_test_Object_TiePoints.txt");
//MatToText(test.RMSE, "out_test_RMSE.txt");
MatToText(test.SI_All_v, "out_test_Residuals_TiePoints.txt");


cout << "residuals" << endl << test.SI_All_v << endl;

cout << "---------------BEGINNING SPACE INTERSECTION FOR THE 16 TIE POINTS-------------" << endl;

LSA lab5;
lab5.SpaceInt_data2("in_tie_points.txt", "in_output27_x_radians.txt", "in_output28_x_radians.txt", "in_int_info.txt");
lab5.SI_Run_LSA();

MatToText(lab5.SI_All_X, "out_Object_TiePoints.txt");
MatToText(lab5.RMSE, "out_RMSE.txt");
MatToText(lab5.SI_All_v, "out_Residuals_TiePoints.txt");
MatToText(lab5.SI_Std, "Standard_Deviation.txt");
cout << "residuals" << endl << lab5.SI_All_v << endl;

return 0;
};
