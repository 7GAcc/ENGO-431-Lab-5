#pragma once

#include "LSA.h"

#include <vector>
#include <cmath>
#include <Eigen/Dense>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <string>

using namespace Eigen;
using namespace std;

//General reading function
MatrixXd Matrix_readIn(int col, string filename);

class LSA {
    public:

    //Constructor
    LSA();

    int p = 4;          //number of control points
    int n = 2*p;        //number of observations, x and y of GCPs
    int u = 6;          //number of unknowns, xc, yc, zc, omega, phi, kappa
    int N_iter = 0;     //number of iterations
//    double bx = 92;     //base distance 92 mm
    double c;           //calibrated focal length //152.15;  //153.358;
    MatrixXd Con;       //Constants, GCP coordinates, xyz
    MatrixXd l_o;       //observations
    MatrixXd x_o;       //estimate of parameters
    MatrixXd x;         //vector of parameters
    MatrixXd delta;     //adjustment vector
    MatrixXd w;         //misclosure
    MatrixXd A;         //design matrix
    MatrixXd Cx;        //variance covariance of estimate
    MatrixXd v;         //residuals
    double SI_AverageHeight;
    /*---------------------------------------------------------------------------------*/
    int num_of_pairs;                // Number of points per image
    MatrixXd UVW_L, UVW_R;
    MatrixXd SIobs;                  // Space intersection observations 
    MatrixXd SI_EOP_L, SI_EOP_R;     // External orientation paramters
    MatrixXd SI_x0;       // Estimated paramters
	MatrixXd SI_info;                // stores the value of c in mm, image format size in mm assuming x and y are equal, stores the standard deviation in m, and tolerance in m
    MatrixXd SI_P;                   // Weight matrix
    MatrixXd SI_Moi_L, SI_Moi_R;
    int SI_n;                        //
    MatrixXd SI_A;                  //Design matrix for space intersection
    MatrixXd SI_w;                  //Misclosure vector for space intersection
    MatrixXd b;
    MatrixXd SI_v;
    MatrixXd SI_x;
    MatrixXd SI_delta;
    MatrixXd SI_CxHat;
    MatrixXd SI_CxHat_Std;
    MatrixXd SI_R;
    //
    MatrixXd SI_All_X;
    MatrixXd SI_All_v;
    MatrixXd RMSE;
    MatrixXd SI_Std;
    //
    void SI_Collin_Eqs(int point);
    void SI_xHat_Calc(int point);
    void SpaceInt_data(string obs_file, string EOP_file, string x0_file, string otherinfo_file);
    void SpaceInt_data2(string obs_file, string x_L_file, string x_R_file, string info_file);
    void SI_Design(int point);       //Computes edesign matrix for a specific pair of points
    void SI_LSA_Single(int point);   //Runs least square adjustment for a specific point (used in SI_Run_LSA()) and updates respective variables
    void SI_Run_LSA();
    void SI_Approx_x0();
    /*---------------------------------------------------------------------------------*/

//    MatrixXd vL;        //vector left image, xL, yL, -c
//    MatrixXd vR;        //vector right image, xR, yR, -c
//    MatrixXd vR_;       //vR rotated by MR^T
//    MatrixXd B;         //bx, by, bz
//    MatrixXd corr;      //correlation matrix
//    double lambda;      //space intersection parameter
//    double miu;         //space intersection parameter
//    MatrixXd Lm;        //Left model space coordinates
//    MatrixXd Rm;        //Right model space coordinates
//    MatrixXd model;     //model space coordinates
//    MatrixXd pY;        //y parallax
    MatrixXd m;         //m rotation matrix R3kappa*R2phi*R1omega

    /*Function: adjustment
    Runs an iterative least squares adjustment
    Inputs: None
    Outputs: None
    */
    void adjustment();

    /*Function: designA
    Populates design matrix A and misclosure vector w
    Inputs: None
    Outputs: None
    */
    void designA();

    /*Function: data
    Reads text file containing xy observations and text file containing known xyz of control points
    Inputs: string of filename x2
    Outputs: None
    */
    void data(string file, string const_file);
    

    /*Function: initialx
    Populates initial estimates of Xc Yc Zc omega phi kappa
    Inputs: None
    Outputs: None
    */
    void initialx();

    
    

    /*Function: variances
    Populates variance matrix and correlation matrix
    Assumes uniform variances for all observations
    Inputs: None
    Outputs: None
    */
//    void variances();

    /*Function spaceIntersection
    Performs space intersection. Calculates y parallax and model space coordinates
    Inputs: name of input text file, number of points
    Outputs: None
    */
//    void spaceIntersection();

    /*Function spaceIntersection, overload
    use when you want to perform space intersection on more points than were used to compute RO parameters
    Inputs: name of input text file, number of points
    Outputs: None
    */
//    void spaceIntersection(string in, int N_points);

    /*Function: output
    Prints parameters, correlation, parallax, model coordinates to textfiles
    Inputs: string to identify object
    Outputs: None
    */
//    void output(string name);
};


    /*Function: print_mat
    Prints matrix to console
    Inputs: MatrixXd to be printed and string for a label
    Outputs: None
    */
    void print_mat(MatrixXd mat, string name);

    /* Function: MatToText
    * Prints Matrix to text file
    */
    void MatToText(MatrixXd M, string filename);

    /*Function: R
    Returns roation matrix about specified axis by specified angle in radians.
    Inputs: int 1, 2, or 3 for axis of rotation. Angle
    Outputs: rotation matrix
    */
    MatrixXd R(int i, double angle);
