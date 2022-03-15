#include "LSA.h"

double M_PI = atan(1) * 4;




LSA::LSA()
{

}

/*---------------------------------------------------------------------------------*/
void LSA::data(string file, string const_file)
{
    //observations
    ifstream infile;
    infile.open(file);
    if (infile.fail())
    {
        cout << "input file failed";
    }
    l_o.resize(n, 1);
    for(int i = 0; i<l_o.rows(); ++i)
    {
        for(int j = 0; j<l_o.cols(); ++j)
        {
            infile >> l_o(i,j);
        }
    }
    print_mat(l_o, "l_o");


    //constants, control points
    ifstream in;
    in.open(const_file);
    if (in.fail())
    {
        cout << "input file failed";
    }
    Con.resize(p,3);
    for(int i = 0; i<Con.rows(); ++i)
    {
        for(int j = 0; j<Con.cols(); ++j)
        {
            in >> Con(i,j);
        }
    }
    print_mat(Con, "Constants");

}

/*---------------------------------------------------------------------------------*/

void LSA::initialx()
{
    //Initial values
    //2D similarity tranformation
    //FIX *************
    x_o.resize(u,1);
//    double tx = 6349.55;        //from lab3?
//    double ty = 3964.65;
//    double a = 1;         //no idea what a and b are
//    double b = 0;
//    double lambda = 7.58564;//sqrt(a*a+b*b);//7.58564;
    double Zave = (276.42 + 280.05 + 266.47 + 248.10)/4;
    double tx, ty, tz, a, b, lambda, K;
    SI_AverageHeight = Zave;
    MatrixXd st_A, st_xo, st_w, st_delta, st_x;
    st_A.resize(n,4);
    st_xo.resize(4,1);
    st_xo.setZero();
    st_x.resize(4,1);
    st_x.setZero();
    st_w.resize(n,1);
    st_delta.resize(4,1);
    st_delta << 5, 5, 5, 5;

    while(st_delta(0,0) > 0.00000001 &&st_delta(1,0) > 0.00000001&&st_delta(2,0) > 0.00000001&&st_delta(3,0) > 0.00000001)
    {
        for(int i = 0; i<p; ++i)
        {
            st_A(2*i,0) = l_o(2*i,0);
            st_A(2*i,1) = -l_o(2*i+1,0);
            st_A(2*i,2) = 1;
            st_A(2*i,3) = 0;
            st_A(2*i+1,0) = l_o(2*i+1,0);
            st_A(2*i+1,1) = l_o(2*i,0);
            st_A(2*i+1,2) = 0;
            st_A(2*i+1,3) = 1;

            st_w(2*i,0) = st_xo(0,0)*l_o(2*i,0)-st_xo(1,0)*l_o(2*i+1,0) +st_xo(2,0) - Con(i,0);
            st_w(2*i+1,0) = st_xo(1,0)*l_o(2*i,0)+st_xo(0,0)*l_o(2*i+1,0) +st_xo(3,0) - Con(i,1);
        }
        st_delta = -(st_A.transpose()*st_A).inverse()*st_A.transpose()*st_w;
        st_x = st_xo + st_delta;
        st_xo = st_x;
        print_mat(st_x, "st_x");
    }

    a = st_x(0,0);
    b = st_x(1,0);
    tx = st_x(2,0);
    ty = st_x(3,0);
    K = atan2(b,a);
    lambda = sqrt(a*a+b*b);
    tz = c*lambda+Zave;
    x_o << tx, ty, tz, 0, 0, K;
    //x_o << 1,1,1,1,1,1;
    print_mat(x_o, "x_o");
}

/*---------------------------------------------------------------------------------*/

void LSA::adjustment()
{
    N_iter++;
	cout << "~~~~~~~~ ITERATION #" << N_iter << " ~~~~~~~~~~\n";

	this->designA();

	delta = -(A.transpose()*A).inverse()*A.transpose()*w;
	x = x_o + delta;
	v = A*delta+w;
	print_mat(delta, "delta");
	print_mat(x, "x");
	cout << "\n omega, phi, kappa (degrees)\n" << x_o(3,0)*180/M_PI << "\n"<< x_o(4,0)*180/M_PI << "\n"<< x_o(5,0)*180/M_PI << "\n";
	print_mat(v, "v");

	for (int i = 0; i < delta.rows(); ++i)
    {
        if (this->delta(i, 0) > 0.00001 || this->delta(i, 0) < -0.00001)
        {
            x_o = x;
            adjustment();
        }
    }
}


/*---------------------------------------------------------------------------------*/

void LSA::designA()
{
    A.resize(n,u);
    w.resize(n,1);
    //m.resize(3,3);
    double omega = x_o(3,0);
    double phi = x_o(4,0);
    double kappa = x_o(5,0);
    m = R(3,kappa)*R(2, phi)*R(1,omega);
    MatrixXd OC(3,1);
    OC << x_o(0,0), x_o(1,0), x_o(2,0);    //PC of image


    for (int i = 0; i<p; ++i)
    {
        MatrixXd OP(3,1);
        OP << Con(i,0), Con(i,1), Con(i,2); //GCP coordinates
        MatrixXd OP_diff = OP-OC;
        MatrixXd UVW = m*OP_diff;
        //misclosure
        w(2*i,0) = -0.006 - c*UVW(0,0)/UVW(2,0) - l_o(2*i,0);
        w(2*i+1,0) = 0.006 - c*UVW(1,0)/UVW(2,0) - l_o(2*i+1,0);

        //Parital derivatives for A
        A(2*i,0) = -c/(UVW(2,0)*UVW(2,0))*(m(2,0)*UVW(0,0)-m(0,0)*UVW(2,0));
        A(2*i,1) = -c/(UVW(2,0)*UVW(2,0))*(m(2,1)*UVW(0,0)-m(0,1)*UVW(2,0));
        A(2*i,2) = -c/(UVW(2,0)*UVW(2,0))*(m(2,2)*UVW(0,0)-m(0,2)*UVW(2,0));
        A(2*i,3) = -c/(UVW(2,0)*UVW(2,0))*(OP_diff(1,0)*(UVW(0,0)*m(2,2)-UVW(2,0)*m(0,2))-OP_diff(2,0)*(UVW(0,0)*m(2,1)-UVW(2,0)*m(0,1)));
        A(2*i,4) = -c/(UVW(2,0)*UVW(2,0))*(OP_diff(0,0)*(-UVW(2,0)*sin(phi)*cos(kappa)-UVW(0,0)*cos(phi))+OP_diff(1,0)*(UVW(2,0)*sin(omega)*cos(phi)*cos(kappa)-UVW(0,0)*sin(omega)*sin(phi))+OP_diff(2,0)*(-UVW(2,0)*cos(omega)*cos(phi)*cos(kappa)+UVW(0,0)*cos(omega)*sin(phi)));
        A(2*i,5) = -c*UVW(1,0)/UVW(2,0);

        A(2*i+1,0) = -c/(UVW(2,0)*UVW(2,0))*(m(2,0)*UVW(1,0)-m(1,0)*UVW(2,0));
        A(2*i+1,1) = -c/(UVW(2,0)*UVW(2,0))*(m(2,1)*UVW(1,0)-m(1,1)*UVW(2,0));
        A(2*i+1,2) = -c/(UVW(2,0)*UVW(2,0))*(m(2,2)*UVW(1,0)-m(1,2)*UVW(2,0));
        A(2*i+1,3) = -c/(UVW(2,0)*UVW(2,0))*(OP_diff(1,0)*(UVW(1,0)*m(2,2)-UVW(2,0)*m(1,2))-OP_diff(2,0)*(UVW(1,0)*m(2,1)-UVW(2,0)*m(1,1)));
        A(2*i+1,4) = -c/(UVW(2,0)*UVW(2,0))*(OP_diff(0,0)*(UVW(2,0)*sin(phi)*sin(kappa)-UVW(1,0)*cos(phi))+OP_diff(1,0)*(-UVW(2,0)*sin(omega)*cos(phi)*sin(kappa)-UVW(1,0)*sin(omega)*sin(phi))+OP_diff(2,0)*(UVW(2,0)*cos(omega)*cos(phi)*sin(kappa)+UVW(1,0)*cos(omega)*sin(phi)));
        A(2*i+1,5) = c*UVW(0,0)/UVW(2,0);

    }
    print_mat(w, "w");
    print_mat(A, "A");
}

/*---------------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------------*/

void LSA::SpaceInt_data2(string obs_file, string x_L_file, string x_R_file, string info_file) {
    SIobs = Matrix_readIn(2, obs_file);
    num_of_pairs = SIobs.rows() / 2;
    //EOPs must be set in text file
    SI_EOP_L.resize(6, 1);
    SI_EOP_R.resize(6, 1);
    SI_EOP_L = Matrix_readIn(1, x_L_file);
    SI_EOP_R = Matrix_readIn(1, x_R_file);
    //SI_EOP_L.resize(6, 1);
    //SI_EOP_R.resize(6, 1);
    //SI_EOP_L= ; //EOPs for left image (angles in radians)
    //SI_EOP_R= ; //EOPs for right image (angles in radians)

    //Getting the M matrices from EOPS (EOPS NEEDED)
    SI_Moi_L = R(3, SI_EOP_L(5, 0)) * R(2, SI_EOP_L(4, 0)) * R(1, SI_EOP_L(3, 0));
    SI_Moi_R = R(3, SI_EOP_R(5, 0)) * R(2, SI_EOP_R(4, 0)) * R(1, SI_EOP_R(3, 0));
    
    SI_info.resize(3, 1);
    SI_info = Matrix_readIn(1, info_file);
    SI_info.conservativeResize(4, 1);
    //4x1 matrix containing  c in mm, average height, standard deviation in m, tolerance in m, respectively
    c = SI_info(0, 0);
    SI_AverageHeight = SI_info(1, 0);
    SI_info(3,0) = SI_info(1,0) / (SI_info(0,0)/1000) / 10 * SI_info(2,0);
    SI_Approx_x0();

    SI_P.resize(4, 4);
    SI_P.fill(0);
    for (int i = 0; i < SI_P.rows(); i++) {
        SI_P(i, i) = 1;
    }
    SI_P = 1 / pow(SI_info(2), 2) * SI_P;

}

void LSA::SI_Approx_x0() {
    SI_x0.resize(num_of_pairs, 3);
    VectorXd Current_x0(3);
    MatrixXd a;
    for (int i = 0; i < num_of_pairs; i++) {
        //Calculating x0 for each point
        a.resize(4, 3);
        b.resize(4, 1);
        double xij_L = SIobs(i, 0);
        double yij_L = SIobs(i, 1);
        double xij_R = SIobs(i + num_of_pairs, 0);
        double yij_R = SIobs(i + num_of_pairs, 1);

        //Making approximate a matrix
        for (int j = 0; j < 3; j++) {
            //for a11,a12,a13
            a(0, j) = xij_L * SI_Moi_L(2, j) + SI_info(0,0) * SI_Moi_L(0, j);
            //for a21,a22,a23
            a(1, j) = yij_L * SI_Moi_L(2, j) + SI_info(0,0) * SI_Moi_L(1, j);
            //for a31,a32,a33
            a(2, j) = xij_R * SI_Moi_R(2, j) + SI_info(0,0) * SI_Moi_R(0, j);
            //for a41,a42,a43
            a(3, j) = yij_R * SI_Moi_R(2, j) + SI_info(0,0) * SI_Moi_R(1, j);
        }
        //making approximate b matrix
        b(0, 0) = a(0, 0) * SI_EOP_L(0, 0) + a(0, 1) * SI_EOP_L(1, 0) + a(0, 2) * SI_EOP_L(2, 0);
        b(1, 0) = a(1, 0) * SI_EOP_L(0, 0) + a(1, 1) * SI_EOP_L(1, 0) + a(1, 2) * SI_EOP_L(2, 0);
        b(2, 0) = a(2, 0) * SI_EOP_R(0, 0) + a(2, 1) * SI_EOP_R(1, 0) + a(2, 2) * SI_EOP_R(2, 0);
        b(3, 0) = a(3, 0) * SI_EOP_R(0, 0) + a(3, 1) * SI_EOP_R(1, 0) + a(3, 2) * SI_EOP_R(2, 0);
        //Getting x0 for point pair
        Current_x0 = (a.transpose() * a).inverse() * a.transpose() * b;
        SI_x0.row(i) = Current_x0;
    }

}


void LSA::SpaceInt_data(string obs_file, string EOP_file, string x0_file, string otherinfo_file) {

    SIobs = Matrix_readIn(2, obs_file);
    num_of_pairs = SIobs.rows() / 2;
    MatrixXd EOPs = Matrix_readIn(2, EOP_file);
    SI_EOP_L = EOPs.col(0);
    SI_EOP_R = EOPs.col(1);
    for (int i = 3; i < 6; i++) {
        SI_EOP_L(i, 0) *= M_PI / 180;
        SI_EOP_R(i, 0) *= M_PI / 180;
    }
    SI_x0 = Matrix_readIn(3, x0_file);
    
    SI_info = Matrix_readIn(1, otherinfo_file);
    //std::cout << "Obs" << endl << SIobs << endl;
    //std::cout << "EOPs Left" << endl << SI_EOP_L << endl;
    //cout << "EOPs Right" << endl << SI_EOP_R << endl;
    //cout << "x0: " << endl << SI_x0 << endl;
    //cout << "SI info" << endl << SI_info << endl;
    SI_P.resize(4,4);
    SI_P.fill(0);
    for (int i = 0; i < SI_P.rows(); i++) {
        SI_P(i, i) = 1;
    }
    SI_P = 1 / pow(SI_info(2), 2) * SI_P;
    SI_Moi_L = R(3, SI_EOP_L(5, 0)) * R(2, SI_EOP_L(4, 0)) * R(1, SI_EOP_L(3, 0));
    SI_Moi_R = R(3, SI_EOP_R(5, 0)) * R(2, SI_EOP_R(4, 0)) * R(1, SI_EOP_R(3, 0));
    //std::cout << "Weight" << endl << SI_P << endl;
}

/*---------------------------------------------------------------------------------*/

void LSA::SI_Run_LSA() {
    
    SI_n = 2 * 2;
    SI_Std.resize(num_of_pairs, 3);
    SI_All_X.resize(num_of_pairs,3);
    SI_All_v.resize(num_of_pairs*2, 2);
    for (int i = 0; i < num_of_pairs; i++) {
        SI_LSA_Single(i);
    }
   
    cout << "All object coordinates: " << endl << SI_All_X << endl;
    cout << "Mean Height of surface: " << endl << SI_All_X.col(2).sum() / SI_All_X.rows() << endl;
    RMSE.resize(2, 2);
    RMSE.fill(0);
    double t1, t2, t3, t4;
    MatrixXd VL(num_of_pairs, 2);
    MatrixXd VR(num_of_pairs, 2);
    for (int i = 0; i < num_of_pairs; i++) {
        VL(i, 0) = SI_All_v(i, 0);
        VL(i, 1) = SI_All_v(i, 1);
        VR(i, 0) = SI_All_v(i + num_of_pairs, 0);
        VR(i, 1) = SI_All_v(i + num_of_pairs, 1);
    }
    t1 = VL.col(0).sum() / VL.rows();
    t2 = VL.col(1).sum() / VL.rows();
    t3 = VR.col(0).sum() / VR.rows();
    t4 = VR.col(1).sum() / VR.rows();
    for (int i = 0; i < VL.rows(); i++) {
        RMSE(0, 0) += pow(VL(i, 0) - t1, 2);
        RMSE(0, 1) += pow(VL(i, 1) - t2, 2);
        RMSE(1, 0) += pow(VR(i, 0) - t3, 2);
        RMSE(1, 1) += pow(VR(i, 1) - t4, 2);
    }


    RMSE(0, 0) = sqrt(RMSE(0, 0) / VL.rows());
    RMSE(0, 1) = sqrt(RMSE(0, 1) / VL.rows());
    RMSE(1, 0) = sqrt(RMSE(1, 0) / VL.rows());
    RMSE(1, 1) = sqrt(RMSE(1, 1) / VL.rows());
    //RMSE is in micrometers
    RMSE *= pow(10, 3);

}

/*---------------------------------------------------------------------------------*/

void LSA::SI_LSA_Single(int point) {
    cout << "-----For point " << point + 1 << " -----" << endl;
    do {
        SI_Design(point);
        SI_Collin_Eqs(point);
        SI_xHat_Calc(point);
        SI_x0.row(point) = SI_x.col(0);

    } while (SI_delta(0, 0) >= SI_info(3) || SI_delta(1, 0) >= SI_info(3) || SI_delta(2, 0) >= SI_info(3));
    //Storing parameters
    SI_All_X.row(point) = SI_x.col(0);
    //Storing all residuals
    SI_All_v(point, 0) = SI_v(0, 0);
    SI_All_v(point, 1) = SI_v(1, 0);
    SI_All_v(point+num_of_pairs, 0) = SI_v(2, 0);
    SI_All_v(point+num_of_pairs, 1) = SI_v(3, 0);
}

/*---------------------------------------------------------------------------------*/

void LSA::SI_Design(int point) {
    SI_A.resize(4, 3);
    VectorXd OPi = SI_x0.row(point);
    VectorXd OCj_L(3), OCj_R(3);
    OCj_L << SI_EOP_L(0, 0), SI_EOP_L(1, 0), SI_EOP_L(2, 0);
    OCj_R << SI_EOP_R(0, 0), SI_EOP_R(1, 0), SI_EOP_R(2, 0);
    UVW_L = SI_Moi_L * (OPi - OCj_L);
    UVW_R = SI_Moi_R * (OPi - OCj_R);
    double c_W2_L = c / pow(UVW_L(2, 0), 2);
    double c_W2_R = c / pow(UVW_R(2, 0), 2);
    b.resize(4, 1);
    for (int i = 0; i < 3; i++) {
        SI_A(0, i) = c_W2_L * (SI_Moi_L(2, i) * UVW_L(0, 0) - SI_Moi_L(0, i) * UVW_L(2, 0));
        SI_A(1, i) = c_W2_L * (SI_Moi_L(2, i) * UVW_L(1, 0) - SI_Moi_L(1, i) * UVW_L(2, 0));
        SI_A(2, i) = c_W2_R * (SI_Moi_R(2, i) * UVW_R(0, 0) - SI_Moi_R(0, i) * UVW_R(2, 0));
        SI_A(3, i) = c_W2_R * (SI_Moi_R(2, i) * UVW_R(1, 0) - SI_Moi_R(1, i) * UVW_R(2, 0));
    }
    b(0, 0) = SI_A(0, 0) * OCj_L(0) + SI_A(0, 1) * OCj_L(1) + SI_A(0, 2) * OCj_L(2);
    b(1, 0) = SI_A(1, 0) * OCj_L(0) + SI_A(1, 1) * OCj_L(1) + SI_A(1, 2) * OCj_L(2);
    b(2, 0) = SI_A(2, 0) * OCj_R(0) + SI_A(2, 1) * OCj_R(1) + SI_A(2, 2) * OCj_R(2);
    b(3, 0) = SI_A(3, 0) * OCj_R(0) + SI_A(3, 1) * OCj_R(1) + SI_A(3, 2) * OCj_R(2);
    cout << "Design Matrix" << endl << SI_A << endl;
}

/*---------------------------------------------------------------------------------*/

void LSA::SI_Collin_Eqs(int point) {
    SI_w.resize(4, 1);
    double xp = -0.006;
    double yp = 0.006;
    SI_w(0, 0) = (-c * UVW_L(0, 0) / UVW_L(2, 0)) - SIobs(point, 0);
    SI_w(1, 0) = (-c * UVW_L(1, 0) / UVW_L(2, 0)) - SIobs(point, 1);
    SI_w(2, 0) = (-c * UVW_R(0, 0) / UVW_R(2, 0)) - SIobs(point + num_of_pairs, 0);
    SI_w(3, 0) = (-c * UVW_R(1, 0) / UVW_R(2, 0)) - SIobs(point + num_of_pairs, 1);
    cout << "Misclosure: " << endl << SI_w << endl;
}

/*---------------------------------------------------------------------------------*/

void LSA::SI_xHat_Calc(int point) {
    //SI_x = (SI_A.transpose() * SI_A).inverse() * SI_A.transpose() * b;
    SI_x.resize(3, 1);
    SI_delta = -(SI_A.transpose() * SI_P * SI_A).inverse() * SI_A.transpose() * SI_P * SI_w;
    SI_x = SI_x0.row(point).transpose() + SI_delta;
    SI_CxHat = (SI_A.transpose() * SI_P * SI_A).inverse();
    SI_v = SI_A * SI_delta + SI_w;
    SI_CxHat_Std.resize(3, 1);
    for (int i = 0; i < SI_CxHat.rows(); i++) {
        SI_CxHat_Std(i, 0) = sqrt(SI_CxHat(i, i));
    }
    MatrixXd SI_Cv = SI_P.inverse() - SI_A * SI_CxHat * SI_A.transpose();
    SI_R = SI_Cv * SI_P;

    SI_Std.row(point) = SI_CxHat_Std.col(0);
    
    cout << "xHat" << endl << SI_x << endl;
    cout << "Delta: " << endl << SI_delta << endl;
    cout << "Residuals: " << endl << SI_v << endl;
    cout << "RMSE" << endl << RMSE << endl;
    cout << "Standard deviations: " << endl << SI_CxHat_Std << endl;
    cout << "Redundancy: " << endl << SI_R.diagonal() << endl;
    cout << "Sum of Redundancies= " << SI_R.diagonal().sum() << endl;
}
/*---------------------------------------------------------------------------------*/

MatrixXd Matrix_readIn(int col, string filename) {
    ifstream infile;
    MatrixXd result(1, col);
    infile.open(filename, ifstream::in);
    int rows = 1;
    while (!infile.eof()) {
        double a;
        result.conservativeResize(rows, col);
        for (int i = 0; i < col; i++) {
            infile >> a;
            result(rows - 1, i) = a;
        }
        rows++;
    }
    infile.close();
    return result;
}


/*---------------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------------*/


void print_mat(MatrixXd mat, string name)
{
    cout << "\n" << name << "\n";
    for(int i = 0; i<mat.rows(); ++i)
    {
        for(int j = 0; j<mat.cols(); ++j)
        {
            cout << mat(i,j) << "  ";
        }
        cout << "\n";
    }
}

/*---------------------------------------------------------------------------------*/

MatrixXd R(int i, double angle)
{
    MatrixXd mat;
    mat.resize(3,3);
    if(i == 1)
    {
        mat << 1, 0, 0,
            0, cos(angle), sin(angle),
            0, -sin(angle), cos(angle);
    }
    else if(i == 2)
    {
        mat << cos(angle), 0, -sin(angle),
            0, 1, 0,
            sin(angle), 0, cos(angle);
    }
    else if(i == 3)
    {
        mat << cos(angle), sin(angle), 0,
            -sin(angle), cos(angle), 0,
            0, 0, 1;
    }
    else
    {
        cout << "Invalid R matrix";
    }
    return mat;
}

void MatToText(MatrixXd M, string filename) {
    ofstream outfile;
    outfile.open(filename, ofstream::out); //Output file with specified name created
    double val;
    for (int i = 0; i < M.rows(); i++) {
        for (int j = 0; j < M.cols(); j++) {
            val = M(i, j);
            outfile << val << "\t"; //Stores i-th row and j-th column in text file
        }
        outfile << endl;
    }
    outfile.close();
}