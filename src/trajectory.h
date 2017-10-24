#ifndef TK_SPLINE_H
#define TK_SPLINE_H

#include <fstream>
#include <math.h>
#include <iostream>
#include <vector>
#include "Eigen-3.3/Eigen/Core"
#include "Eigen-3.3/Eigen/QR"

using namespace std;
using namespace Eigen;

// ------------
// Declarations
// ------------

// Define templates min jerk and snap trajectory
Eigen::VectorXd MJT(Eigen::VectorXd, double);
Eigen::VectorXd MST(Eigen::VectorXd, double);

// --------------
// Implementation
// --------------

// For using Eigen to solve linear system to equation, see following
// https://eigen.tuxfamily.org/dox/group__TutorialLinearAlgebra.html

Eigen::VectorXd MJT(Eigen::VectorXd bounds, double t) {

	/*
	Solve for minimum jerk trajectory trajectory.
	bounds: Boundary conditions at t_i=0 & t_f=t
            [s_i, s_i_dot, s_i_double_dot,
             s_f, s_f_dot, s_f_double_dot]
	t: duration of maneuver
	return value: length 6 array of coefficients in 5th degree polynomial.
	*/	
	// use getXY to transform points from (s,d) to (x,y)

    Eigen::MatrixXd A(6,6);
    Eigen::VectorXd b(6);
    Eigen::VectorXd x(6);

    double t2 = t*t;
    double t3 = t2*t;
    double t4 = t3*t;
    double t5 = t4*t;
	A << 1, 0,   0,    0,     0,     0,
		 0, 1,   0,    0,     0,     0,
		 0, 0,   1,    0,     0,     0,
		 1, t,  t2,   t3,    t4,    t5,
		 0, 1, 2*t, 3*t2,  4*t3,  5*t4,
		 0, 0,   2,  6*t, 12*t2, 20*t3;

	x = A.colPivHouseholderQr().solve(b);

    return x;
	
}


Eigen::VectorXd MST(Eigen::VectorXd bounds, double t) {

    /*
	Solve for minimum snap trajectory trajectory.
	bounds: Boundary conditions at t_i=0 & t_f=t
            [s_i, s_i_dot, s_i_double_dot, s_i_triple_dot,
             s_f, s_f_dot, s_f_double_dot, s_f_triple_dot]
	t: duration of maneuver
	return value: length 6 array of coefficients in 5th degree polynomial.
	*/	
	// use getXY to transform points from (s,d) to (x,y)

	Eigen::MatrixXd A(8,8);
    Eigen::VectorXd b(8);
    Eigen::VectorXd x(8);

    double t2 = t*t;
    double t3 = t2*t;
    double t4 = t3*t;
    double t5 = t4*t;
	double t6 = t5*t;
    double t7 = t6*t;
    A << 1, 0,   0,    0,     0,     0,      0,      0,
		 0, 1,   0,    0,     0,     0,      0,      0,
		 0, 0,   2,    0,     0,     0,      0,      0,
         0, 0,   0,    6,     0,     0,      0,      0,
		 1, t,  t2,   t3,    t4,    t5,     t6,     t7,
		 0, 1, 2*t, 3*t2,  4*t3,  5*t4,   6*t5,   7*t6,
		 0, 0,   2,  6*t, 12*t2, 20*t3;  30*t4,  42*t5;
         0, 0,   0,    6,  24*t, 60*t2; 120*t3, 210*t4;

	x = A.colPivHouseholderQr().solve(b);

    return x;

}

#endif /* TK_SPLINE_H */