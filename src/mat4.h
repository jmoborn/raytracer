#ifndef __MAT4_H__
#define __MAT4_H__
#define _USE_MATH_DEFINES

#include <math.h>
#include <algorithm>
#include <iostream>

/*

4x4 Matrix class for use in storing
homogeneous transformations.

Author: Jeremy Oborn

*/
class mat4
{
public:
	mat4();
	~mat4(){};

	void translate(double tx, double ty, double tz);
	void scale(double sx, double sy, double sz);
	void rotateX(double theta);
	void rotateY(double theta);
	void rotateZ(double theta);

	double data[4][4];

	double & operator()(const int i, const int j);
	double operator()(const int i, const int j) const;

	mat4 & operator*=(const mat4& m);
	mat4 & operator*=(const double c);

	friend mat4 operator*(const mat4& m1, mat4& m2);
	friend mat4 operator*(const mat4&m, const double c);
};

#endif