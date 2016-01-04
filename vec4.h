#ifndef __VEC4_H__
#define __VEC4_H__

#include <math.h>
#include <iostream>
#include <limits>

#include "mat4.h"

/*

Vector class of four components meant for use in storing
points, rays, normals etc. in homogeneous coordinates.

Author: Jeremy Oborn

*/
class vec4 {
public:
	vec4();
	vec4(double x, double y, double z);
	vec4(double x, double y, double z, double w);
	~vec4(){};

	double length() const;
	double dot(vec4& v);
	vec4 cross(vec4& v);
	vec4 reflect(vec4& normal);
	vec4 refract(vec4& normal, double n1, double n2, double &cos2t);
	void normalize();
	void clamp(double lower, double upper);

	double x;
	double y;
	double z;
	double w;

	double operator==(const vec4& v);
	double operator!=(const vec4& v);
	double operator<(const vec4& v);
	double operator<=(const vec4& v);
	double operator>(const vec4& v);
	double operator>=(const vec4& v);

	vec4 & operator=(const vec4& v);
	vec4 & operator+=(const vec4& v);
	vec4 & operator-=(const vec4& v);
	vec4 & operator=(const double c);
	vec4 & operator+=(const double c);
	vec4 & operator-=(const double c);
	vec4 & operator*=(const double c);
	vec4 & operator*=(const mat4& m);

	friend vec4 operator-(const vec4& v1, const vec4& v2);
	friend vec4 operator+(const vec4& v1, const vec4& v2);
	friend vec4 operator-(const vec4& v, const double c);
	friend vec4 operator-(const double c, const vec4& v);
	friend vec4 operator+(const vec4& v, const double c);
	friend vec4 operator+(const double c, const vec4& v);
	friend vec4 operator*(const mat4& m, const vec4& v);
	friend vec4 operator*(const vec4& v, const double c);
	friend vec4 operator*(const double c, const vec4& v);
	friend vec4 operator*(const vec4& v1, const vec4& v2);

	friend std::ostream& operator<<(std::ostream &out, vec4& v1);
};

class vec2 {
public:
	vec2();
	vec2(double u, double v);
	~vec2(){};

	double u;
	double v;

	vec2 & operator*=(const double c);

	friend vec2 operator+(const vec2& v1, const vec2& v2);
	friend vec2 operator*(const vec2& v, const double c);
};

/*
TODO: actual wavelength spectrum representation
*/
class wavelength {
public:
	wavelength(int samples=3, vec4 color=vec4(1,1,1), double nwl=0.5);
	~wavelength(){};
	int get_num_samples();
	vec4 get_color_sample(int sample_num, double random,  double &nwl);
	wavelength get_sample(int sample_num, double random);
	vec4 get_total_color();
	vec4 get_color_from_mult(double mult);
	double get_wavelength();
	double get_range_low();
	double get_range_high();
private:
	int num_samples;
	vec4 total_color;
	double n_wavelength;
	double range_low;
	double range_high;
};

class ray {
public:
	ray(vec4 origin=vec4(0,0,0), vec4 direction=vec4(0,0,0));
	~ray(){};

	ray inherit();
	vec4 end();

	vec4 o;
	vec4 d;
	double t;
	// double refract_obj;
	vec4 hit_norm;
	vec4 hit_color;
	vec2 hit_uv;
	int hit_mtl;
	double prior_ior;
	int debug;
	wavelength spectrum;
	short stack[5];
};

#endif