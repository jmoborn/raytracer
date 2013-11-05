#ifndef __SPHERE_H__
#define __SPHERE_H__

#include <iostream>
#include <string>
#include <math.h>
#include <algorithm>

#include "vec4.h"
#include "object.h"

class sphere : public object{
public:
	sphere(double radius, vec4 center);
	sphere(double radius, vec4 center, vec4 color, double diff = 1.0);
	~sphere(){};
	bool intersect(ray& r);
	vec4 diffuse();
	vec4 reflect();
	vec4 get_normal(const vec4& p);
	double r;
	vec4 c;
	vec4 color;
	double diffuse_mult;
	double reflect_mult;
};

#endif