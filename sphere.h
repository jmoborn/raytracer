#ifndef __SPHERE_H__
#define __SPHERE_H__

#include <iostream>
#include <string>
#include <math.h>
#include "vec4.h"

class sphere{
public:
	sphere(double radius, vec4 center);
	~sphere(){};
	bool intersect(ray& r);
	vec4 shade();
	double r;
	vec4 c;
};

#endif