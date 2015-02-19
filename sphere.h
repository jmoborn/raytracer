#ifndef __SPHERE_H__
#define __SPHERE_H__

#include <iostream>
#include <string>
#include <math.h>
#include <algorithm>

#include "vec4.h"
#include "material.h"
#include "object.h"

class sphere : public object{
public:
	sphere(double radius, vec4 center, material m=material());
	~sphere(){};
	bool intersect(ray& r);
	vec4 get_normal(const vec4& p);
    vec2 get_uv(const vec4& p);
	double r;
	vec4 c;
};

#endif