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
	sphere(double radius, vec4 center, int material_index);
	~sphere(){};
	bool intersect(ray& r);
	void print_info(std::ostream &out);
	vec4 get_normal(const vec4& p);
    vec2 get_uv(const vec4& p);
	double r;
	vec4 c;
};

#endif