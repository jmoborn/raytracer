#ifndef __OBJECT_H__
#define __OBJECT_H__

#include <iostream>
#include <string>
#include "vec4.h"
#include "material.h"

class object{
public:

	virtual vec4 get_normal(const vec4& p)=0;

	/*
	Calculate if a given ray intersects with the object.
	Returns true if this is the case.
	The given ray's t, hit_norm (not guaranteed to be normalized), and hit_color parameters 
		are only modifiedif the new intersection point is less than t and greater than 0.
	*/
	virtual bool intersect(ray& r)=0;
	virtual vec4 diffuse();
	virtual vec4 reflect();
	virtual vec4 refract();

	material shader;
};

#endif