#ifndef __OBJECT_H__
#define __OBJECT_H__

#include <iostream>
#include <string>
#include "vec4.h"

class object{
public:


	virtual vec4 get_normal(const vec4& p)=0;

	/*
	Calculate if a given ray intersects with the object.
	Returns true if this is the case.
	The given ray's t, hit_norm, and hit_color parameters are only modified
		 if the new intersection point is less than t and greater than 0.
	*/
	virtual bool intersect(ray& r)=0;
	virtual vec4 diffuse()=0;
	virtual vec4 reflect()=0;
};

#endif