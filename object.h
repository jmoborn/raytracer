#ifndef __OBJECT_H__
#define __OBJECT_H__

#include <iostream>
#include <string>
#include "vec4.h"
#include "material.h"

class object{
public:
	/*
	Calculate if a given ray intersects with the object.
	Returns true if this is the case.
	The given ray's t and hit_norm (not guaranteed to be normalized) parameters 
		are only modified if the new intersection point is less than t and greater than 0.
	*/
	virtual bool intersect(ray& r)=0;
	virtual void print_info(std::ostream &out)=0;
	int get_mtl_idx();
protected:
	int mtl_idx;

};

#endif