#ifndef __OBJECT_H__
#define __OBJECT_H__

#include <iostream>
#include <string>
#include "vec4.h"

class object{
public:
	virtual vec4 get_normal(const vec4& p)=0;//{return vec4(0,0,0);}
	virtual bool intersect(ray& r)=0;//{return vec4(0,0,0);};
	virtual vec4 get_color()=0;//{return vec4(0,0,0);};
};

#endif