#ifndef __MATERIAL_H__
#define __MATERIAL_H__

#include "vec4.h"

class material {
public:
	material(vec4 c=vec4(0.5,0.5,0.5), vec4 ref_color=vec4(0.5,0.5,0.5), double diff=1.0, 
			 double ref=0.0, double spec=32.0, double fract=0.0, double ior=1.0);
	~material();
	vec4 diffuse_color;
	vec4 reflect_color;
	double diffuse;
	double reflect;
	double specular;
	double refract;
	double ior;

};

#endif