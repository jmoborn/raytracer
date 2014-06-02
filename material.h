#ifndef __MATERIAL_H__
#define __MATERIAL_H__

#include "vec4.h"
#include "pixelmap.h"

class material {
public:
	material(vec4 c=vec4(0.5,0.5,0.5), vec4 ref_color=vec4(0.5,0.5,0.5), double diff=1.0, 
			 double ref=0.0, double spec=32.0, double fract=0.0, double ior=1.0);
	~material();
	add_diffuse_map(pixelmap &p);
	vec4 diffuse_color;
	vec4 reflect_color;
	double diffuse;
	bool has_diff_map;
	pixelmap diff_map;
	double reflect;
	double specular;
	double refract;
	double ior;

};

#endif