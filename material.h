#ifndef __MATERIAL_H__
#define __MATERIAL_H__

#include "vec4.h"
#include "pixelmap.h"

class material {
public:
	// material(vec4 c=vec4(1.0,1.0,1.0), vec4 ref_color=vec4(0.5,0.5,0.5), double diff=1.0, 
	// 		 double ref=0.0, double spec=32.0, double fract=0.0, double ior=1.0);
	material(vec4 diff_c=vec4(1.0,1.0,1.0), vec4 refl_c=vec4(0.0,0.0,0.0), 
			 vec4 refr_c=vec4(0.0,0.0,0.0), vec4 emit_c=vec4(1.0,1.0,1.0), 
			 double diff=1.0, double refl=0.0, double refr=0.0, double emit=0.0, double ior=1.0);
	~material();
	void add_diffuse_map(pixelmap &p);
	vec4 get_map_color(vec2 uv);
	vec4 get_diffuse_color(vec2 uv=vec2(0,0));
	vec4 get_reflect_color();
	vec4 get_refract_color();
	double get_dispersion_ior(vec4 color);
	vec4 diffuse_color;
	vec4 reflect_color;
	vec4 refract_color;
	vec4 emit_color;
	double diffuse;
	double reflect;
	double refract;
	double emit;
	double ior;
	bool has_diff_map;
	pixelmap diff_map;
	
	double energy;

};

#endif