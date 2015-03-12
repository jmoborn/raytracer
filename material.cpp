#include "material.h"

material::material(vec4 c, vec4 ref_color, double diff, double ref, double spec, double fract, double ior)
{
	this->diffuse_color = c;
	this->reflect_color = ref_color;
	this->diffuse = diff;
	this->reflect = ref;
	this->specular = spec;
	this->refract = fract;
	this->ior = ior;
	this->has_diff_map = false;
	this->emit = 0.0;

	double energy = diffuse+reflect+refract;
	if(energy!=1.0)
	{
		diffuse /= energy;
		reflect /= energy;
		refract /= energy;
	}
}

material::~material(){}

void material::add_diffuse_map(pixelmap &p)
{
	this->has_diff_map = true;
	this->diff_map = p;
}

vec4 material::get_map_color(vec2 uv)
{
	return this->diff_map.getpixel(uv.u, uv.v);
}
