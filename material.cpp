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
}

material::~material(){}

material::add_diffuse_map(pixelmap &p)
{
	this->has_diff_map = true;
	this->diff_map = p;
}
