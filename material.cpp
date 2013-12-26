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
}

material::~material(){}

