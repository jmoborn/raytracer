#include "object.h"

vec4 object::diffuse()
{
	return shader.diffuse_color*shader.diffuse;
}

vec4 object::reflect()
{
	return shader.reflect_color*shader.reflect;
}

vec4 object::refract()
{
	return shader.diffuse_color*shader.refract;
}