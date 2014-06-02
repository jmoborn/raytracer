#include "object.h"

vec4 object::diffuse()
{
	if(shader.has_diff_map)
	{
		vec4 map_color;
		return map_color*shader.diffuse_color*shader.diffuse
	}
	else
	{
		return shader.diffuse_color*shader.diffuse;
	}
}

vec4 object::reflect()
{
	return shader.reflect_color*shader.reflect;
}

vec4 object::refract()
{
	return shader.diffuse_color*shader.refract;
}