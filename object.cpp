#include "object.h"

vec4 object::diffuse(vec2 uv)
{
	if(shader.has_diff_map)
	{
		// std::cout << uv.u << ", " << uv.v << std::endl;
		vec4 map_color = shader.get_map_color(uv);
		return map_color*shader.diffuse_color*shader.diffuse;
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