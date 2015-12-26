#include "sphere.h"

sphere::sphere(double radius, vec4 center, int material_index)
{
	this->r = radius;
	this->c = center;
	// this->shader = m;
	this->mtl_idx = material_index;
}

/*
ray v is a paramaterized equation of a line with origin o and normalized direction d
v: v.o + t*v.d
only positive values of t are considered hits
*/
bool sphere::intersect(ray& v)
{
	// double B = 2*(v.d.x*v.o.x - v.d.x*c.x + v.d.y*v.o.y - 
	// 			  v.d.y*c.y + v.d.z*v.o.z - v.d.z*c.z);
	// double C = v.o.x*v.o.x - 2*v.o.x*c.x + c.x*c.x + v.o.y*v.o.y - 
	// 			2*v.o.y*c.y + c.y*c.y + v.o.z*v.o.z - 2*v.o.z*c.z + c.z*c.z - r*r;
	vec4 diff = v.o - c;
	double B = (diff.dot(v.d))*2;
	double C = diff.dot(diff) - r*r;
	double discriminant = B*B - 4*C;
	if(discriminant<0)
	{
		return false;
	}
	if (discriminant==0)
	{
		if((-B/2)<v.t)
		{
			v.t = -B/2;
			vec4 end = v.end();
			v.hit_norm = get_normal(end);
			v.hit_uv = get_uv(end);
			v.hit_mtl = mtl_idx;
			// v.hit_color = diffuse(v.hit_uv);
		}
		return true;
	}
	double root = (-B - sqrt(discriminant))/2;
	if(root <= 0)
	{
		root = (-B + sqrt(discriminant))/2;
		if(root <= 0)
		{
			return false;
		}
		else
		{
			if(root < v.t)
			{
				v.t = root;
				vec4 end = v.end();
				v.hit_norm = get_normal(end);
				v.hit_uv = get_uv(end);
				v.hit_mtl = mtl_idx;
				// v.hit_color = diffuse(v.hit_uv);
			}
			return true;
		}
	}
	else
	{
		if(root < v.t)
		{
			v.t = root;
			vec4 end = v.end();
			v.hit_norm = get_normal(end);
			v.hit_uv = get_uv(end);
			v.hit_mtl = mtl_idx;
			// v.hit_color = diffuse(v.hit_uv);
		}
		return true;
	}
}

vec4 sphere::get_normal(const vec4& p)
{
	return (p - this->c)*(1/this->r);
}

vec2 sphere::get_uv(const vec4& p)
{
	vec4 p_rel = p - c;
	double u = acos(p_rel.z/(sqrt(p_rel.x*p_rel.x+p_rel.y*p_rel.y+p_rel.z*p_rel.z)))/(M_PI);
	double v = (atan(p_rel.y/p_rel.z)+M_PI)/(2.0*M_PI);
	return vec2(u, v);
// std::cout << "u: " << u << " v: " << v << std::endl;
}