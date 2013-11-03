#include "sphere.h"

sphere::sphere(double radius, vec4 center)
{
	this->r = radius;
	this->c = center;
	this->color = (0.5, 0.5, 0.5);
}

sphere::sphere(double radius, vec4 center, vec4 color)
{
	this->r = radius;
	this->c = center;
	this->color = color;
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
		v.t = -B/2;
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
			v.t = root;
			return true;
		}
	}
	else
	{
		v.t = root;
		return true;
	}
	
	// if(root >0)
	// {
	// 	v.t = root;
	// 	return true;
	// }
}

vec4 sphere::get_normal(const vec4& p)
{
	return (p - this->c)*(1/this->r);
}

vec4 sphere::get_color()
{
	return color;
}