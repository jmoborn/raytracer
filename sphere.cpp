#include "sphere.h"

sphere::sphere(double radius, vec4 center)
{
	this->r = radius;
	this->c = center;
}
/*
ray v is a paramaterized equation of a line with origin o and normalized direction d
v: v.o + t*v.d
only positive values of t are considered hits
*/
bool sphere::intersect(ray& v)
{
	double B = 2*(v.d.x*v.o.x - v.d.x*c.x + v.d.y*v.o.y - 
				  v.d.y*c.y + v.d.z*v.o.z - v.d.z*c.z);
	double C = v.o.x*v.o.x - 2*v.o.x*c.x + c.x*c.x + v.o.y*v.o.y - 
				2*v.o.y*c.y + c.y*c.y + v.o.z*v.o.z - 2*v.o.z*c.z + c.z*c.z - r*r;
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
	if(root > 0)
	{
		v.t = root;
		return true;
	}
	root = (-B + sqrt(discriminant))/2;
	if(root >0)
	{
		v.t = root;
		return true;
	}
	return false;
}

vec4 sphere::shade()
{
	return vec4(0,0,0);
}