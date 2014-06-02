#include "vec4.h"

vec4::vec4()
{
	this->x = 0.0;
	this->y = 0.0;
	this->z = 0.0;
	this->w = 1.0;
}

vec4::vec4(double x, double y, double z)
{
	this->x = x;
	this->y = y;
	this->z = z;
	this->w = 1.0;
}

vec4::vec4(double x, double y, double z, double w)
{
	if(w!=0)
	{
		this->x = x/w;
		this->y = y/w;
		this->z = z/w;
		this->w = 1.0;
	}
	else
	{
		this->x = x;
		this->y = y;
		this->z = z;
		this->w = w;
	}
}

double vec4::length() const
{
	return sqrt(x*x + y*y + z*z);
}

double vec4::dot(vec4& v)
{
	return x*v.x + y*v.y + z*v.z;
}

vec4 vec4::cross(vec4& v)
{
	double xx = y*v.z - v.y*z;
	double xy = z*v.x - v.z*x;
	double xz = x*v.y - v.x*y;
	return vec4(xx, xy, xz);
}

void vec4::normalize()
{
	double l_inv = 1/length();
	x *= l_inv;
	y *= l_inv;
	z *= l_inv;
}

void vec4::clamp(double limit)
{
	if(x>limit)x=limit;
	if(y>limit)y=limit;
	if(z>limit)z=limit;
}

double vec4::operator==(const vec4& v)
{
	return (x==v.x && y==v.y && z==v.z);
}
double vec4::operator!=(const vec4& v)
{
	return (x!=v.x || y!=v.y || z!=v.z);
}
double vec4::operator<(const vec4& v)
{
	return this->length()<v.length();
}
double vec4::operator<=(const vec4& v)
{
	return this->length()<=v.length();
}
double vec4::operator>(const vec4& v)
{
	return this->length()>v.length();
}
double vec4::operator>=(const vec4& v)
{
	return this->length()>=length();
}
vec4 & vec4::operator=(const vec4& v)
{
	x = v.x;
	y = v.y;
	z = v.z;
	w = v.w;
	return *this;
}
vec4 & vec4::operator+=(const vec4& v)
{
	x += v.x;
	y += v.y;
	z += v.z;
	return *this;
}
vec4 & vec4::operator-=(const vec4& v)
{
	x -= v.x;
	y -= v.y;
	z -= v.z;
	return *this;
}
vec4 & vec4::operator=(const double c)
{
	x = y = z = c;
	w = 1.0;
	return *this;
}
vec4 & vec4::operator+=(const double c)
{
	x += c;
	y += c;
	z += c;
	return *this;
}
vec4 & vec4::operator-=(const double c)
{
	x -= c;
	y -= c;
	z -= c;
	return *this;
}
vec4 & vec4::operator*=(const double c)
{
	x *= c;
	y *= c;
	z *= c;
	return *this;
}
vec4 & vec4::operator*=(const mat4& m)
{	
	double mx = m(0,0)*x + m(0,1)*y + m(0,2)*z + m(0,3)*w;
	double my = m(1,0)*x + m(1,1)*y + m(1,2)*z + m(1,3)*w;
	double mz = m(2,0)*x + m(2,1)*y + m(2,2)*z + m(2,3)*w;
	double mw = m(3,0)*x + m(3,1)*y + m(3,2)*z + m(3,3)*w;
	x = mx/mw;
	y = my/mw;
	z = mz/mw;
	w = 1.0;
	return *this;
}

vec4 operator-(const vec4& v1, const vec4& v2)
{
	return vec4(v1.x-v2.x, v1.y-v2.y, v1.z-v2.z);
}
vec4 operator+(const vec4& v1, const vec4& v2)
{
	return vec4(v1.x+v2.x, v1.y+v2.y, v1.z+v2.z);
}
vec4 operator-(const vec4& v, double c)
{
	return vec4(v.x-c, v.y-c, v.z-c);
}
vec4 operator+(const vec4& v, double c)
{
	return vec4(v.x+c, v.y+c, v.z+c);
}
vec4 operator*(const mat4& m, const vec4& v)
{
	double mvx = m(0,0)*v.x + m(0,1)*v.y + m(0,2)*v.z + m(0,3)*v.w;
	double mvy = m(1,0)*v.x + m(1,1)*v.y + m(1,2)*v.z + m(1,3)*v.w;
	double mvz = m(2,0)*v.x + m(2,1)*v.y + m(2,2)*v.z + m(2,3)*v.w;
	double mvw = m(3,0)*v.x + m(3,1)*v.y + m(3,2)*v.z + m(3,3)*v.w;
	return vec4(mvx, mvy, mvz, mvw);
}
vec4 operator*(const vec4& v, const double c)
{
	return vec4(v.x*c, v.y*c, v.z*c);
}
vec4 operator*(const vec4& v1, const vec4& v2)
{
	return vec4(v1.x*v2.x, v1.y*v2.y, v1.z*v2.z);
}

vec2::vec2()
{
	u = 0.0;
	v = 0.0;
}

vec2::vec2(double u, double v)
{
	this->u = u;
	this->v = v;
}

ray::ray(vec4 origin, vec4 direction)
{
	this->o = origin;
	this->d = direction;
	this->t = std::numeric_limits<double>::infinity();
	this->refract_obj = -1;
}

vec4 ray::end()
{
	return this->o + (this->d*this->t);
}