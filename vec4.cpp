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

vec4 vec4::reflect(vec4& normal)
{
	return normal*(dot(normal)*2) - *this;
}

vec4 vec4::refract(vec4& normal, double n1, double n2, double &cos2t)
{
	double eta = n1 / n2;
	double c1 = -dot(normal);
	cos2t = 1 - eta * eta * (1 - c1 * c1);
	if (cos2t<0)
		return vec4(0,0,0);//reflect(normal);
	double c2 = eta*c1-sqrt(cos2t);
	vec4 dir = eta*(*this) + c2*normal;
	return dir;


	// bool into = n.dot(nl)>0;                // Ray from outside going in? 
 //   double nc=1, nt=1.5, nnt=into?nc/nt:nt/nc, ddn=r.d.dot(nl), cos2t; 
 //   if ((cos2t=1-nnt*nnt*(1-ddn*ddn))<0)    // Total internal reflection 
 //     return obj.e + f.mult(radiance(reflRay,depth,Xi)); 
 //   Vec tdir = (r.d*nnt - n*((into?1:-1)*(ddn*nnt+sqrt(cos2t)))).norm(); 
 //   double a=nt-nc, b=nt+nc, R0=a*a/(b*b), c = 1-(into?-ddn:tdir.dot(n)); 
 //   double Re=R0+(1-R0)*c*c*c*c*c,Tr=1-Re,P=.25+.5*Re,RP=Re/P,TP=Tr/(1-P); 
 //   return obj.e + f.mult(depth>2 ? (erand48(Xi)<P ?   // Russian roulette 
 //     radiance(reflRay,depth,Xi)*RP:radiance(Ray(x,tdir),depth,Xi)*TP) : 
 //     radiance(reflRay,depth,Xi)*Re+radiance(Ray(x,tdir),depth,Xi)*Tr); 
}

void vec4::normalize()
{
	double l_inv = 1/length();
	x *= l_inv;
	y *= l_inv;
	z *= l_inv;
}

void vec4::clamp(double lower, double upper)
{
	if(x>upper)x=upper;
	if(y>upper)y=upper;
	if(z>upper)z=upper;
	if(x<lower)x=lower;
	if(y<lower)y=lower;
	if(z<lower)z=lower;
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
vec4 operator-(double c, const vec4& v)
{
	return vec4(v.x-c, v.y-c, v.z-c);
}
vec4 operator+(const vec4& v, double c)
{
	return vec4(v.x+c, v.y+c, v.z+c);
}
vec4 operator+(double c, const vec4& v)
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
vec4 operator*(const double c, const vec4& v)
{
	return vec4(v.x*c, v.y*c, v.z*c);
}
vec4 operator*(const vec4& v1, const vec4& v2)
{
	return vec4(v1.x*v2.x, v1.y*v2.y, v1.z*v2.z);
}
std::ostream& operator<<(std::ostream &out, vec4& v1)
{
	out << "(" << v1.x << ", " << v1.y << ", " << v1.z << ")";
	return out;
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

vec2 & vec2::operator*=(const double c)
{
	u *= c;
	v *= c;
	return *this;
}

vec2 operator+(const vec2& v1, const vec2& v2)
{
	return vec2(v1.u+v2.u, v1.v+v2.v);
}
vec2 operator*(const vec2& v, const double c)
{
	return vec2(v.u*c, v.v*c);
}

wavelength::wavelength(int samples, vec4 color, double nwl)
{
	num_samples = samples;
	// n_wavelength = mult;
	total_color = color;
	if(samples>1)
	{
		range_low = 0.0;
		range_high = 1.0;
	}
	else
	{
		range_low = nwl;
		range_high = nwl;
	}
	n_wavelength = nwl;
	// if(samples==3)
	// {
	// 	total_color = vec4(1,1,1);
	// 	range_low = 0.0;
	// 	range_high = 1.0;
	// }
	// else
	// {
	// 	// total_color = get_color_from_mult(mult);
	// 	range_low = 0.0;
	// 	range_hight = mult;
	// }
}

int wavelength::get_num_samples()
{
	return num_samples;
}

vec4 wavelength::get_color_sample(int sample_num, double random, double &nwl)
{
	if(sample_num<num_samples)
	{
		if(num_samples>1)
		{
			// if(sample_num==0) return vec4(1,0,0);
			// else if(sample_num==1) return vec4(0,1,0);
			// else return vec4(0,0,1);

			double interval = (range_high - range_low)/(double)num_samples;
			// double interval = 0.333333333333;
			nwl = interval*(sample_num + random);
			// if(nwl<0.33333333)nwl=0.0;
			// else if(nwl<0.66666666)nwl=0.5;
			// else nwl=1.0;
			return get_color_from_mult(nwl);
		}
		else
		{
			nwl = n_wavelength;
			return total_color;
		}
	}
	return vec4(0,0,0);
}

wavelength wavelength::get_sample(int sample_num, double random)
{
	if(sample_num<num_samples)
	{
		if(num_samples>1)
		{
			double interval = (range_high - range_low)/(double)num_samples;
			// double interval = 0.333333333333;
			double nwl = interval*(sample_num + random);
			vec4 sample_color = get_color_from_mult(nwl);
			return wavelength(1, sample_color, nwl);
		}
		else
		{
			return wavelength(num_samples, total_color, n_wavelength);
		}
	}
	else
	{
		return wavelength();
	}
}

vec4 wavelength::get_total_color()
{
	return total_color;
}

vec4 wavelength::get_color_from_mult(double mult)
{
	// if(mult<=0.3333333333) return vec4(0,0,1);
	// else if(mult<=0.6666666666) return vec4(0,1,0);
	// else return vec4(1,0,0);

	double b = mult<0.5 ? -6.0*mult+3.0 :  6.0*mult-5.0;
	double g = mult<0.5 ?  6.0*mult-1.0 : -6.0*mult+5.0;
	double r = mult<0.5 ? -6.0*mult+1.0 :  6.0*mult-3.0;
	vec4 color(r, g, b);
	color.clamp(0.0, 1.0);
	// color.normalize();
	// std::cout << "mult: " << mult << " " << color.x << " " << color.y << " " << color.z << std::endl;
	return color;

	// double tmp_1 = 1.0;
	// double tmp_2 = 0.0;
	// double tmp_rgb[3];
	// tmp_rgb[0] = mult+0.333;
	// tmp_rgb[1] = mult;
	// tmp_rgb[2] = mult-0.333;
	// double rgb[3];
	// for(int i=0; i<3; i++)
	// {
	// 	if(tmp_rgb[i]<0)tmp_rgb[i]+=1.0;
	// 	if(tmp_rgb[i]>1)tmp_rgb[i]-=1.0;

	// 	if(6*tmp_rgb[i] < 1) rgb[i] = tmp_2+(tmp_1-tmp_2)*6*tmp_rgb[i];
	// 	else if(2*tmp_rgb[i]<1) rgb[i] = tmp_1;
	// 	else if(3*tmp_rgb[i]<2) rgb[i] = tmp_2+(tmp_1-tmp_2)*(0.666-tmp_rgb[i])*6;
	// 	else rgb[i] = tmp_2; 
	// }
	// return vec4(rgb[0], rgb[1], rgb[2]);
}

/*
return wavelength in micrometers
*/
double wavelength::get_wavelength()
{
	return n_wavelength;

	// if(num_samples!=3)
	// {
	// 	if(total_color.x>0.999)
	// 	{
	// 		return .625+.075*mult;
	// 	}
	// 	if(total_color.y>0.999)
	// 	{
	// 		return .520+.045*mult;
	// 	}
	// 	if(total_color.y>0.999)
	// 	{
	// 		return .440+.060*mult;
	// 	}
	// }
	// else
	// {
	// 	return .400+.300*mult;
	// }
}

double wavelength::get_range_low()
{
	return range_low;
}

double wavelength::get_range_high()
{
	return range_high;
}

ray::ray(vec4 origin, vec4 direction)
{
	this->o = origin;
	this->d = direction;
	this->t = std::numeric_limits<double>::infinity();
	this->prior_ior = 1;
	// this->refract_bounces = 0;
	this->debug = 0;
	for(int i=0; i<5; i++)
	{
		this->stack[i] = -1;
	}
}

vec4 ray::end()
{
	return this->o + (this->d*this->t);
}

ray ray::inherit()
{
	ray result = ray(vec4(0,0,0), vec4(0,1,0));
	result.prior_ior = this->prior_ior;
	result.spectrum = this->spectrum;
	result.debug = this->debug;
	for(int i=0; i<5; i++)
	{
		result.stack[i] = this->stack[i];
	}
	return result;
}
