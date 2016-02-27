#include "material.h"

material::material(vec4 diff_c, vec4 refl_c, vec4 refr_c, vec4 emit_c, double diff, double refl, double refr, double emit, double ior)
{
	this->diffuse_color = diff_c;
	this->reflect_color = refl_c;
	this->refract_color = refr_c;
	this->emit_color = emit_c;
	this->diffuse = diff;
	this->reflect = refl;
	this->refract = refr;
	this->emit = emit;
	this->ior = ior;
	this->has_diff_map = false;

	this->energy = diffuse+reflect+refract+emit;
}

material::~material(){}

vec4 material::get_diffuse_color(vec2 uv)
{
	if(has_diff_map)
	{
		vec4 map_color = get_map_color(uv);
		return map_color*diffuse_color*diffuse + map_color*emit_color*emit;
	}
	else
	{
		return diffuse_color*diffuse + emit_color*emit;
	}
}

vec4 material::get_reflect_color()
{
	return reflect_color*reflect;
}

vec4 material::get_refract_color()
{
	return diffuse_color*refract; // TODO: change to refract_color
}

void material::add_diffuse_map(pixelmap &p)
{
	this->has_diff_map = true;
	this->diff_map = p;
}

vec4 material::get_map_color(vec2 uv)
{
	return this->diff_map.getpixel(uv.u, uv.v);
}

double material::get_ior(double mult, int type)
{
	// double wavelength = .4 +.3*mult;

	// coefficients for borosilicate crown glass (BK7)
	// double B1 = 1.03961212;
	// double B2 = 0.231792344;
	// double B3 = 1.01046945;
	// double C1 = .00600069867;
	// double C2 = .0200179144;
	// double C3 = .0103560653;

	// double lambda2 = wavelength*wavelength;
	// sellmeier equation
	// double n2 = 1 + (B1*lambda2)/(lambda2-C1) + (B2*lambda2)/(lambda2-C2) + (B3*lambda2)/(lambda2-C3);
	// return sqrt(n2);

	// bezier spline control
	// float P0 =  0.f;
	// float P1 = 0.01f;
	// float P2 = 1.0f;
	// float P3 = 1.0f;
	// float t = mult;
	// float comp_t = 1.f - mult;
	// float val = P0*comp_t*comp_t + P1*2.f*comp_t*t + P2*t*t;

	float range = 0.06;

	return ior + (mult-0.5)*range;
}

double material::get_diffuse()
{
	return diffuse;
}

double material::get_reflect()
{
	return reflect;
}

double material::get_refract()
{
	return refract;
}

double material::get_emit()
{
	return emit;
}

double material::get_energy()
{
	return energy;
}
