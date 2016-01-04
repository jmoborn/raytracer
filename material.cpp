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
	// if(energy!=1.0)
	// {
	// 	diffuse /= energy;
	// 	reflect /= energy;
	// 	refract /= energy;
	// }
	// if(this->emit>0.0)
	// {
	// 	this->diffuse_color = this->emit_color;
	// 	this->diffuse = this->emit;
	// }
}

material::~material(){}

vec4 material::get_diffuse_color(vec2 uv)
{
	if(has_diff_map)
	{
		// std::cout << uv.u << ", " << uv.v << std::endl;
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

double material::get_ior(double mult)
{
	
	// wavelength spectrum(1, color);
	// double wavelength = spectrum.get_wavelength();

	// double wavelength = .4 +.3*mult;

	// double B1 = 1.03961212;
	// double B2 = 0.231792344;
	// double B3 = 1.01046945;
	// double C1 = .00600069867;
	// double C2 = .0200179144;
	// double C3 = .0103560653;

	// double lambda2 = wavelength*wavelength;

	// double n2 = 1 + (B1*lambda2)/(lambda2-C1) + (B2*lambda2)/(lambda2-C2) + (B3*lambda2)/(lambda2-C3);
	// // std::cout << "ior: " << sqrt(n2) << std::endl;
	// return sqrt(n2);


	// double wavelength = 0.01*(2*color.x - color.z + color.y);
	// if(wavelength>0) 
	// std::cout << "get_ior: " << ior << " + " << wavelength << " = " << ior + wavelength << std::endl;

	// double wavelength = 0.06*mult*(color.x - color.z);
	// return ior + wavelength;

	// std::cout << "mult: " << mult << std::endl;
	// return ior + (mult-0.5)*0.000009;
	return ior + (mult-0.5)*0.06;
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
