#include "raytracer.h"

const double raytracer::ray_tolerance = 0.000001;
const short raytracer::DIFF = 0;
const short raytracer::REFL = 1;
const short raytracer::REFR = 2;

raytracer::raytracer(std::string scenefile)
{
	max_depth = 5;
	ambience = 0.5;
	focal_length = 3.43;
	// f/1.4 f/2 f/2.8 f/4 f/5.6 f/8
	// f_stop = 4.0;
	lense_radius = 0.033;

	rgb_total = vec4(0,0,0);
	rgb_denied = vec4(0,0,0);
	
	load_scene(scenefile);
}

raytracer::~raytracer()
{
	for(int i=0; i<objs.size(); i++)
	{
		free(objs[i]);
	}
	for(int i=0; i<mtls.size(); i++)
	{
		free(mtls[i]);
	}
}

void raytracer::load_scene(std::string& scenefile)
{
	FILE *objfile = fopen(scenefile.c_str(), "r");
	if(objfile==NULL)
	{
		std::cout << "Error loading file " << scenefile << std::endl;
		return;
	}
	char line[256];
	while (fgets(line, sizeof(line), objfile))
	{
		std::stringstream ss;
		ss << line;
		std::string tok;
		ss >> tok;
		if(tok[0]=='#')
		{
			continue;
		}
		// TODO: environment info -- path depth, environment light
		if(tok=="c")
		{
			std::string width_str, height_str, fov_str, samples_str;
			ss >> width_str >> height_str >> fov_str >> samples_str;
			this->image_width = atoi(width_str.c_str());
			this->image_height = atoi(height_str.c_str());
			this->samples1D = atoi(fov_str.c_str());
			this->fov = atof(samples_str.c_str());
		}
		if(tok=="m")
		{
			vec4 diff = read_vector(ss);
			vec4 refl = read_vector(ss);
			vec4 refr = read_vector(ss);
			vec4 emit = read_vector(ss);
			std::string sdiff, srefl, srefr, semit, sior;
			ss >> sdiff >> srefl >> srefr >> semit >> sior;
			material *mat = new material(diff, refl, refr, emit, atof(sdiff.c_str()), atof(srefl.c_str()), 
									   atof(srefr.c_str()), atof(semit.c_str()), atof(sior.c_str()));
			char peek = ss.peek();
			while(peek==' ')
			{
				ss.get();
				peek = ss.peek();
			}
			if(ss.peek()!='\n')
			{
				std::string map;
				ss >> map;
				pixelmap p(map);
				mat->add_diffuse_map(p);
			}

			this->mtls.push_back(mat);
		}
		if(tok=="o")
		{
			std::string objfile, smtl;
			ss >> objfile >> smtl;
			mesh *m = new mesh(objfile, atoi(smtl.c_str()));
			this->objs.push_back(m);
			// TODO: support for mesh lights
			// if(m->shader.emit>0.0)
			// {
			// 	this->lights.push_back(m);
			// }
		}
		if(tok=="s")
		{
			std::string sradius, smtl;
			ss >> sradius >> smtl;
			vec4 s_pos = read_vector(ss);
			sphere *s = new sphere(atof(sradius.c_str()), s_pos, atoi(smtl.c_str()));
			this->objs.push_back(s);
			//check for emission
			if(mtls[s->get_mtl_idx()]->get_emit()>0.0)
			{
				this->lights.push_back(s);
			}
		}
		if(tok=="l")/*deprecated*/
		{
			std::string sradius;
			ss >> sradius;
			vec4 l_pos = read_vector(ss);
			vec4 emit = read_vector(ss);
		}
	}
}
void raytracer::print_scene_info(std::ostream& out)
{
	for(int i=0; i<mtls.size(); i++)
	{
		out << "material " << i << std::endl;
		out << "  diffuse: " << mtls[i]->get_diffuse() << ", reflect: " << mtls[i]->get_reflect() <<
		       ", refract: " << mtls[i]->get_refract() << ", emit: " << mtls[i]->get_emit() << std::endl;
	}

	for(int i=0; i<objs.size(); i++)
	{
		objs[i]->print_info(out);
	}

	out << "camera" << std::endl;
	out << "  resolution: " << image_width << " x " << image_height << std::endl;
	out << "  samples: " << samples1D << " x " << samples1D << std::endl;

}

vec4 raytracer::read_vector(std::stringstream& ss)
{
	std::string x, y, z;
	ss >> x >> y >> z;
	return vec4(atof(x.c_str()), atof(y.c_str()), atof(z.c_str()));
}

int raytracer::intersect_scene(ray& v)
{
	int closest_obj = -1;
	double min_dist = std::numeric_limits<double>::infinity();
	for(int o=0; o<objs.size(); o++)
	{
		if(objs[o]->intersect(v))
		{
			if(v.t<min_dist)
			{
				min_dist = v.t;
				closest_obj = o;
			}
		}
	}
	return closest_obj;
}

int raytracer::intersect_shadow(ray& v)
{
	int hit = 0;
	double light_dist = v.t - ray_tolerance;
	for(int o=0; o<objs.size(); o++)
	{
		bool hit_o = objs[o]->intersect(v);
		if(hit_o&&v.t<light_dist)
		{
			hit = 1;
			break;
		}
	}
	return hit;
}

vec4 raytracer::shade(ray& v)
{
	vec4 color;
	for(int l=0; l<lights.size(); l++)
	{
		vec4 pt = v.end();
		vec4 sw = lights[l]->c - pt;
		vec4 up_y(0,1,0);
		vec4 up_x(1,0,0);
		vec4 su = (fabs(sw.x)>0.1?up_y:up_x).cross(sw);
		su.normalize();
		vec4 sv = sw.cross(su);
		vec4 l_dist = pt - lights[l]->c;
		double cos_a_max = sqrt(1-lights[l]->r*lights[l]->r/(l_dist).dot(l_dist));
		double eps1 = randd();
		double eps2 = randd();
		double cos_a = 1-eps1+eps1*cos_a_max;
		double sin_a = sqrt(1-cos_a*cos_a);
		double phi = 2*M_PI*eps2;
		vec4 l_dir = su*cos(phi)*sin_a + sv*sin(phi)*sin_a + sw*cos_a;
		l_dir.normalize();

		vec4 l_org(pt);
		l_org += (l_dir*ray_tolerance);
		ray vs(l_org, l_dir);

		this->lights[l]->intersect(vs);
		int hit = intersect_shadow(vs);

		double shadow_mult = 0.0;
		if(!hit)
		{
			//std::cout << "hit light" << std::endl;
			double omega = 2*M_PI*(1-cos_a_max);
			shadow_mult = omega;//*(1.0/PI);
			// shadow_mult = 1.0;
		}
		
		vec4 light(sw);
		light.normalize();
		// pt += (light*ray_tolerance);
		ray s(pt, light);
		// double shadow_mult = trace_shadow(s, self, l);
		double n_dot_l = v.hit_norm.dot(light);
		vec4 reflect = v.hit_norm*(2*n_dot_l) - light;
		vec4 eye = v.d*(-1);
		reflect.normalize();
		// color += v.hit_color*ambience + v.hit_color * lights[l]->diffuse() * std::max(0.0, n_dot_l);
		vec4 light_color = mtls[lights[l]->get_mtl_idx()]->get_diffuse_color();
		// color += v.hit_color * lights[l]->diffuse() * std::max(0.0, n_dot_l);
		color += v.hit_color * light_color * std::max(0.0, n_dot_l);
		// double phong = objs[self]->shader.specular;
		// if(phong!=0)
		// 	color  += lights[l]->diffuse() * pow(std::max(0.0, eye.dot(reflect)), phong);
		color *= shadow_mult;

	}
	return color;


	// vec4 color;
	// for(int l=0; l<lights.size(); l++)
	// {
	// 	vec4 pt = v.end();
	// 	vec4 light = lights[l]->c - pt;
	// 	light.normalize();
	// 	// pt += (light*r.);
	// 	ray s(pt, light);
	// 	double shadow_mult = trace_shadow(s, self, l);
	// 	double n_dot_l = v.hit_norm.dot(light);
	// 	vec4 reflect = v.hit_norm*(2*n_dot_l) - light;
	// 	vec4 eye = v.d*(-1);
	// 	reflect.normalize();
	// 	color += v.hit_color*ambience + v.hit_color * lights[l]->diffuse() * std::max(0.0, n_dot_l);
	// 	double phong = objs[self]->shader.specular;
	// 	if(phong!=0)
	// 		color  += lights[l]->diffuse() * pow(std::max(0.0, eye.dot(reflect)), phong);
	// 	color *= shadow_mult;
	// 	// if(depth==0)
	// 	// 	color += ambience;
	// }
	// return color;
}

vec4 raytracer::trace_path(ray& v, int depth)
{
	if(v.debug) std::cout << std::endl;
	vec4 color;
	if(depth>max_depth) 
	{
		if(v.debug) std::cout << "MAX BOUNCES " << depth << std::endl;
		// if(v.prior_ior!=0) color += ambience;
		// if(v.prior_ior!=0 && depth < max_depth+4)
		// {

		// }
		// else
		// {
			return color; // TODO: russian roulette
		// }
	}

	int closest_obj = intersect_scene(v);
	
	if(closest_obj!=-1)
	{
		v.hit_color = mtls[v.hit_mtl]->get_diffuse_color(v.hit_uv);
		if(mtls[v.hit_mtl]->get_emit())
		{
			// if(v.stack[0]==DIFF && v.stack[1]==REFR) v.hit_color *= 6.0;
			if(v.debug) std::cout << "EMIT " << depth << std::endl;
			// std::cout << "we hit the light: " << objs[closest_obj]->diffuse().x << " " << objs[closest_obj]->diffuse().y << " " << objs[closest_obj]->diffuse().z << std::endl;
			return v.hit_color;//*objs[closest_obj]->shader.emit;
			// return objs[closest_obj]->diffuse();
		}
		v.hit_norm.normalize();
		color += shade(v);

		vec4 org = v.end();
		vec4 inc = v.d*(-1);
		vec4 refl_color, fract_color;
		double n1, n2;
		ray r = v.inherit();
		bool into = true;
		if(v.hit_norm.dot(v.d)>0)
		{
			v.hit_norm *= -1;
			into = false;
			// if(v.prior_ior<=1.000001) 
			// {
			// 	std::cout << "error$@! wl: " << v.spectrum.get_range_low() << " - " << v.spectrum.get_range_high() << std::endl;
			// }
		}
		// if(v.refract_obj == closest_obj)
		// {
		// 	v.hit_norm *= -1;
		// 	r.refract_obj = -1;
		// }

		double rng = randd()*mtls[v.hit_mtl]->get_energy();
		if (rng < mtls[v.hit_mtl]->get_refract())
		{
// std::cout << "ERROR" << std::endl;
			// refraction
			// int spectrum = randi(0,3);
			// vec4 dispersion_color = vec4(1,1,1);
			// if(spectrum==0) dispersion_color = vec4(1,0,0);
			// else if(spectrum==1) dispersion_color = vec4(0,1,0);
			// else dispersion_color = vec4(0,0,1);

			// for(int i=0; i<v.spectrum.get_num_samples(); i++)
			// {

			// vec4 dispersion_color = v.spectrum.get_color_sample(i);
			// wavelength spectrum(1, dispersion_color);
			// double min_ior = mtls[v.hit_mtl]->get_ior(v.spectrum.get_range_low());
			// double max_ior = mtls[v.hit_mtl]->get_ior(v.spectrum.get_range_high()); // TODO

			n1 = v.prior_ior;
			n2 = (!into) ? 1.0 : mtls[v.hit_mtl]->get_ior(v.spectrum.get_wavelength());//ior;
			// n2 = (!into) ? 1.0 : min_ior;


			vec4 inc_neg = inc*-1;
			double cos2t;
			// inc_neg.refract(v.hit_norm, n1, n2, cos2t);

			// if(cos2t>=0 && into && v.spectrum.get_num_samples()>1)
			// {
			// 	n2 = max_ior;
			// 	inc_neg.refract(v.hit_norm, n1, n2, cos2t);
			// }

			// if(cos2t<0) rgb_denied += v.spectrum.get_total_color();
			double R = 0.0; // always refract TODO: fresnel reflection
			double fresnel_rng = randd();

			// double R0s = ((n1 - n2)*(n1 - n2))/((n1+n2)*(n1+n2));
			// double R0 = R0s;//*R0s;
			// double R = R0 + (1-R0)*pow((1-v.hit_norm.dot(inc)), 5.0);
			// double fresnel_rng = randd();
			// double iof = n1/n2;
			// double c1 = -v.hit_norm.dot(v.d);
			// double cos2t = 1-iof*iof*(1-c1*c1);

			if (fresnel_rng > R )//&& cos2t >= 0)
			{
				for(int i=0; i<v.spectrum.get_num_samples(); i++)
				{
					if(v.debug) std::cout << "REFRACT " << depth << std::endl;
				ray r_refract = v.inherit();
				double rand_dispersion = randd();
				// double norm_wl;
				// vec4 dispersion_color = v.spectrum.get_color_sample(i, rand_dispersion, norm_wl);
				wavelength spectrum = v.spectrum.get_sample(i, rand_dispersion);
				// wavelength spectrum(1, dispersion_color, norm_wl);

				// dispersion_color = v.spectrum.get_color_sample(i+1);
				// norm_wl = i*0.3333333;
				// if (v.spectrum.get_num_samples()>1)
				// rgb_total += dispersion_color;
				n2 = (!into) ? 1.0 : mtls[v.hit_mtl]->get_ior(spectrum.get_wavelength());
				// std::cout << "ior " << n2 << " next color: " << dispersion_color.x << " " << dispersion_color.y << " " << dispersion_color.z << std::endl;
				// R0s = ((n1 - n2)*(n1 - n2))/((n1+n2)*(n1+n2));
				// R0 = R0s;//*R0s;
				//  R = R0 + (1-R0)*pow((1-v.hit_norm.dot(inc)), 5.0);
				//  fresnel_rng = randd();
				//  iof = n1/n2;
				//  c1 = -v.hit_norm.dot(v.d);
				//  cos2t = 1-iof*iof*(1-c1*c1);
				// //refract
				// // vec4 fract;
				// double c2 = sqrt(cos2t);
				// // if(n2!=1 || n2!=1.3)
				// // if(n2>1.3)
				// // std::cout << "n2: " << n2 << " i: " << i << " color: " << dispersion_color.x << " " << dispersion_color.y << " " << dispersion_color.z << std::endl;
				// vec4 dir = v.d*iof + v.hit_norm * (iof*c1 - c2);

				vec4 dir = inc_neg.refract(v.hit_norm, n1, n2, cos2t);
				// if(cos2t<0) std::cout << "total internal reflection" << std::endl;
	
				if(cos2t<0) 
				{
					n2 = v.prior_ior;
					dir = inc.reflect(v.hit_norm);
				}
				dir.normalize();
// 				if (v.debug) std::cout << v.refract_bounces << std::endl;
// 				if (v.refract_bounces==(max_depth-2))
// 				{
// 					if(v.debug) std::cout << "too many internal bounces" << std::endl;
// 					fract = v.d;
// 				}
// if(v.debug)
// {
// 	std::cout << "REFRACT" << std::endl;
// 	std::cout << "c1: " << c1 << std::endl;
// 	// std::cout << "c2: " << c2 << std::endl;
// 	std::cout << "c2 = " << "sqrt(1 - " << iof << "*" << iof << "*(1 - " << c1 << "*" << c1 << "))" << std::endl;
// 	std::cout << "v.d: " << v.d.x << " " << v.d.y << " " << v.d.z << std::endl;
// 	std::cout << "origin: " << org.x << " " << org.y << " " << org.z << std::endl;
// }
				org += (dir*ray_tolerance);
				r_refract.o = org;
				r_refract.d = dir;
				r_refract.prior_ior = n2;
				r_refract.spectrum = spectrum;
				r_refract.stack[depth] = REFR;
				// r_refract.debug = v.debug;
				// ray f(org, fract);
				// f.refract_obj = (v.refract_obj==closest_obj) ? -1.0 : closest_obj; 
				// f.refract_bounces = v.refract_bounces+1;
				// f.debug = v.debug;
				color += spectrum.get_total_color()*mtls[v.hit_mtl]->get_refract_color()*trace_path(r_refract, depth+1);				
				// TODO: attenuation

				}
			}
			else
			{
				if(v.debug) std::cout << "FRESNEL REFLECT " << depth << std::endl;
				//reflect
				// vec4 dir = v.hit_norm*(v.hit_norm.dot(inc)*2) - inc;
				vec4 dir = inc.reflect(v.hit_norm);
				dir.normalize();
				// if(depth>2) dir = v.d;
				// dir = v.d;
				// depth--;
				org += (dir*ray_tolerance);
				r.o = org;
				r.d = dir;
				r.stack[depth] = REFL;
				// r.spectrum = v.spectrum;
				// r.debug = v.debug;
				// r.prior_ior = v.prior_ior;//mtls[v.hit_mtl]->get_ior(v.spectrum.get_wavelength());
				// r.refract_obj = (v.refract_obj==closest_obj) ? -1.0 : closest_obj;
				color += mtls[v.hit_mtl]->get_refract_color()*trace_path(r, depth+1);
				// break;
			}
		// }
		}
		else if(rng < mtls[v.hit_mtl]->get_refract() + mtls[v.hit_mtl]->get_reflect())
		{
if(v.debug) std::cout << "REFLECT " << depth << std::endl;
// r.debug = v.debug;
// std::cout << "ERROR" << std::endl;
			// reflect
			// vec4 dir = v.hit_norm*(v.hit_norm.dot(inc)*2) - inc;
			vec4 dir = inc.reflect(v.hit_norm);
			org += (dir*ray_tolerance);
			vec4 up(0.0,1.0,0.0);
			vec4 axis_x = dir.cross(up);
			vec4 axis_y = axis_x.cross(dir);
			double dist_radius = 0.0; // TODO: per material glossy reflection control
			
			double cur_theta = randd()*2.0*M_PI;
			double cur_radius = randd()*dist_radius;
			double cur_x = cur_radius*cos(cur_theta);
			double cur_y = cur_radius*sin(cur_theta);
			dir += (axis_x*cur_x);
			dir += (axis_y*cur_y);
			dir.normalize();
			r.o = org;
			r.d = dir;
			r.stack[depth] = REFL;
			// r.spectrum = v.spectrum;
			// r.debug = v.debug;
			// r.prior_ior = into ? mtls[v.hit_mtl]->get_ior(v.spectrum.get_wavelength()) : 1.0;
			// r.prior_ior = v.prior_ior;
			color += mtls[v.hit_mtl]->get_reflect_color()*trace_path(r, depth+1);

			// for(int i=0; i<reflect_samples; i++)
			// {
			// 	ray cur_r(r);
			// 	double cur_theta = randd()*2.0*PI;
			// 	double cur_radius = randd()*dist_radius;
			// 	double cur_x = cur_radius*cos(cur_theta);
			// 	double cur_y = cur_radius*sin(cur_theta);
			// 	vec4 cur_dir(dir);
			// 	cur_dir += (axis_x*cur_x);
			// 	cur_dir += (axis_y*cur_y);
			// 	cur_dir.normalize();
			// 	cur_r.o = reflect_org;
			// 	cur_r.d = cur_dir;
			// 	total_refl_color += objs[closest_obj]->reflect()*trace_ray(cur_r, depth+1, closest_obj);
			// }
			// refl_color = total_refl_color*(1.0/reflect_samples);
		}
		else
		{
if(v.debug) std::cout << "DIFF " << depth << std::endl;
// r.debug = v.debug;
			// diffuse
			double r1 = 2*M_PI*randd();
			double r2 = randd(), r2s = sqrt(r2);
			vec4 w_axis(v.hit_norm);
			vec4 up_y(0,1,0);
			vec4 up_x(1,0,0);
			vec4 up = fabs(w_axis.x)>0.1?up_y:up_x;
			vec4 u_axis = up.cross(w_axis);
			u_axis.normalize();
			vec4 v_axis = w_axis.cross(u_axis);
			vec4 dir = (u_axis*cos(r1)*r2s + v_axis*sin(r1)*r2s + w_axis*sqrt(1-r2));
			dir.normalize();

			// if (v.hit_norm.dot(dir)<0.0) std::cout << "ERROR!" << std::endl;
			
			//diffuse
			// vec4 dir(v.hit_norm);
			// vec4 up(0.0,1.0,0.0);
			// if(dir == up || dir == (up*-1)) up.x = 1.0;
			// vec4 axis_x = dir.cross(up);
			// vec4 axis_z = axis_x.cross(dir);
			// axis_x.normalize();
			// axis_z.normalize();

			// double rand1 = randd();
			// double rand2 = randd();
			// double theta = acos(-sqrt(rand1));
			// double phi = 2.0*PI*rand2;
			// // double theta = rand1*PI/4.0;
			// double cur_x = sin(theta)*cos(phi);
			// double cur_z = sin(theta)*sin(phi);
			// double cur_y = cos(theta);
			// vec4 hemi(cur_x, cur_y, cur_z);
			// hemi.normalize();
			// // if(hemi.length()!=1.0) std::cout << "ERROR: not unit hemisphere" << std::endl;
			// vec4 straight(0.0, 1.0, 0.0);
			// double proj = hemi.dot(straight);
			// dir *= proj;
			// dir += (axis_x*hemi.x);
			// dir += (axis_z*hemi.z);
			// dir.normalize();
			// dir *= -1.0;

// std::cout << "ray: " << dir.x << " " << dir.y << " " << dir.z << std::endl;
// std::cout << "u_axis" << u_axis.x << " " << u_axis.y << " " << u_axis.z << std::endl;
// std::cout << "w_axis" << w_axis.x << " " << w_axis.y << " " << w_axis.z << std::endl;
// std::cout << "v_axis" << v_axis.x << " " << v_axis.y << " " << v_axis.z << std::endl;
// std::cout << " axis_x: " << axis_x.x << " " << axis_x.y << " " << axis_x.z << std::endl;
// std::cout << " axis_z: " << axis_z.x << " " << axis_z.y << " " << axis_z.z << std::endl;
// std::cout << "up: " << up.x << " " << up.y << " " << up.z << std::endl;
// std::cout << "norm: " << v.hit_norm.x << " " << v.hit_norm.y << " " << v.hit_norm.z << std::endl;
			org += (dir*ray_tolerance);
			r.o = org;
			r.d = dir;
			r.stack[depth] = DIFF;
			// r.spectrum = v.spectrum;
			// r.debug = v.debug;
			// r.prior_ior = v.prior_ior;
			// r.prior_ior = into ? mtls[v.hit_mtl]->get_ior(v.spectrum.get_wavelength()) : 1.0;
			color += v.hit_color*trace_path(r, depth+1);//*(pow(0.5, depth));
		}
	}
	else
	{
		if(v.debug) std::cout << "ENV LIGHT " << depth << std::endl;
		color += ambience;//vec4(0.385, 0.457, 0.55);
// if(v.debug) 
// {
// 	std::cout << "NOTHING" << std::endl;
// 	std::cout << "origin: " << v.o.x << " " << v.o.y << " " << v.o.z << std::endl;
// 	std::cout << "direction: " << v.d.x << " " << v.d.y << " " << v.d.z << std::endl;
// }
	}
	return color;
}

// returns a random double between 0 and 1
double raytracer::randd()
{
	// return drand48();
	int tid = 0;
	#ifdef _OPENMP
	tid = omp_get_thread_num();
	#endif
	return ((double)rand_r(&(rand_seed[tid])))/((double)RAND_MAX);
}

// returns a random double between -0.5 and 0.5
double raytracer::randd_negative()
{
	int tid = 0;
	#ifdef _OPENMP
	tid = omp_get_thread_num();
	#endif
	return ((double)(rand_r(&(rand_seed[tid])) - RAND_MAX/2))/((double)RAND_MAX);
}

// returns a random integer between [lo, hi]
int raytracer::randi(int lo, int hi)
{
	int tid = 0;
	#ifdef _OPENMP
	tid = omp_get_thread_num();
	#endif
	return rand_r(&(rand_seed[tid])) % (lo-hi) + lo;
}

void raytracer::init_rand(int num_threads)
{
	srand(time(NULL));
	rand_seed = new unsigned int[num_threads];
	for(int i=0; i<num_threads; i++)
	{
		rand_seed[i] = rand();
	}
}

void raytracer::render_image(std::string& outfile)
{
	// define camera properties
	// int image_width = 1168;
	// int image_height = 657;
	vec4 look_from(0,0,1);
	double r_fov = fov * (M_PI/180.0);// radians

	pixelmap pic(image_width, image_height);
	double aspect_ratio = (double)image_height/(double)image_width;
	double samples = samples1D*samples1D;

	// calculate image plane size in world space
	double scene_width = tan(r_fov/2)*2;
	double scene_height = tan((r_fov*aspect_ratio)/2)*2;
	double pixel_width = scene_width/image_width;
	double pixel_slice = pixel_width/samples1D;

	int total_pixels = image_width*image_height;
	double last_report = 0;

	int thread_count = 1;
	#ifdef _OPENMP 
	thread_count = omp_get_max_threads();
	#endif
	std::cout << "threads: " << thread_count << std::endl;
	init_rand(thread_count);

	#pragma omp parallel for schedule(dynamic, 2)
	for(int i=0; i<image_width; i++)
	{
		for(int j=0; j<image_height; j++)
		{
			int ii = i - image_width/2;
			int jj = j - image_height/2;
			double x = ii*pixel_width + pixel_slice/2;
			double y = jj*pixel_width + pixel_slice/2;
			
			vec4 color(0,0,0);
			// for each sample
			for(int m=0; m<samples1D; m++)
			{
				for(int n=0; n<samples1D; n++)
				{
					double xoff = randd_negative()*pixel_slice;
					double yoff = randd_negative()*pixel_slice;
					vec4 dir(x+(pixel_slice*m)+xoff,y+(pixel_slice*n)+yoff,0);
					dir -= look_from;
					dir.normalize();

					vec4 look_at = look_from + dir*focal_length;

					vec4 org(look_from);
					double rr = sqrt(randd());
					double rtheta = 2*M_PI*randd();
					org.x = rr*cos(rtheta)*lense_radius;
					org.y = rr*sin(rtheta)*lense_radius;
					dir = look_at - org;
					dir.normalize();

					ray v(org, dir);
// if(i==696&&j==346)
// {
// 	v.debug=1;
// }
					color += trace_path(v, 0);
				}
			}
			color *= 1/samples;
			color.clamp(0.0, 1.0);


			// double x = (double)i/(double)image_width;
			// wavelength wl;
			// vec4 color = wl.get_color_from_mult(x);
			// rgb_total += color;
			// if(x<.34&&x>.32) std::cout << color.x << " " << color.y << " " << color.z << std::endl;
			if(i==835&&j==346) std::cout << color.x << " " << color.y << " " << color.z << std::endl;
			// else color = vec4(0,0,0);

			pic.setpixel(i, j, color);
			// double decimal = ((double)((i+1)*(j+1)))/((double)(total_pixels));
			// // double decimal = ((double)((image_width - i)*(j+1)))/((double)(total_pixels));

			// if(decimal-last_report > 0.1)
			// {
			// 	last_report = decimal;
			// 	std::cout << (int)(decimal*100.0) << "\%" << std::endl;
			// 	pic.writeppm(outfile);
			// }
			rand_total += randd();
			
		}
	}
	pic.writeppm(outfile);
	// std::cout << rgb_total.x << " " << rgb_total.y << " " << rgb_total.z << std::endl;
	std::cout << rgb_denied.x << " " << rgb_denied.y << " " << rgb_denied.z << std::endl;
	std::cout << rand_total/((double)total_pixels) << std::endl;
}

double when()
{
    struct timeval tp;
    gettimeofday(&tp, NULL);
    return ((double) tp.tv_sec + (double) tp.tv_usec * 1e-6);
}

// g++ *.cpp -O3 -fopenmp -o raytracer
int main(int argc, char *argv[])
{
	std::string scenefile;
	std::string outfile = "test.ppm";
	if(argc < 2)
	{
		std::cout << "No scene file parameter: Loading scenes/test.scene" << std::endl;
		scenefile = "scenes/test.scene";
	}
	else
	{
		scenefile = argv[1];
	}
	if(argc < 3)
	{
		std::cout << "No output file parameter: saving image to test.ppm" << std::endl;
		outfile = "test.ppm";
	}
	else
	{
		outfile = argv[2];
	}
	raytracer r(scenefile);
	r.print_scene_info(std::cout);

	// start timer
	double start = when();

	r.render_image(outfile);

	double time_gone_by = when() - start;
	std::cout << "TOTAL TIME: " << time_gone_by << " seconds" << std::endl;
	
}