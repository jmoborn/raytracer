#include "raytracer.h"

const double raytracer::ray_tolerance = 0.0000000001;
const short raytracer::DIFF = 0;
const short raytracer::REFL = 1;
const short raytracer::REFR = 2;
const short raytracer::USER = 3;

raytracer::raytracer(std::string scenefile)
{
	max_depth = 5;
	min_depth = 5;
	ambience = 0.5;
	focal_length = 3.431;
	lense_radius = 0.033;
	
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

	out << "lights: " << lights.size() << std::endl;

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
	ray vs(v.o, v.d);
	for(int o=0; o<objs.size(); o++)
	{
		bool hit_o = objs[o]->intersect(vs);
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
		vec4 l_dist = pt - lights[l]->c;
		double cos_a_max = sqrt(1-lights[l]->r*lights[l]->r/(l_dist).dot(l_dist));

		vec4 l_dir = lights[l]->c - pt;
		vec4 up_y(0,1,0);
		vec4 up_x(1,0,0);
		vec4 up = fabs(l_dir.x)>0.1?up_y:up_x;
		vec4 axis_x = l_dir.cross(up);
		axis_x.normalize();
		vec4 axis_y = axis_x.cross(l_dir);
		axis_y.normalize();
		double hyp = 1.0/cos_a_max;
		double dist_radius = sqrt(hyp*hyp-1);
		double rr = sqrt(randd());
		double rtheta = 2*M_PI*randd();

		double cur_x = rr*cos(rtheta);
		double cur_y = rr*sin(rtheta);
		l_dir += (axis_x*cur_x);
		l_dir += (axis_y*cur_y);
		l_dir.normalize();

		vec4 l_org(pt);
		l_org += (l_dir*ray_tolerance);
		ray vl(l_org, l_dir);
		this->lights[l]->intersect(vl);
		double light_dist = vl.t - ray_tolerance;

		ray vs(l_org, l_dir);
		int hit = intersect_scene(vs);
		if(vs.t>light_dist) hit = -1;

		double shadow_mult = 0.0;
		if(hit<0)
		{
			double omega = 2*M_PI*(1-cos_a_max);
			shadow_mult = omega;
		}

		double n_dot_l = v.hit_norm.dot(l_dir);
		vec4 light_color = mtls[lights[l]->get_mtl_idx()]->get_diffuse_color();
		color += v.hit_color * light_color * std::max(0.0, n_dot_l) * shadow_mult; // * M_1_PI;
	}
	return color;
}

vec4 raytracer::trace_path(ray& v, int depth)
{
	if(v.debug) std::cout << std::endl;
	vec4 color;
	if(depth>min_depth) 
	{
		if(v.debug) std::cout << "MAX BOUNCES " << depth << std::endl;
		return color;
	}

	int closest_obj = intersect_scene(v);
	
	if(closest_obj!=-1)
	{
		v.hit_color = mtls[v.hit_mtl]->get_diffuse_color(v.hit_uv);
		if(mtls[v.hit_mtl]->get_emit())
		{
			if(v.debug) std::cout << "EMIT " << depth << std::endl;			
			return v.hit_color;
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
		}

		double rng = randd()*mtls[v.hit_mtl]->get_energy();
		if (rng < mtls[v.hit_mtl]->get_refract())
		{
			n1 = v.prior_ior;
			n2 = (!into) ? 1.0 : mtls[v.hit_mtl]->get_ior(v.spectrum.get_wavelength());
			vec4 inc_neg = inc*-1;

			for(int i=0; i<v.spectrum.get_num_samples(); i++)
			{
				if(v.debug) std::cout << "REFRACT " << depth << std::endl;
				// refract
				ray r_refract = v.inherit();
				double rand_dispersion = randd();
				
				wavelength new_spectrum = v.spectrum.get_sample(i, rand_dispersion);
				n2 = (!into) ? 1.0 : mtls[v.hit_mtl]->get_ior(new_spectrum.get_wavelength(), v.stack[0]);
				double cos2t;
				double fresnel;
				vec4 dir = inc_neg.refract(v.hit_norm, n1, n2, cos2t, fresnel);
	
				if(cos2t<0 || randd()<fresnel) 
				{
					n2 = v.prior_ior;
					dir = inc.reflect(v.hit_norm);
				}
				dir.normalize();
				org += (dir*ray_tolerance);
				r_refract.o = org;
				r_refract.d = dir;
				r_refract.prior_ior = n2;
				r_refract.spectrum = new_spectrum;
				r_refract.stack[depth] = REFR;
				color += new_spectrum.get_total_color()*mtls[v.hit_mtl]->get_refract_color()*trace_path(r_refract, depth+1);				
				// TODO: attenuation
			}
		}
		else if(rng < mtls[v.hit_mtl]->get_refract() + mtls[v.hit_mtl]->get_reflect())
		{
			if(v.debug) std::cout << "REFLECT " << depth << std::endl;
			// reflect
			vec4 dir = inc.reflect(v.hit_norm);
			org += (v.hit_norm*ray_tolerance);
			vec4 up(0.0,1.0,0.0);
			vec4 axis_x = dir.cross(up);
			axis_x.normalize();
			vec4 axis_y = axis_x.cross(dir);
			axis_y.normalize();
			double dist_radius = 0.0;//0.09; // TODO: per material glossy reflection control
			
			double cur_theta = randd()*2.0*M_PI;
			double cur_radius = sqrt(randd())*dist_radius;
			double cur_x = cur_radius*cos(cur_theta);
			double cur_y = cur_radius*sin(cur_theta);
			dir += (axis_x*cur_x);
			dir += (axis_y*cur_y);
			dir.normalize();
			r.o = org;
			r.d = dir;
			r.stack[depth] = REFL;
			color += mtls[v.hit_mtl]->get_reflect_color()*trace_path(r, depth+1);
		}
		else
		{
			if(v.debug) std::cout << "DIFF " << depth << std::endl;
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

			org += (v.hit_norm*ray_tolerance);
			r.o = org;
			r.d = dir;
			r.stack[depth] = DIFF;
			color += v.hit_color*trace_path(r, depth+1);
		}
	}
	else
	{
		if(v.debug) std::cout << "ENV LIGHT " << depth << std::endl;
		color += ambience;
	}
	return color;
}

// returns a random double between 0 and 1
double raytracer::randd()
{
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
					color += trace_path(v, 0);
				}
			}
			color *= 1/samples;
			color.clamp(0.0, 1.0);

			pic.setpixel(i, j, color);

			// output progress
			// if(i%(image_width/4)==0 && j==0)
			// {
			// 	std::cout << (double)i/(double)image_width << std::endl;
			// 	pic.writeppm(outfile);
			// }
			
		}
	}
	pic.writeppm(outfile);
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
	std::string outfile;
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
		std::cout << "No output file parameter: saving image to renders/test.ppm" << std::endl;
		outfile = "renders/test.ppm";
	}
	else
	{
		outfile = argv[2];
	}
	raytracer r(scenefile);
	r.print_scene_info(std::cout);

	double start = when();

	r.render_image(outfile);

	double time_gone_by = when() - start;
	std::cout << "TOTAL TIME: " << time_gone_by << " seconds" << std::endl;
	
}