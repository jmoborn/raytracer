#include "raytracer.h"

raytracer::raytracer()
{
	srand(time(NULL));
	max_depth = 5;
	shadow_samples = 1;
	reflect_samples = 1;
	refract_samples = 1;
	samples = 3;
	ambience = 0.5;
	ray_tolerance = 0.001;
	std::string scenefile = "clock.scene";
	load_scene(scenefile);
}

raytracer::~raytracer()
{
	
}

void raytracer::load_scene(std::string& scenefile)
{
	FILE *objfile = fopen(scenefile.c_str(), "r");
	if(objfile==NULL)
	{
		std::cout << "Error loading file " << scenefile << std::endl;
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
		if(tok=="m")
		{
			std::cout << "creating material" << std::endl;
			vec4 diff = read_vector(ss);
			vec4 refl = read_vector(ss);
			std::string sdiff, sref, sspec, sfract, sior;
			ss >> sdiff >> sref >> sspec >> sfract >> sior;
			material *mat = new material(diff, refl, atof(sdiff.c_str()), atof(sref.c_str()), 
									   atof(sspec.c_str()), atof(sfract.c_str()), atof(sior.c_str()));
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
			peek = ss.peek();
			while(peek==' ')
			{
				ss.get();
				peek = ss.peek();
			}
			if(ss.peek()!='\n')
			{
				std::string emit_str;
				ss >> emit_str;
				double emit = atof(emit_str.c_str());
				mat->emit = emit;
			}

			this->mtls.push_back(mat);
		}
		if(tok=="o")
		{
			std::string objfile, smtl;
			ss >> objfile >> smtl;
			mesh *m = new mesh(objfile, *this->mtls[atoi(smtl.c_str())]);
			this->objs.push_back(m);
		}
		if(tok=="s")
		{
			std::cout << "creating sphere" << std::endl;
			std::string sradius, smtl;
			ss >> sradius >> smtl;
			vec4 s_pos = read_vector(ss);
			sphere *s = new sphere(atof(sradius.c_str()), s_pos, *this->mtls[atoi(smtl.c_str())]);
			std::cout << "radius " << s->r << std::endl;
			std::cout << "position " << s->c.x << " " << s->c.y << " " << s->c.z << std::endl;
			vec4 tcol = s->diffuse();
			std::cout << "diffuse " << tcol.x << " " << tcol.y << " " << tcol.z << std::endl;
			this->objs.push_back(s);
		}
		if(tok=="l")
		{
			std::cout << "creating light" << std::endl;
			std::string sradius;
			ss >> sradius;
			vec4 l_pos = read_vector(ss);
			vec4 emit = read_vector(ss);
			sphere *l = new sphere(atof(sradius.c_str()), l_pos);
			std::cout << "radius " << l->r << std::endl;
			std::cout << "position " << l->c.x << " " << l->c.y << " " << l->c.z << std::endl;
			l->shader.diffuse_color = emit;
			vec4 tcol = l->diffuse();
			std::cout << "diffuse " << tcol.x << " " << tcol.y << " " << tcol.z << std::endl;
			l->shader.emit = tcol.length();
			this->lights.push_back(l);
			this->objs.push_back(l);
		}
	}
}

vec4 raytracer::read_vector(std::stringstream& ss)
{
	std::string x, y, z;
	ss >> x >> y >> z;
	return vec4(atof(x.c_str()), atof(y.c_str()), atof(z.c_str()));
}

vec4 raytracer::shade(ray& v, int depth, int self)
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
		double phi = 2*PI*eps2;
		vec4 l_dir = su*cos(phi)*sin_a + sv*sin(phi)*sin_a + sw*cos_a;
		l_dir.normalize();

		vec4 l_org(pt);
		l_org += (l_dir*ray_tolerance);
		ray vs(l_org, l_dir);
		int hits = 0;

		this->lights[l]->intersect(vs);
		double light_dist = vs.t;
		for(int o=0; o<objs.size(); o++)
		{
			bool hit_o = objs[o]->intersect(vs);
			if(hit_o&&vs.t<light_dist)
			{
				hits++;
				break;
			}
		}
		double shadow_mult = 0.0;
		if(!hits)
		{
			//std::cout << "hit light" << std::endl;
			double omega = 2*PI*(1-cos_a_max);
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
		color += v.hit_color * lights[l]->diffuse() * std::max(0.0, n_dot_l);
		double phong = objs[self]->shader.specular;
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
	// 	// pt += (light*ray_tolerance);
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

vec4 raytracer::trace_path(ray& v, int depth, int self)
{
	vec4 color;
	if(depth>max_depth) return color;
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
	if(closest_obj!=-1)
	{
		if(objs[closest_obj]->shader.emit)
		{
			// std::cout << "we hit the light: " << objs[closest_obj]->diffuse().x << " " << objs[closest_obj]->diffuse().y << " " << objs[closest_obj]->diffuse().z << std::endl;
			return v.hit_color;//*objs[closest_obj]->shader.emit;
			// return objs[closest_obj]->diffuse();
		}
		v.hit_norm.normalize();
		color += shade(v, depth, closest_obj);

		vec4 org = v.end();
		vec4 inc = v.d*(-1);
		vec4 refl_color, fract_color;
		double n1, n2;
		ray r;
		bool into = true;
		if(v.hit_norm.dot(v.d)>0)
		{
			v.hit_norm *= -1;
			into = false;
		}
		// if(v.refract_obj == closest_obj)
		// {
		// 	v.hit_norm *= -1;
		// 	r.refract_obj = -1;
		// }

		double rng = randd()*objs[closest_obj]->shader.energy;
		if (rng < objs[closest_obj]->shader.refract)
		{
// std::cout << "ERROR" << std::endl;
			//refraction
			n1 = (self == -1) ? 1.0 : objs[self]->shader.ior;
			n2 = (!into) ? 1.0 : objs[closest_obj]->shader.ior;

			double R0s = ((n1 - n2)*(n1 - n2))/((n1+n2)*(n1+n2));
			double R0 = R0s*R0s;
			double R = R0 + (1-R0)*pow((1-v.hit_norm.dot(inc)), 5.0);
			double fresnel_rng = randd();
			double iof = n1/n2;
			double c1 = -v.hit_norm.dot(v.d);
			double cos2t = 1-iof*iof*(1-c1*c1);
			if (fresnel_rng > R && cos2t >= 0)
			{

				//refraction
				vec4 fract;
				double c2 = sqrt(cos2t);
				fract = v.d*iof + v.hit_norm * (iof*c1 - c2);
				fract.normalize();
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
				org += (fract*ray_tolerance);
				ray f(org, fract);
				// f.refract_obj = (v.refract_obj==closest_obj) ? -1.0 : closest_obj; 
				// f.refract_bounces = v.refract_bounces+1;
				// f.debug = v.debug;
				color += objs[closest_obj]->refract()*trace_path(f, depth+1, closest_obj);
			}
			else
			{
				//reflection
				vec4 dir = v.hit_norm*(v.hit_norm.dot(inc)*2) - inc;
				org += (dir*ray_tolerance);
				dir.normalize();
				r.o = org;
				r.d = dir;
				r.refract_obj = (v.refract_obj==closest_obj) ? -1.0 : closest_obj;
				color += objs[closest_obj]->refract()*trace_path(r, depth+1, closest_obj);
			}
		}
		else if(rng < objs[closest_obj]->shader.refract + objs[closest_obj]->shader.reflect)
		{
// if(v.debug) std::cout << "SPEC" << std::endl;
// r.debug = v.debug;
// std::cout << "ERROR" << std::endl;
			//reflect
			vec4 dir = v.hit_norm*(v.hit_norm.dot(inc)*2) - inc;
			org += (dir*ray_tolerance);
			vec4 up(0.0,1.0,0.0);
			vec4 axis_x = dir.cross(up);
			vec4 axis_y = axis_x.cross(dir);
			double dist_radius = 0.0;
			
			double cur_theta = randd()*2.0*PI;
			double cur_radius = randd()*dist_radius;
			double cur_x = cur_radius*cos(cur_theta);
			double cur_y = cur_radius*sin(cur_theta);
			dir += (axis_x*cur_x);
			dir += (axis_y*cur_y);
			dir.normalize();
			r.o = org;
			r.d = dir;
			color += objs[closest_obj]->reflect()*trace_path(r, depth+1, closest_obj);

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
// if(v.debug) std::cout << "DIFF" << std::endl;
// r.debug = v.debug;
			double r1 = 2*PI*randd();
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
			color += v.hit_color*trace_path(r, depth+1, closest_obj);//*(pow(0.5, depth));
		}
	}
	else
	{
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

vec4 raytracer::trace_ray(ray& v, int depth, int self)
{
	int closest_obj = -1;
	double min_dist = std::numeric_limits<double>::infinity();
	vec4 color;
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
	if(closest_obj!=-1)
	{
		v.hit_norm.normalize();
		color += shade(v, depth, closest_obj);
		//do we have reflection or refraction?
		if(((objs[closest_obj]->shader.reflect > 0)||(objs[closest_obj]->shader.refract > 0))&&(depth<max_depth))
		{
			vec4 org = v.end();
			vec4 inc = v.d*(-1);
			vec4 refl_color, fract_color;
			double n1, n2;
			ray r;
			if(v.refract_obj == closest_obj)
			{
				v.hit_norm *= -1;
				r.refract_obj = -1;
			}
			//calculate reflection ray
			if(objs[closest_obj]->shader.reflect > 0)
			{
				vec4 dir = v.hit_norm*(v.hit_norm.dot(inc)*2) - inc;
				vec4 reflect_org = org + (dir*ray_tolerance);
				vec4 up(0.0,1.0,0.0);
				vec4 axis_x = dir.cross(up);
				vec4 axis_y = axis_x.cross(dir);
				double dist_radius = 0.75;
				vec4 total_refl_color = vec4(0.0,0.0,0.0);
				for(int i=0; i<reflect_samples; i++)
				{
					ray cur_r(r);
					double cur_theta = randd()*2.0*PI;
					double cur_radius = randd()*dist_radius;
					double cur_x = cur_radius*cos(cur_theta);
					double cur_y = cur_radius*sin(cur_theta);
					vec4 cur_dir(dir);
					cur_dir += (axis_x*cur_x);
					cur_dir += (axis_y*cur_y);
					cur_dir.normalize();
					cur_r.o = reflect_org;
					cur_r.d = cur_dir;
					total_refl_color += objs[closest_obj]->reflect()*trace_ray(cur_r, depth+1, closest_obj);
				}
				refl_color = total_refl_color*(1.0/reflect_samples);
			}
			//calculate refraction ray
			if(objs[closest_obj]->shader.refract > 0)
			{
				n1 = (self == -1) ? 1.0 : objs[self]->shader.ior;
				n2 = (v.refract_obj==closest_obj) ? 1.0 : objs[closest_obj]->shader.ior;
				double iof = n1/n2;
				double c1 = -v.hit_norm.dot(v.d);
				double c2 = sqrt(1-iof*iof*(1-c1*c1));
				vec4 fract = v.d*iof + v.hit_norm * (iof*c1 - c2);
				fract.normalize();
				vec4 refract_org = org + (fract*ray_tolerance);
				ray f(refract_org, fract);
				f.refract_obj = (v.refract_obj==closest_obj) ? -1.0 : closest_obj; 
				fract_color = objs[closest_obj]->refract()*trace_ray(f, depth+1, closest_obj);
			}
			//fresnel effect with schlick's approximation
			if(objs[closest_obj]->shader.reflect > 0&&objs[closest_obj]->shader.refract > 0)
			{
				double R0s = ((n1 - n2)*(n1 - n2))/((n1+n2)*(n1+n2));
				double R0 = R0s*R0s;
				double R = R0 + (1-R0)*pow((1-v.hit_norm.dot(inc)), 5.0);
				refl_color *= R;
				fract_color *= (1-R);
			}
			color += (refl_color + fract_color);
		}
	}
	else
	{
		color += vec4(0.2, 0.2, 0.2);
	}
	return color;
}

//returns a random double between 0 and 1
double raytracer::randd()
{
	return ((double)rand())/((double)RAND_MAX);
}

//returns a random double between -0.5 and 0.5
double raytracer::randd_negative()
{
	return ((double)(rand() - RAND_MAX/2))/((double)RAND_MAX);
}

double raytracer::trace_shadow(ray& v, int self, int light)
{
	int hits = 0;
	for(int i=0; i<this->shadow_samples; i++)
	{
		ray vs = v;
		double x = randd_negative();
		double y = randd_negative();
		double z = randd_negative();
		vec4 offset(x,y,z);
		offset.normalize();
		offset *= this->lights[light]->r;
		vec4 new_dir = (this->lights[light]->c + offset) - vs.o;
		new_dir.normalize();
		vs.o += (new_dir*ray_tolerance);
		vs.d = new_dir;
		this->lights[light]->intersect(vs);
		double light_dist = vs.t;
		for(int o=0; o<objs.size(); o++)
		{
			bool hit_o = objs[o]->intersect(vs);
			if(hit_o&&vs.t<light_dist)
			{
				hits++;
				break;
			}
		}
	}
	
	return ((double)(this->shadow_samples - hits))/((double)this->shadow_samples);
}

const double raytracer::PI;

int main()
{
	raytracer r;

	//TODO: read this info from the scene file
	//define camera properties
	int image_width = 1168;
	int image_height = 657;
	// int image_width = 880;
	// int image_height = 660;
	// int image_width = 640;
	// int image_height = 360;
	// int image_width = 512;
	// int image_height = 512;
	vec4 look_from(0,0,1);
	double fov = 54.0;//degrees
	// double fov = 40.0;//degrees
	fov *= (r.PI/180.0);//radians
	pixelmap pic(image_width, image_height);
	double aspect_ratio = (double)image_height/(double)image_width;
	int samples1D = 70;
	double samples = samples1D*samples1D;
	// int samples1D = sqrt(samples);

	//calculate image plane size in world space
	double scene_width = tan(fov/2)*2;
	double scene_height = tan((fov*aspect_ratio)/2)*2;
	double pixel_width = scene_width/image_width;
	double pixel_slice = pixel_width/samples1D;

	//start timer
	clock_t start = clock();

	int total_pixels = image_width*image_height;
	double last_report = 0;

	//for each pixel
	for(int i=0; i<image_width; i++)
	// for(int i=image_width-1; i>=0; i--)
	{
		for(int j=0; j<image_height; j++)
		{
			int ii = i - image_width/2;
			int jj = j - image_height/2;
			double x = ii*pixel_width + pixel_slice/2;
			double y = jj*pixel_width + pixel_slice/2;
			
			vec4 color(0,0,0);
			//for each sample
			for(int m=0; m<samples1D; m++)
			{
				for(int n=0; n<samples1D; n++)
				{
					double xoff = r.randd_negative()*pixel_slice;
					double yoff = r.randd_negative()*pixel_slice;
					vec4 dir(x+(pixel_slice*m)+xoff,y+(pixel_slice*n)+yoff,0);
					dir -= look_from;
					dir.normalize();
					ray v(look_from, dir);
v.debug = 0;
// if(i==248 &&j==(image_height-155)&&m==0&&n==0) 
// {
// 	v.debug=1;
// 	std::cout << "debugging" << std::endl;
// }
// if(i==70&&j==(image_height -169)&&m==0&&n==0)
// {
// 	v.debug = 1;
// std::cout << " debugging" << std::endl;
// }
// else
// {
// 	v.debug = 0;
// }
// std::cout << "before trace" << std::endl;
					color += r.trace_path(v, 0, -1);
// std::cout << "after trace" << std::endl;
				}
			}
			//TODO: change back to 1!!!
			color *= 1/samples;
			color.clamp(1.0);
// std::cout << i << " " << j << std::endl;
// if(i==1100&&j==600)
// {
// std::cout <<"color: " << color.x << " " << color.y << " " << color.z << std::endl;
// }
			pic.setpixel(i, j, color);
			double decimal = ((double)((i+1)*(j+1)))/((double)(total_pixels));
			// double decimal = ((double)((image_width - i)*(j+1)))/((double)(total_pixels));

			if(decimal-last_report > 0.1)
			{
				last_report = decimal;
				std::cout << (int)(decimal*100.0) << "\%" << std::endl;
				pic.writeppm("clock.ppm");
			}
			
		}
	}
	pic.writeppm("clock.ppm");

	clock_t time_gone_by = clock() - start;
	std::cout << "TOTAL TIME: " << (double)time_gone_by / ((double)CLOCKS_PER_SEC) << " seconds" << std::endl;
	
}