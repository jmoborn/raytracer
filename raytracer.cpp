#include "raytracer.h"

raytracer::raytracer()
{
	srand(time(NULL));
	max_depth = 3;
	shadow_samples = 32;
	samples = 4;
	ambience = 0.1;
	ray_tolerance = 0.001;
	std::string scenefile = "example.scene";
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
			sphere *l = new sphere(atof(sradius.c_str()), l_pos);
			std::cout << "radius " << l->r << std::endl;
			std::cout << "position " << l->c.x << " " << l->c.y << " " << l->c.z << std::endl;
			vec4 tcol = l->diffuse();
			std::cout << "diffuse " << tcol.x << " " << tcol.y << " " << tcol.z << std::endl;
			this->lights.push_back(l);
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
		vec4 light = lights[l]->c - pt;
		light.normalize();
		pt += (light*ray_tolerance);
		ray s(pt, light);
		double shadow_mult = trace_shadow(s, self, l);
		double n_dot_l = v.hit_norm.dot(light);
		vec4 reflect = v.hit_norm*(2*n_dot_l) - light;
		vec4 eye = v.d*(-1);
		reflect.normalize();
		color += v.hit_color*ambience + v.hit_color * lights[l]->diffuse() * std::max(0.0, n_dot_l);
		double phong = objs[self]->shader.specular;
		if(phong!=0)
			color  += lights[l]->diffuse() * pow(std::max(0.0, eye.dot(reflect)), phong);
		color *= shadow_mult;
		if(depth==0)
			color += ambience;
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
				r.o = reflect_org;
				r.d = dir;
				refl_color = objs[closest_obj]->reflect()*trace_ray(r, depth+1, closest_obj);
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
		color += (.2, .2, .2);
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
	ray vl = v;
	this->lights[light]->intersect(vl);
	int light_dist = vl.t;
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
		vs.d = new_dir;
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

int main()
{
	raytracer r;
	
	//TODO: read this info from the scene file
	//define camera properties
	int image_width = 512;
	int image_height = 512;
	vec4 look_from(0,0,1);
	double fov = 45.0;//degrees
	fov *= (r.PI/180.0);//radians
	pixelmap pic(image_width, image_height);
	double aspect_ratio = (double)image_height/(double)image_width;
	double samples = 4;
	int samples1D = sqrt(samples);

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
					color += r.trace_ray(v, 0, -1);
				}
			}
			color *= 1/samples;
			color.clamp(1.0);
			pic.setpixel(i, j, color);
			double decimal = ((double)((i+1)*(j+1)))/((double)(total_pixels));
			if(decimal-last_report > 0.1)
			{
				last_report = decimal;
				std::cout << (int)(decimal*100.0) << "\%" << std::endl;
			}
			
		}
	}
	pic.writeppm("test.ppm");

	clock_t time_gone_by = clock() - start;
	std::cout << "TOTAL TIME: " << (double)time_gone_by / ((double)CLOCKS_PER_SEC) << " seconds" << std::endl;
}