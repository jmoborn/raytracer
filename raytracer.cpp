#include "raytracer.h"

raytracer::raytracer()
{
	srand(time(NULL));
	max_depth = 3;
	shadow_samples = 16;
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
		if((objs[closest_obj]->shader.reflect > 0)||(objs[closest_obj]->shader.refract > 0)&&(depth<max_depth))
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
			if(objs[closest_obj]->shader.reflect > 0&&(depth<max_depth))
			{
				vec4 dir = v.hit_norm*(v.hit_norm.dot(inc)*2) - inc;
				vec4 reflect_org = org + (dir*ray_tolerance);
				r.o = reflect_org;
				r.d = dir;
				refl_color = objs[closest_obj]->reflect()*trace_ray(r, depth+1, closest_obj);
			}
			//calculate refraction ray
			if(objs[closest_obj]->shader.refract > 0&&(depth<max_depth))
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
	
	//define camera properties
	int image_width = 512;
	int image_height = 512;//384;
	vec4 look_from(0,0,1);
	double fov = 45.0;//degrees
	fov *= (r.PI/180.0);//radians
	//REMOVE LINE BELOW!!!
	// fov = 1.10714872*2;
	pixelmap pic(image_width, image_height);
	double aspect_ratio = (double)image_height/(double)image_width;
	double samples = 16;
	int samples1D = sqrt(samples);

	//calculate image plane size in world space
	double scene_width = tan(fov/2)*2;
	double scene_height = tan((fov*aspect_ratio)/2)*2;
	double pixel_width = scene_width/image_width;
	double pixel_slice = pixel_width/samples1D;

	//********************Original Test Scene********************
	// //create objects
	// sphere s(0.5, vec4(-.5,-.5,-4), vec4(0,0,1), 0.25);
	// sphere s1(0.33, vec4(.014, -.6, -3.775), vec4(1,1,1), 0.25);
	// sphere tiny(0.2, vec4(.5,.5, -4), vec4(0,1,0), 0.25);
	// sphere tinier(0.11, vec4(.362, .459, -3.418), vec4(1,1,0), 0.25);
	// sphere ref(0.2, vec4(-1.27, -.913, -3.254), vec4(1,0,0), 0.25);
	// mesh cube("cube.obj", vec4(1.0, 0, 1.0));
	// r.objs.push_back(&s);
	// r.objs.push_back(&s1);
	// r.objs.push_back(&tiny);
	// r.objs.push_back(&tinier);
	// r.objs.push_back(&ref);
	// r.objs.push_back(&cube);

	// //create lights
	// sphere l1(1.0, vec4(0,0,-3.5), vec4(1.0,1.0,1.0));
	// sphere l2(0.5, vec4(-1,-.75,-3),  vec4(.3,.3,.3));
	// r.lights.push_back(&l1);
	// r.lights.push_back(&l2);

	//********************Scenell********************
	// sphere s1(.2, vec4(0,.1,-1.2), vec4(.1, .1, .1), .25);
	// mesh l_tri("l_tri.obj", vec4(1,1,0));
	// mesh r_tri("r_tri.obj", vec4(0,0,1));
	// r.objs.push_back(&s1);
	// r.objs.push_back(&l_tri);
	// r.objs.push_back(&r_tri);

	// sphere l1(1.0, vec4(0,100,0), vec4(1,1,1));
	// r.lights.push_back(&l1);

	//********************Diffuse********************
	// sphere s1(.05, vec4(.35,0,-.1), vec4(1,1,1), 1.0);
	// sphere s2(.075, vec4(.2,0,-.1), vec4(1,0,0), 1.0);
	// sphere s3(.3, vec4(-.6,0,0), vec4(0,1,0), 1.0);
	// mesh m1("tri1.obj", vec4(0,0,1));
	// mesh m2("tri2.obj", vec4(1,1,0));
	// r.objs.push_back(&s1);
	// r.objs.push_back(&s2);
	// r.objs.push_back(&s3);
	// r.objs.push_back(&m1);
	// r.objs.push_back(&m2);
	// sphere l1(1.0, vec4(100,0,0), vec4(1,1,1));
	// r.lights.push_back(&l1);

	// material white_reflect(vec4(1,1,1), vec4(1,1,1), .15, 1.0, 32.0, 0.0, 1.0);
	// material blue_reflect(vec4(0,0,0.5), vec4(1,1,1), 1.0, 1.0, 32.0, 0.0, 1.0);
	// material orange_reflect(vec4(.75,.25,0), vec4(1,1,1), 1.0, 1.0, 32.0, 0.0, 1.0);
	// material yellow_reflect(vec4(1,1,0), vec4(1,1,1), 1.0, 1.0, 32.0, 0.0, 1.0);
	// material red_diffuse(vec4(.5,0,0), vec4(1,1,1), 1.0, 0.20, 32.0, 0.0, 1.0);
	// material green_diffuse(vec4(0,.5,0), vec4(1,1,1), 1.0, 0.20, 32.0, 0.0, 1.0);
	// material blue_diffuse(vec4(0,0,.5), vec4(1,1,1), 1.0, 0.20, 32.0, 0.0, 1.0);
	// material white_small_reflect(vec4(1,1,1), vec4(1,1,1), 1.0, .5, 32.0, 0.0, 1.0);

	// sphere s1(1.5, vec4(3.441, -.121, -12.337), white_reflect);
	// mesh pyramid("pyramid.obj", orange_reflect);
	// mesh cube("cube.obj", blue_reflect);
	// mesh torus("torus.obj", yellow_reflect);
	// mesh ring1("ring1.obj", white_reflect);
	// mesh test("normal_test.obj", white_reflect);
	// mesh ring2("ring2.obj", white_reflect);
	// mesh ring3("ring3.obj", white_reflect);
	// mesh ring4("ring4.obj", white_reflect);
	// mesh ring5("ring5.obj", white_reflect);
	// mesh ring6("ring6.obj", white_reflect);
	// mesh ring7("ring7.obj", white_reflect);
	// mesh ring8("ring8.obj", white_reflect);
	// mesh floorplane("floor.obj", white_small_reflect);
	// mesh ceiling("ceiling.obj", white_small_reflect);
	// mesh back("back.obj", red_diffuse);
	// mesh left("left.obj", blue_diffuse);
	// mesh right("right.obj", green_diffuse);
	// r.objs.push_back(&s1);
	// r.objs.push_back(&pyramid);
	// r.objs.push_back(&cube);
	// r.objs.push_back(&torus);
	// r.objs.push_back(&ring1);
	// r.objs.push_back(&test);
	// r.objs.push_back(&ring2);
	// r.objs.push_back(&ring3);
	// r.objs.push_back(&ring4);
	// r.objs.push_back(&ring5);
	// r.objs.push_back(&ring6);
	// r.objs.push_back(&ring7);
	// r.objs.push_back(&ring8);
	// r.objs.push_back(&floorplane);
	// r.objs.push_back(&ceiling);
	// r.objs.push_back(&back);
	// r.objs.push_back(&left);
	// r.objs.push_back(&right);

	// sphere l1(1.0, vec4(-1,1,-1));
	// r.lights.push_back(&l1);

	//DEBUG SCENE LOADER ******************************
	// material white_reflect(vec4(1,1,1), vec4(1,1,1), .15, 1.0, 32.0, 0.0, 1.0);
	// sphere s(1.5, vec4(3.441, -.121, -12.337), white_reflect);
	// std::cout << "radius " << s.r << std::endl;
	// std::cout << "position " << s.c.x << " " << s.c.y << " " << s.c.z << std::endl;
	// vec4 tcol = s.diffuse();
	// std::cout << "diffuse " << tcol.x << " " << tcol.y << " " << tcol.z << std::endl;
	// r.objs.push_back(&s);

	// sphere l(1.0, vec4(-1,1,-1));
	// std::cout << "radius " << l.r << std::endl;
	// std::cout << "position " << l.c.x << " " << l.c.y << " " << l.c.z << std::endl;
	// vec4 lcol = l.diffuse();
	// std::cout << "diffuse " << lcol.x << " " << lcol.y << " " << lcol.z << std::endl;
	// r.lights.push_back(&l);

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
	std::cout << "TOTAL TIME: " << (double)time_gone_by / ((double)CLOCKS_PER_SEC) << std::endl;
}