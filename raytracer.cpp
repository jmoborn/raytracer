#include "raytracer.h"

raytracer::raytracer()
{
	max_depth = 4;
}

raytracer::~raytracer()
{
	
}

vec4 raytracer::shade(ray& v)
{
	vec4 color;
	for(int l=0; l<lights.size(); l++)
	{
		vec4 pt = v.end();
		vec4 light = lights[l]->c - pt;
		light.normalize();
		v.hit_norm.normalize();
		double n_dot_l = v.hit_norm.dot(light);
		vec4 reflect = v.hit_norm*(2*n_dot_l) - light;
		vec4 eye = v.d*(-1);
		reflect.normalize();
		color += v.hit_color * lights[l]->diffuse() * std::max(0.0, n_dot_l);
		color += lights[l]->diffuse() * pow(std::max(0.0, eye.dot(reflect)), 32.0);
	}
	return color;
}

vec4 raytracer::trace_ray(ray& v, int depth)
{
	int closest_obj = -1;
	double min_dist = std::numeric_limits<double>::infinity();
	vec4 color;
	//for each object
	for(int o=0; o<objs.size(); o++)
	{
		// if(o==5) std::cout << v.hit_color.x << std::endl;
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
		if(closest_obj==5) std::cout << shade(v).x << std::endl;
		color += shade(v);
		// std::cout << "depth: " << depth << std::endl;
		// std::cout << "max depth: " << max_depth << std::endl;
		// std::cout << objs[closest_obj]->reflect().length() << std::endl;
		if((objs[closest_obj]->reflect().length())&&(depth<max_depth))
		{
			//calculate reflection ray
			vec4 org = v.end();
			vec4 inc = v.d*(-1);
			vec4 dir = v.hit_norm*(v.hit_norm.dot(inc)*2) - inc;
			ray r(org, dir);
			color += objs[closest_obj]->reflect()*trace_ray(r, depth+1);
			//calculate refraction ray
			// double iof = 1.5;
			// double c1 = -v.hit_norm.dot(inc);
			// double c2 = sqrt(1-iof*iof*(1-c1*c1));
			// vec4 fract = inc*iof + v.hit_norm * (iof*c1 - c2);
			// ray f(org, fract);
			// color += trace_ray(f, depth+1)*0.5;
		}
	}
	return color;
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

	//create objects
	sphere s(0.5, vec4(-.5,-.5,-4), vec4(0,0,1), 0.25);
	// sphere s1(0.33, vec4(.014, -.6, -3.775), vec4(1,0,0), 0.25);
	sphere s1(0.33, vec4(.014, -.6, -3.775), vec4(1,1,1), 0.25);
	sphere tiny(0.2, vec4(.5,.5, -4), vec4(0,1,0), 0.25);
	sphere tinier(0.11, vec4(.362, .459, -3.418), vec4(1,1,0), 0.25);
	sphere ref(0.2, vec4(-1.27, -.913, -3.254), vec4(1,0,0), 0.25);
	mesh cube("cube.obj", vec4(1.0, 0, 1.0));
	r.objs.push_back(&s);
	r.objs.push_back(&s1);
	r.objs.push_back(&tiny);
	r.objs.push_back(&tinier);
	r.objs.push_back(&ref);
	// r.objs.push_back(&cube);

	//create lights
	// vec4 l(0,0,-3.5);
	sphere l1(1.0, vec4(0,0,-3.5), vec4(1.0,1.0,1.0));
	sphere l2(0.5, vec4(-1,-.75,-3),  vec4(.3,.3,.3));
	r.lights.push_back(&l1);
	r.lights.push_back(&l2);
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
					vec4 dir(x+(pixel_slice*m),y+(pixel_slice*n),0);
					dir -= look_from;
					// std::cout << i << ", " << j << ": " << dir.x << " " << dir.y << " " << dir.z << std::endl;
					dir.normalize();
					ray v(look_from, dir);
					color += r.trace_ray(v, 0);
					// if(cube.intersect(v)) color += vec4(1,1,1);
					// int closest_obj = -1;
					// double min_dist = std::numeric_limits<double>::infinity();
					// //for each object
					// for(int o=0; o<r.objs.size(); o++)
					// {
						
					// 	if(r.objs[o]->intersect(v))
					// 	{
					// 		if(v.t<min_dist)
					// 		{
					// 			min_dist = v.t;
					// 			closest_obj = o;
					// 		}
					// 	}
					// }
					// if(closest_obj!=-1)
					// {
					// 	color += r.shade(v);
					// 	if(r.objs[closest_obj]->reflect)
					// 	{

					// 	}
						// for(int l=0; l<r.lights.size(); l++)
						// {
						// 	// std::cout << "we had an intersection" << std::endl;
						// 	vec4 pt = v.end();
						// 	// vec4 normal = r.objs[closest_obj]->get_normal(pt);
						// 	vec4 light = r.lights[l]->c - pt;
						// 	// std::cout << pt.x << " " << pt.y << " " << pt.z << std::endl;
						// 	light.normalize();
						// 	v.hit_norm.normalize();
						// 	double n_dot_l = v.hit_norm.dot(light);
						// 	vec4 reflect = v.hit_norm*(2*n_dot_l) - light;
						// 	vec4 eye = v.d*(-1);
						// 	reflect.normalize();
						// 	color += v.hit_color * r.lights[l]->get_color() * std::max(0.0, n_dot_l);
						// 	color += r.lights[l]->get_color() * pow(std::max(0.0, eye.dot(reflect)), 32.0);
						// }
					// }
				}
			}
			color *= 1/samples;
			color.clamp(1.0);
			pic.setpixel(i, j, color);
		}
	}
	pic.writeppm("test.ppm");
}