#include "raytracer.h"

raytracer::raytracer()
{
	v = vec4(1,0,0);
	m.rotateY(60.0);
	v *= m;
	g = mesh();
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
		color += v.hit_color * lights[l]->get_color() * std::max(0.0, n_dot_l);
		color += lights[l]->get_color() * pow(std::max(0.0, eye.dot(reflect)), 32.0);
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
	sphere s(0.5, vec4(-.5,-.5,-4), vec4(0,0,1));
	sphere s1(0.33, vec4(.014, -.6, -3.775), vec4(1,0,0));
	sphere tiny(0.2, vec4(.5,.5, -4), vec4(0,1,0));
	sphere tinier(0.11, vec4(.362, .459, -3.418), vec4(1,1,0));
	r.objs.push_back(&s);
	r.objs.push_back(&s1);
	r.objs.push_back(&tiny);
	r.objs.push_back(&tinier);

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
					int closest_obj = -1;
					double min_dist = std::numeric_limits<double>::infinity();
					//for each object
					for(int o=0; o<r.objs.size(); o++)
					{
						
						if(r.objs[o]->intersect(v))
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
						color += r.shade(v);
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
					}
				}
			}
			color *= 1/samples;
			color.clamp(1.0);
			pic.setpixel(i, j, color);
		}
	}
	pic.writeppm("test.ppm");
}