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

	sphere s(0.5, vec4(-.5,-.5,-4), vec4(0,0,1));
	sphere tiny(0.1, vec4(.5,.5, -4), vec4(0,1,0));
	r.objs.push_back(&s);
	r.objs.push_back(&tiny);

	vec4 l(0,0,-3.5);

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
					//for each object
					for(int o=0; o<r.objs.size(); o++)
					{
						if(r.objs[o]->intersect(v))
						{
							vec4 pt = v.end();
							vec4 normal = r.objs[o]->get_normal(pt);
							vec4 light = l - pt;
							// std::cout << pt.x << " " << pt.y << " " << pt.z << std::endl;
							light.normalize();
							normal.normalize();
							color += r.objs[o]->get_color() * std::max(0.0, normal.dot(light));
						}
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