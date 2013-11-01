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
	std::cout << "I don't work yet . . ." << std::endl;
	sphere sp(1.5, vec4(0,0,-1));
	ray test1(vec4(0,0,1), vec4(0,0,-1));
	if(sp.intersect(test1))
		std::cout << "-2: " << test1.t << std::endl;
	// ray test2(vec4(-4,0,0), vec4(1,0,0));
	// if(s.intersect(test2))
	// 	std::cout << "-4: " <<test2.t << std::endl;
	// ray test3(vec4(-4,0,0), vec4(0,0,1));
	// if(s.intersect(test3))
	// 	std::cout << "4: " <<test3.t << std::endl;
	
	int image_width = 512;
	int image_height = 512;
	double fov = 45.0;
	pixelmap pic(image_width, image_height);
	double aspect_ratio = image_height/image_width;
	double scene_width = tan(fov);
	double scene_height = tan(fov*aspect_ratio);
	double pixel_width = scene_width/image_width;

	sphere s(0.5, vec4(0,0,-1));
	sphere tiny(0.1, vec4(0,.75, -1));

	// ray test1(vec4(-1,0,0), vec4(0.99813, -0.06111, -0.06111));
	// if(s.intersect(test1))
	// 	std::cout << test1.t << std::endl;

	for(int i=0; i<image_width; i++)
	{
		for(int j=0; j<image_height; j++)
		{
			int ii = i - image_width/2;
			int jj = j - image_height/2;
			double x = ii*pixel_width + pixel_width/2;
			double y = jj*pixel_width + pixel_width/2;
			vec4 org(0,0,1);
			vec4 dir(x,y,0);
			dir -= org;
			dir.normalize();
			ray v(org, dir);
			// pic.setpixel(i, j, vec4(abs(x/scene_width), abs(y/scene_height), 0));
			if(s.intersect(v))
			{
				pic.setpixel(i, j, vec4(0,0,1));
			}
			if(tiny.intersect(v))
			{
				pic.setpixel(i, j, vec4(0,1,0));
			}
		}
	}
	pic.writeppm("test.ppm");

}