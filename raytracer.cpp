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
	// sphere sp(1.5, vec4(0,0,-1));
	// ray test1(vec4(0,0,1), vec4(0,0,-1));
	// if(sp.intersect(test1))
	// 	std::cout << "-2: " << test1.t << std::endl;
	// ray test2(vec4(-4,0,0), vec4(1,0,0));
	// if(s.intersect(test2))
	// 	std::cout << "-4: " <<test2.t << std::endl;
	// ray test3(vec4(-4,0,0), vec4(0,0,1));
	// if(s.intersect(test3))
	// 	std::cout << "4: " <<test3.t << std::endl;
	
	int image_width = 512;
	int image_height = 512;
	
	vec4 look_from(0,0,1);
	double fov = 45.0;

	pixelmap pic(image_width, image_height);
	double aspect_ratio = image_height/image_width;
	double samples = 16;
	int samples1D = sqrt(samples);

	double scene_width = tan(fov);
	double scene_height = tan(fov*aspect_ratio);
	double pixel_width = scene_width/image_width;
	double pixel_slice = pixel_width/samples1D;

	sphere s(0.5, vec4(0,0,-1));
	sphere tiny(0.1, vec4(0,.75, -1));

	vec4 l(1,-1,1);
	l.normalize();

	// ray test1(vec4(-1,0,0), vec4(0.99813, -0.06111, -0.06111));
	// if(s.intersect(test1))
	// 	std::cout << test1.t << std::endl;

	for(int i=0; i<image_width; i++)
	{
		for(int j=0; j<image_height; j++)
		{
			int ii = i - image_width/2;
			int jj = j - image_height/2;
			double x = ii*pixel_width + pixel_slice/2;
			double y = jj*pixel_width + pixel_slice/2;
			
			vec4 color(0,0,0);
			for(int m=0; m<samples1D; m++)
			{
				for(int n=0; n<samples1D; n++)
				{
					vec4 org(0,0,1);
					vec4 dir(x+(pixel_slice*m),y+(pixel_slice*n),0);
					dir -= look_from;
					dir.normalize();
					ray v(look_from, dir);
					if(s.intersect(v))
					{
						vec4 normal = s.normal(v.end());
						normal.normalize();
						color += vec4(0,0,1) * std::max(0.0, normal.dot(l));
					}
					if(tiny.intersect(v))
					{
						vec4 normal = s.normal(v.end());
						normal.normalize();
						// if(normal.dot(l)>0) std::cout << normal.dot(l) << std::endl;
						color += vec4(0,1,0) * std::max(0.0, normal.dot(l));
					}
				}
			}
			color *= 1/samples;
			color += .1;
			color.clamp(1.0);
			pic.setpixel(i, j, color);
			// vec4 org(0,0,1);
			// vec4 dir(x,y,0);
			// dir -= org;
			// dir.normalize();
			// ray v(org, dir);
			// if(s.intersect(v))
			// {
			// 	pic.setpixel(i, j, vec4(0,0,1));
			// }
			// if(tiny.intersect(v))
			// {
			// 	pic.setpixel(i, j, vec4(0,1,0));
			// }
		}
	}
	pic.writeppm("test.ppm");
}