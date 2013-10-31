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
	sphere s(1.5, vec4(0,0,0));
	ray test1(vec4(-2,0,0), vec4(1,0,0));
	if(s.intersect(test1))
		std::cout << "-2: " << test1.t << std::endl;
	ray test2(vec4(-4,0,0), vec4(1,0,0));
	if(s.intersect(test2))
		std::cout << "-4: " <<test2.t << std::endl;
	ray test3(vec4(-4,0,0), vec4(0,0,1));
	if(s.intersect(test3))
		std::cout << "4: " <<test3.t << std::endl;

}