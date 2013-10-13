#include <string>
#include <iostream>

#include "mat4.h"
#include "vec4.h"
#include "mesh.h"

class raytracer {
public:
	raytracer();
	~raytracer();
private:
	vec4 v;
	mat4 m;
	mesh g;
};