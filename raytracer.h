#include <string>
#include <iostream>
#include <vector>
#include <math.h>

#include "mat4.h"
#include "vec4.h"
#include "sphere.h"
#include "mesh.h"
#include "pixelmap.h"

class raytracer {
public:
	raytracer();
	~raytracer();
private:
	vec4 v;
	mat4 m;
	mesh g;
};