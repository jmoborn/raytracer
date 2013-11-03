#ifndef __MESH_H__
#define __MESH_H__

#include <iostream>
#include <stdio.h>
#include <sstream>
#include <stdlib.h>
#include <string>
#include <vector>
#include "vec4.h"
#include "object.h"

class face {
public:
	face(std::vector<int> pts, std::vector<int> txts);
	~face(){};

	std::vector<int> pnts;
	std::vector<int> txts;

private:
	bool isQuad;
};

class mesh : public object {
public:
	mesh();
	mesh(std::string filepath);
	~mesh(){};
	bool intersect(ray& r);
	vec4 get_normal(const vec4& p);
	vec4 get_color();

	std::vector<vec4> verts;
	std::vector<face> faces;
	std::vector<vec2> texts;

private:
	void readobj(std::string& filepath);

};

#endif