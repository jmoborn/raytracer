#ifndef __IMAGE_H__
#define __IMAGE_H__

#include <iostream>
#include <fstream>
#include <string>
#include <math.h>
#include "vec4.h"

class pixelmap{
public:
	pixelmap(int height, int width);
	~pixelmap();

	void setpixel(int x, int y, vec4 c);
	vec4 getpixel(int x, int y);
	//pixels should be set from bottom to top, left to right
	//however ppms are left to right, top to bottom
	//this function takes care of flipping them
	void writeppm(std::string filename);
	
	int width;
	int height;

private:
	vec4 **pixels;

};

#endif