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
	void writeppm(std::string filename);
	
	int width;
	int height;

private:
	vec4 **pixels;

};

#endif