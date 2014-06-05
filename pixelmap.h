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
	pixelmap(std::string filename);
	pixelmap(const pixelmap &copy);
	pixelmap();
	~pixelmap();

	void setpixel(int x, int y, vec4 c);
	vec4 getpixel(int x, int y);

	//returns the bilinearly interpolated color
	vec4 getpixel(double u, double v);
	
	//pixels should be set from bottom to top, left to right
	//however ppms are left to right, top to bottom
	//this function takes care of flipping them
	void writeppm(std::string filename);
	
	int width;
	int height;

	pixelmap& operator=(const pixelmap other);

private:
	char skipcomment(std::ifstream &file);
	void readppm(std::string filename);
	vec4 **pixels;

};

#endif