#include "pixelmap.h"

pixelmap::pixelmap(int w, int h)
{
	this->width = w;
	this->height = h;
	this->pixels = new vec4*[width];
	for(int i=0; i<width; i++)
	{
		pixels[i] = new vec4[height];
		for(int j=0; j<height; j++)
		{
			pixels[i][j] = vec4(0,0,0);
		}
	}
}

pixelmap::~pixelmap()
{
	delete this->pixels;
}

void pixelmap::setpixel(int x, int y, vec4 c)
{
	this->pixels[x][y] = c;
}

vec4 pixelmap::getpixel(int x, int y)
{
	return this->pixels[x][y];
}

void pixelmap::writeppm(std::string filename)
{
	std::ofstream file;
  	file.open (filename.c_str());
  	file << "P3\n" << width << " " << height << "\n" << 255 << "\n";
  	for(int i=0; i<width; i++)
  	{
  		for(int j=0; j<height; j++)
  		{
  			file << round(pixels[i][j].x*255) << " " 
  				 << round(pixels[i][j].y*255) << " " 
  				 << round(pixels[i][j].z*255) << "  ";
  		}
  		file << "\n";
  	}
  	file.close();
}