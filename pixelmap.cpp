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

pixelmap::pixelmap(std::string filename)
{
	readppm(filename);
}

pixelmap::~pixelmap()
{
	for(int i=0; i<width; i++)
	{
		delete [] pixels[i];
	}
	delete [] this->pixels;
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
  	for(int i=height-1; i>=0; i--)
  	{
  		for(int j=0; j<width; j++)
  		{
  			file << round(pixels[j][i].x*255) << " " 
  				 << round(pixels[j][i].y*255) << " " 
  				 << round(pixels[j][i].z*255) << "  ";
  		}
  		file << "\n";
  	}
  	file.close();
}

// get next non-whitespace character
char pixelmap::skipcomment(std::ifstream &file)
{
	char next = file.peek();
	while(next==' '||next=='\n')
	{
		file.get();
		next = file.peek();
	}
	char s[256];
	if(next=='#') file.getline(s, 256);
	return next;
}

void pixelmap::readppm(std::string filename)
{
	std::ifstream file(filename.c_str());
    //TODO: comments
    std::string filetype;
    file >> filetype;
    std::cout << "filetype: " << filetype << std::endl;
    if(filetype!="P3") return; //wrong file type
    char next = skipcomment(file);
    file >> this->width;
    std::cout << "width: " << width << std::endl;
    next = skipcomment(file);
    file >> this->height;
    std::cout << "height: " << height << std::endl;
    next = skipcomment(file);
    double max;
    file >> max;
    next = skipcomment(file);
    this->pixels = new vec4*[width];
	for(int i=0; i<width; i++)
	{
		pixels[i] = new vec4[height];
		// for(int j=0; j<height; j++)
		// {

		// 	vec4 tmp;
		// 	file >> tmp.x >> tmp.y >> tmp.z;
		// 	double invmax = 1.0/max;
		// 	tmp *= invmax;
		// 	pixels[i][j] = tmp;
		// }
	}

	for(int i=height-1; i>=0; i--)
	{
		for(int j=0; j<width; j++)
		{
			
			vec4 tmp;
			file >> tmp.x >> tmp.y >> tmp.z;
			double invmax = 1.0/max;
			tmp *= invmax;
			pixels[j][i] = tmp;
		}
	}

	file.close();
}