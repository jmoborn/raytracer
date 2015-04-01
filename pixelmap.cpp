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

pixelmap::pixelmap()
{
	this->width = 2;
	this->height = 2;
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

pixelmap::pixelmap(const pixelmap &copy)
{
	// for(int i=0; i<width; i++)
	// {
	// 	delete [] pixels[i];
	// }
	// delete [] this->pixels;
	this->width = copy.width;
	this->height = copy.height;
	this->pixels = new vec4*[width];
	for(int i=0; i<this->width; i++)
	{
		this->pixels[i] = new vec4[this->height];
		for(int j=0; j<height; j++)
		{
			this->pixels[i][j] = copy.pixels[i][j];
		}
	}
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
	this->pixels[x][y].x = c.x;
	this->pixels[x][y].y = c.y;
	this->pixels[x][y].z = c.z;
	// std::cout << "setting pixel: " << x << " " << y << std::endl;
	// std::cout << this->pixels[x][y].x << " " << this->pixels[x][y].y << " " << this->pixels[x][y].z << std::endl;
}

vec4 pixelmap::getpixel(int x, int y)
{
	return this->pixels[x][y];
}

vec4 pixelmap::getpixel_wrap(int x, int y)
{
	x %= width;
	y %= height;
	// if(x>width||x<0) x = x%width;
	// if(y>height||y<0) y = y%height;
	return this->pixels[x][y];

}

vec4 pixelmap::getpixel(double u, double v)
{
	// int i = (u*(double)width);
	// int j = (v*(double)height);
	// return this->pixels[i%width][j%height];

	u = mod(u);
	v = mod(v);
	double x = u*this->width;
	double y = v*this->height;
	double x2 = round(x) + 0.5;
	double x1 = round(x) - 0.5;
	double y2 = round(y) + 0.5;;
	double y1 = round(y) - 0.5;

	vec4 v11 = getpixel_wrap((int)x1,(int)y1);
	vec4 v12 = getpixel_wrap((int)x1,(int)y2);
	vec4 v21 = getpixel_wrap((int)x2,(int)y1);
	vec4 v22 = getpixel_wrap((int)x2,(int)y2);

	return bilerp(x, y, x1, y1, x2, y2, v11, v12, v21, v22);
}

vec4 pixelmap::lerp(double x, double x1, double x2, vec4 v1, vec4 v2)
{
	vec4 w1 = v2*((x - x1)/(x2 - x1));
	vec4 w2 = v1*((x2 - x)/(x2 - x1));
	return w1 + w2;
}

vec4 pixelmap::bilerp(double x, double y, double x1, double y1, double x2, double y2, 
						  vec4 v11, vec4 v12, vec4 v21, vec4 v22)
{
	vec4 wx1 = lerp(x, x1, x2, v11, v21);
	vec4 wx2 = lerp(x, x1, x2, v12, v22);
	vec4 wy  = lerp(y, y1, y2, wx1, wx2);
	return wy;
}

double pixelmap::mod(double u)
{
	double u_floor = floor(u);
	return u - u_floor;
	
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
// if(j==1100&&i==600)
// {
// std::cout << "writing pixel: " << j << " " << i << std::endl;
// std::cout << this->pixels[j][i].x << " " << this->pixels[j][i].y << " " << this->pixels[j][i].z << std::endl;
// }
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
			// if(tmp.x<0.5||tmp.y<0.5||tmp.z<0.5)
			// 	std::cout << tmp.x << ", " << tmp.y << ", " << tmp.z << std::endl;
			pixels[j][i] = tmp;
		}
	}

	file.close();
}

pixelmap& pixelmap::operator=(const pixelmap copy)
{
	// for(int i=0; i<width; i++)
	// {
	// 	delete [] pixels[i];
	// }
	// delete [] this->pixels;
	this->width = copy.width;
	this->height = copy.height;
	this->pixels = new vec4*[width];
	for(int i=0; i<this->width; i++)
	{
		this->pixels[i] = new vec4[this->height];
		for(int j=0; j<height; j++)
		{
			this->pixels[i][j] = copy.pixels[i][j];
		}
	}
}
