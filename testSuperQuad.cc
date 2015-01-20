#include <iostream>
#include <cmath>
#include <stdlib.h>
#include <vector>
using namespace std;

#define PI 3.141592654
#define numQuads 40
#define numFrames 1000
#define windSpd 0.5

const char * header [] = 
{
    "Display \"Animation\"  \"Screen\"  \"rgbdouble\"",
    "Format 640 480",
    "CameraAt 0 0 10", 
    "CameraUp 0 0 1",
    "CameraEye 40 40 40",
    "CameraFOV 50" 
 };

int nlines = sizeof(header) / sizeof(char *);

struct coord
{
	float x;
	float y;
	float z;
};

struct quad
{
	float R;
	float r;
	float n;
	float e;
	float dz;
	coord pos;
	bool isS;
	coord color;
	int surface;
	coord vel;
	float angle1;
	float angle2;
	/*
		phimin = -PI/2
		phimax = PI/2
		theta/thetaMax = 360
	*/
};

quad initQuad();
float checkNorthEast(float ne);
coord getWind(coord c);
quad advFrame(quad q);
bool checkBounds(coord c);
void drawQuad(quad q);
void drawSphere(quad q);
void drawTorus(quad q);
void getSurface(int s);

int main (int argc, char * argv[])
{

  	// Print out header
  
  	for (int i = 0; i < nlines; i++)
    	cout << header[i] << endl;

  	// Loop over frames
	
	vector<quad> quads;
	
	for(int i = 0; i < numQuads; i++)
	{
		quads.push_back(initQuad());
	}
	
  	int frame = 0;

  	for(frame = 0; ; frame++)
	{
      	cout << "FrameBegin " << frame << endl;
      	cout << "WorldBegin" << endl << endl;
      	
      	cout << "\tOptionReal \"Divisions\" 20\n";
      	cout << "\tPointLight 0 0 -3 1 1 1 64\n";
      	cout << "\tPointLight 0 0 40 1 1 1 512\n";
    	cout << "\tPointLight 40 0 40 1 1 1 512\n";
    	cout << "\tPointLight 0 40 40 1 1 1 512\n";
    	
    	cout << "\tKa 0.25\n";
    	cout << "\tKd 0.75\n";
    	cout << "\tKs 0.75\n\n";
      
      	for(int i = 0; i < numQuads; i++)
      	{
      		cout << endl;
      		drawQuad(quads[i]);
      		quads[i] = advFrame(quads[i]);
      	}
		
		cout << endl;
      	cout << "WorldEnd" << endl;
      	cout << "FrameEnd" << endl;
    }
}

quad initQuad()
{
	quad q;
	q.R = (rand() % 200 / 100.0) + 0.75;
	q.r = (rand() % 75 / 100.0) + 0.50;
	q.n = rand() % 300 / 100.0;
	q.e = rand() % 300 / 100.0;
	//checkNorthEast(q.n);
	//checkNorthEast(q.e);
	
	q.pos.x = (rand() % 200 / 100.0) - 1;
	q.pos.y = (rand() % 200 / 100.0) - 1;
	q.pos.z = 0;
	
	q.dz = (rand() % 450 / 1000.0) + 0.05;
	
	q.isS = rand() % 2;
	
	q.surface = rand() % 4;
	
	q.color.x = rand() % 100 / 100.0;
	q.color.y = rand() % 100 / 100.0;
	q.color.z = rand() % 100 / 100.0;
	
	q.vel.x = rand() % 25 / 100.0 - 0.125;
	q.vel.y = rand() % 25 / 100.0 - 0.125;
	
	q.angle1 = rand() % 600 / 100.0;
	q.angle2 = rand() % 600 / 100.0;
	
	return q;
}

float checkNorthEast(float ne)
{
	while(ne == 0)
		ne = rand() % 200 / 100.0;
	
	return ne;
}

coord getWind(coord c)
{
	coord w;
	w.x = windSpd*sin(c.y);
	w.y = windSpd*sin(c.x);
	w.z = c.z;
	
	return w;
}

quad advFrame(quad q)
{
	coord c = q.pos;
	coord w = getWind(c);
	
	c.x += q.vel.x + w.x;
	c.y += q.vel.y + w.y;
	c.z += q.dz;
	
	q.pos = c;
	
	if(checkBounds(q.pos))
		q = initQuad();
	
	return q;
}

bool checkBounds(coord c)
{
	if(abs(c.x) > 30)
		return true;
	if(abs(c.y) > 30)
		return true;
	if(abs(c.z) > 30)
		return true;
	
	return false;
}

void drawQuad(quad q)
{
	cout << "\tColor " << q.color.x << " " << q.color.y << " ";
		cout << q.color.z << endl;
	getSurface(q.surface);
	cout << "\tXformPush\n";
	cout << "\t\tRotate " << "\"X\" " << q.angle1 << endl;
	cout << "\t\tRotate " << "\"Y\" " << q.angle2 << endl;
	cout << "\t\tTranslate " << q.pos.x << " " << q.pos.y << " ";
		cout << q.pos.z << endl;
	
	if(q.isS)
		drawSphere(q);
	else
		drawTorus(q);
	
	cout << "\tXformPop\n";
}

void drawSphere(quad q)
{
	cout << "\t\tSqSphere " << q.R << " " << q.n << " " << q.e << " ";
		cout << q.R*(-1) << " " << q.R << " 360\n";
}

void drawTorus(quad q)
{
	cout << "\t\tSqTorus " << q.R << " " << q.r << " " << q.n << " ";
		cout << q.e << " " << -PI/2 << " " << PI/2 << " 360\n";
}

void getSurface(int s)
{
	string surface = "\"matte\"";
	
	if(s == 0)
		surface = "\"matte\"";
	else if(s == 1)
		surface = "\"metal\"";
	else if(s == 2)
		surface = "\"plastic\"";
	else if(s == 3)
		surface = "\"painted plastic\"";
	
	cout << "\tSurface " << surface << endl;
}