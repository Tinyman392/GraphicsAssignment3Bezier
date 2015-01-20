#include <iostream>
#include <cmath>
#include <stdlib.h>
#include <vector>
using namespace std;

//pi
#define PI 3.141592654
//for testing
string interpolate = "on";
const int u = 5;
const int v = 3;
const int divisions=10;
const int zBnd = 10;
const int xFact = 150;
const int numFrames = 1000;
const float ka = 0.25;
const float kd = 0.75;
const float ks = 0.75;
const string surface = "\"painted plastic\"";
const float brightness = 512;

/*
const char * header [] = 
  {
    "Display \"Animation\"  \"Screen\"  \"rgbdouble\"",
    "Format 640 480",
    "CameraAt 0 0 0", 
    "CameraUp 0 0 1",
    "CameraEye 0 0 40",
    "CameraFOV 50" 
  };
*/

const char * header [] = 
{
    "Display \"Animation\"  \"Screen\"  \"rgbdouble\"\n",
    "Format 640 480",
    "CameraAt 0 0 0", 
    "CameraUp 0 0 1",
    "CameraEye 75 75 25",
    "CameraFOV 20\n" 
};

//get the size of the header array
int nlines = sizeof(header) / sizeof(char *);
 
struct coord
{
	float x;
 	float y;
 	float z;
 	bool down;
};
 
void initPts(vector<coord> &pts);
void addPt(vector<coord> &pts, int j, int i);

int main()
{
	//print header
	for(int i = 0; i < nlines; i++)
	{
		cout << header[i] << endl;
	}
	
	vector<coord> pts;
	
	//initial frame and balls vector
	int frame = 0;
	
	initPts(pts);
	coord c;
	c.x = (rand() % 50 + 50) / 100.0;
	c.y = (rand() % 50 + 50) / 100.0;
	c.z = (rand() % 50 + 50) / 100.0;
	bool downr = rand() % 2;
	bool downg = rand() % 2;
	bool downb = rand() % 2;
 	
 	//for(frame = 0; frame < numFrames; frame++)
	for(frame = 0; ; frame++)
	{
		//start frame and world blocks
		cout << "FrameBegin " << frame << endl;
		cout << "WorldBegin" << endl << endl;
		
		cout << "OptionReal \"Divisions\" " << divisions << endl;
		cout << "OptionBool \"Interpolate\" " << interpolate << endl;
		
		//cout << "PointLight 0 0 0 1 1 1 " << brightness << endl;
		cout << "PointLight -15 15 40 1 1 1 " << brightness << endl;
		cout << "PointLight 15 15 40 1 1 1 " << brightness << endl;
		cout << "PointLight 15 -15 40 1 1 1 " << brightness << endl;
		cout << "PointLight -15 -15 40 1 1 1 " << brightness << "\n\n";
		
		cout << "Surface " << surface << endl;
		
		cout << "Ka " << ka << endl;
		cout << "Kd " << kd << endl;
		cout << "Ks " << ks << endl;
		
		cout << "Color " << c.x << " " << c.y << " " << c.z << endl;
		
		cout << "Patch \"Bezier\" \"P\"\n";
		cout << u*2 << " " << v*2 << endl;
		
		for(int i = 0; i < (u*2+1)*(v*2+1); i++)
		{
			cout << pts[i].x << " " << pts[i].y << " " << pts[i].z << endl;
			
			int upDown = rand() % 2;
			float inc1 = (rand() % 500) * (1.0/xFact);
			float inc2 = (rand() % 500) * (1.0/xFact);
			float inc3 = (rand() % 500) * (1.0/xFact);
			
			float incr = (rand() % 100) / 100000.0;
			float incg = (rand() % 100) / 100000.0;
			float incb = (rand() % 100) / 100000.0;
			
			if(downr)
				c.x += incr;
			if(!downr)
				c.x -= incr;
			
			if(downg)
				c.y += incg;
			if(!downg)
				c.y -= incg;
			
			if(downb)
				c.z += incb;
			if(!downb)
				c.z -= incb;
			
			//cerr << inc << endl;
			
			if(pts[i].down)
			{
				//pts[i].x += inc1;
				//pts[i].y += inc2;
				pts[i].z += inc3;
			}
			else
			{
				//pts[i].x -= inc1;
				//pts[i].y -= inc2;
				pts[i].z -= inc3;
			}
			
			/*
			if(pts[i].x > -10)
				pts[i].x -= 2*inc1;
			if(pts[i].x < 10)
				pts[i].x += 2*inc1;
			
			if(pts[i].y > -10)
				pts[i].y -= 2*inc2;
			if(pts[i].y < 10)
				pts[i].y += 2*inc2;
			*/
			
			if(pts[i].z < -1*zBnd)
			{
				pts[i].down = !pts[i].down;
				pts[i].z += 2*inc3;
			}
			if(pts[i].z > zBnd)
			{
				pts[i].down = !pts[i].down;
				pts[i].z -= 2*inc3;
			}
			
			if(c.x > 1)
			{
				downr = !downr;
				c.x -= 2* incr;
			}
			else if(c.x < .25)
			{
				downr = !downr;
				c.x += 2*incr;
			}
			
			if(c.y > 1)
			{
				downg = !downg;
				c.y -= 2* incg;
			}
			else if(c.y < .25)
			{
				downg = !downg;
				c.y += 2*incg;
			}
			
			if(c.z > 1)
			{
				downb = !downb;
				c.z -= 2* incb;
			}
			else if(c.z < .25)
			{
				downb = !downb;
				c.z += 2*incb;
			}
		}
		
		//end frame and world blocks
		cout << "WorldEnd" << endl;
		cout << "FrameEnd" << endl << endl;
	}
	
	return 0;
}

void initPts(vector<coord> &pts)
{
	for(int i = -1*v; i <= v; i++)
	{
		for(int j = -1*u; j <= u; j++)
		{
			addPt(pts, j, i);
		}
	}
}

void addPt(vector<coord> &pts, int j, int i)
{
	coord c;
	
	c.x = 3*j;
	c.y = 3*i;
	c.z = 0;
	
	c.down = rand() % 2;
	
	pts.push_back(c);
}