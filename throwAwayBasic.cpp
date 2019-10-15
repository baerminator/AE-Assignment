#include <iostream>
#include <vector> 
#include <stdlib.h>
using namespace std;


// structs for testing:
struct point {
    int x, y;
};

// Standard Functions: 

bool isLeftTurn ( point p1, point p2, point p3){
float test = (p3.y - p1.y)*(p2.x - p1.x) -
             (p2.y - p1.y)*(p3.x - p1.x);
return (test > 0 ) ? 1 : 0; 
}

// Implement a struct, which contains the convex hull and number of comparisons.
struct PointPlane {
    int NrComps;
    vector<point> ConveXHull;
    vector<point> AllPoints;

    int GeneratePointList (int RangeX,int RangeY, int n){
        vector<point> points;
        for ( int i = 0; i <= n; i++){
            point NewPoint;
            NewPoint.x = rand() % RangeX + 1;
            NewPoint.y = rand() % RangeY + 1;
            points.push_back(NewPoint);
        };
        this->ConveXHull = points;
        this->AllPoints = points; 
        return 0;
    };
    int GetHullPoints (){
        vector<point>::iterator iter;
        for (iter = this->ConveXHull.begin(); iter != this->ConveXHull.end(); iter++){
            std::cout << (*iter).x << " ";
            std::cout << (*iter).y << " \n";
        }
        return 0;
    }
};
// The general idea is, we create classes which inherit the baseclass 
// and have a additional function, which implements different throwaway tactics.
struct testPlane : PointPlane {
    int TestMig(){
        std::cout  << " DETTE ER EN TEST FOR SATAN! ";
    }
};

int main() 
{   
    
    PointPlane plane;
    plane.GeneratePointList(20,100,10);
    plane.GetHullPoints();
    testPlane testo;
    testo.TestMig();
    return 0;
}
