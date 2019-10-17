#include <iostream>
#include <vector> 
#include <stdlib.h>
#include <cassert>
#include <iterator>

using namespace std;

// structs for testing:
struct point {
    int x, y;
};


typedef vector<point>::iterator I; //Define I type as I


// Standard Functions: 

bool isLeftTurn ( point p1, point p2, point p3){
float test = (p3.y - p1.y)*(p2.x - p1.x) -
             (p2.y - p1.y)*(p3.x - p1.x);
return (test > 0 ) ? 1 : 0; 
}

//Equality operator for our points
bool operator==(const point& p1, const point& p2) {
    return (p1.x == p2.x && p1.y == p2.y);
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
        I iter;
        for (iter = this->ConveXHull.begin(); iter != this->ConveXHull.end(); iter++){
            std::cout << (*iter).x << " ";
            std::cout << (*iter).y << " \n";
        }
        return 0;
    }
};

bool noTurn(point p, point q, point r) {
    int px = p.x;
    int py = p.y;
    int qx = q.x;
    int qy = q.y;
    int rx = r.x;
    int ry = r.y;
    int lhs = (qx - px) * (ry - py);
    int rhs = (rx - px) * (qy - py);

    return lhs == rhs;
}
//Find if point p is on line segment qr
bool onLineSegment(point p, point q, point r) {
    point left = (q.x < r.x) ? q : r;
    point right = (q.x > r.x) ? q : r;
    if (p.x < left.x) {
        return false;
    }
    if (p.x > right.x) {
        return false;
    }
    return noTurn(p, q, r);
}

// Function to find eight extreme points
vector<point> findExtrema(I first, I past) {
    assert(first != past);
    vector<point> max_position(8, (*first)); //Vector to keep extreme points 
    
    //Define starting values for the eight extreme points
    int xmin = (*first).x; 
    int xmax = (*first).x;
    int ymin = (*first).y;
    int ymax = (*first).y;
    int ne = xmin + ymin;
    int se = xmin - ymin;
    int sw = -(xmin + ymin);
    int nw = ymin - xmin;

    I i; 
    for (i = first; i != past; ++i) { //Iterate over each point in the convex hull
        if ((*i).x < xmin) { //Update min x value 
            xmin = (*i).x; 
            max_position[0] = (*i);
        }
        if ((*i).x > xmax) { //Update max x value
            xmax = (*i).x;
            max_position[1] = (*i);
        }
        if ((*i).y < ymin) { //Update min y value
            ymin = (*i).y;
            max_position[2] = (*i);
        }
        if ((*i).y > ymax) { //Update max y value
            ymax = (*i).y;
            max_position[3] = (*i);    
        }
        int d = (*i).x + (*i).y; //Check if north-eastern point
        
        if (d > ne) {
            ne = d; //Update north-eastern point
            max_position[4] = (*i); 
        }
        d = (*i).x - (*i).y; //Check if south-eastern point
        if (d > se) { //Update south-eastern point
            se = d;
            max_position[5] = (*i);
        }
        d = -((*i).x + (*i).y); //Check if south-western point 
        if (d > sw) { //Update south-western point
            sw = d;
            max_position[6] = (*i);
        }
        d = (*i).y - (*i).x; //Chck if north-western point
        if (d > nw) { //Update north-western point
            nw = d;
            max_position[7] = (*i);
        }
    }
    return max_position;

}   
//Remove duplicates from set of points
void removeDuplicates(vector<point> S) {
    assert(S.size() != 0);
    size_t i = 0; 
    size_t j = 1;
    point v = S[0];
    while (j != S.size()) {
        if (not(S[j] == v)) {
            std::swap(S[i], v);
            ++i;
            v = S[j];
        }
        ++j;
    }
    std::swap(S[i], v);
    ++i;
    std:size_t remains = S.size() - i;
    while (remains != 0) {
        S.pop_back();
        --remains;
    }
}

bool betweenPolynomialChains(point p, I upper, I uend, I lower, I lend) {
    I i = upper;
    I last = std::prev(uend);
    do {
        ++i;
    } while (i != last and (*i).x < p.x);
    if (isLeftTurn(*std::prev(i), *i, p)) {
        return false;
    }
    last = std::prev(lend);
    I j = lower;
    do {
        ++j;
    } while (j != last and (*j).x > p.x);
    if (isLeftTurn(*std::prev(j), *j, p)) {
        return false;
    }
    return true;
}

I eliminateInnerPoints(I first, I past, I polygon, I rest) {
    assert(polygon != rest);
    using P = point;
    if (std::distance(polygon, rest) == 1) {
        I mid = std::partition(first, past, [&](P point r) -> bool {
            return not (r == *polygon);
        });
        return mid;
    }
    
}
// The general idea is, we create classes which inherit the baseclass 
// and have a additional function, which implements different throwaway tactics.
struct testPlane : PointPlane {
    int ThrowAwayHull(){
        std::cout  << " DETTE ER EN TEST FOR SATAN! ";
    }
};

int main() 
{   
    
    PointPlane plane;
    plane.GeneratePointList(20,100,10);
    plane.GetHullPoints();
    testPlane testo;
    testo.ThrowAwayHull();
    return 0;
}
