#include <algorithm>
#include <iostream>
#include <ostream>
#include <vector>
#include <tuple>
#include <iostream>
#include <math.h>
#include <vector> 
#include <stdlib.h>
#include <limits>
//#include<boost/shared_ptr.hpp>
//#include<CGAL/Polygon_2.h>
//#include<CGAL/create_straight_skeleton_2.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_2.h>

#include <cassert>
#include <iterator>
#include <cmath>
#include <algorithm>
#include <random>
#include <fstream>
using namespace std;
typedef std::tuple<int, int> point;
typedef vector<point>::iterator I;
typedef vector<point> VP;

//###############################################################################################
//########################## SIMPLE PlANESWEEP ################################################
//###############################################################################################

// maybe implement inplace algorithm <3

bool ccw(const point &a, const point &b, const point &c)
{
    return ((std::get<0>(b) - std::get<0>(a)) * (std::get<1>(c) - std::get<1>(a))) > ((std::get<1>(b) - std::get<1>(a)) * (std::get<0>(c) - std::get<0>(a)));
}

//Performs X sort ( maybe implement lexicographical sorting)
std::vector<point> PlaneSweep(std::vector<point> p)
{
    if (p.size() == 0)
        return std::vector<point>();
    std::sort(p.begin(), p.end(), [](point &a, point &b) {
        if (std::get<0>(a) < std::get<0>(b))
            return true;
        return false;
    });

    std::vector<point> RESULT;

    // lower hull
    for (const auto &pt : p)
    {
        while (RESULT.size() >= 2 && !ccw(RESULT.at(RESULT.size() - 2), RESULT.at(RESULT.size() - 1), pt))
        {
            RESULT.pop_back();
        }
        RESULT.push_back(pt);
    }

    // upper hull
    auto t = RESULT.size() + 1;
    for (auto it = p.crbegin(); it != p.crend(); it = std::next(it))
    {
        auto pt = *it;
        while (RESULT.size() >= t && !ccw(RESULT.at(RESULT.size() - 2), RESULT.at(RESULT.size() - 1), pt))
        {
            RESULT.pop_back();
        }
        RESULT.push_back(pt);
    }

    RESULT.pop_back();
    return RESULT;
}

// ############################################################################################
// ######################### Extrema ##########################################################
// ############################################################################################

// Function to find eight extreme points
vector<point> Extrema_8(I first, I past)
{
    assert(first != past);
    vector<point> max_position(8, (*first)); //Vector to keep extreme points

    //Define starting values for the eight extreme points
    int xmin = get<0>(*first);
    int xmax = get<0>(*first);
    int ymin = get<1>(*first);
    int ymax = get<1>(*first);
    int ne = xmin + ymin;
    int se = xmin - ymin;
    int sw = -(xmin + ymin);
    int nw = ymin - xmin;

    I i;
    for (i = first; i != past; ++i)
    { //Iterate over each point in the convex hull
        if (get<0>(*i) < xmin)
        { //Update min x value
            xmin = get<0>(*i);
            max_position[0] = (*i);
        }
        if (get<0>(*i) > xmax)
        { //Update max x value
            xmax = get<0>(*i);
            max_position[4] = (*i);
        }
        if (get<1>(*i) < ymin)
        { //Update min y value
            ymin = get<1>(*i);
            max_position[6] = (*i);
        }
        if (get<1>(*i) > ymax)
        { //Update max y value
            ymax = get<1>(*i);
            max_position[2] = (*i);
        }
        int d = get<0>(*i) + get<1>(*i); //Check if north-eastern point

        if (d > ne)
        {
            ne = d; //Update north-eastern point
            max_position[3] = (*i);
        }
        d = get<0>(*i) - get<1>(*i); //Check if south-eastern point
        if (d > se)
        { //Update south-eastern point
            se = d;
            max_position[5] = (*i);
        }
        d = -(get<0>(*i) + get<1>(*i)); //Check if south-western point
        if (d > sw)
        { //Update south-western point
            sw = d;
            max_position[7] = (*i);
        }
        d = get<1>(*i) - get<0>(*i); //Check if north-western point
        if (d > nw)
        { //Update north-western point
            nw = d;
            max_position[1] = (*i);
        }
    }
    
    std::cout <<"Max Positions are \n";
    for (I iter = max_position.begin(); iter != max_position.end(); iter ++){
        std::cout << get<0>(*iter) << " " << get<1>(*iter) << "\n";
    }
    // This part removes duplicates:
    max_position.resize(std::distance(max_position.begin(), std::unique(max_position.begin(), max_position.end())));

    return max_position;
}

// ############################################################################
// ######################## SQUARE THROW AWAY #################################
// ############################################################################

vector<point> Extrema_4(I first, I past)
{
    assert(first != past);
    vector<point> max_position(4, (*first)); //Vector to keep extreme points

    //Define starting values for the eight extreme points
    int xmin = get<0>(*first);
    int xmax = get<0>(*first);
    int ymin = get<1>(*first);
    int ymax = get<1>(*first);
    int ne = xmin + ymin;
    int se = xmin - ymin;
    int sw = -(xmin + ymin);
    int nw = ymin - xmin;

    I i;
    for (i = first; i != past; ++i)
    { //Iterate over each point in the convex hull
        if (get<0>(*i) < xmin)
        { //Update min x value
            xmin = get<0>(*i);
            max_position[0] = (*i);
        }
        if (get<0>(*i) > xmax)
        { //Update max x value
            xmax = get<0>(*i);
            max_position[2] = (*i);
        }
        if (get<1>(*i) < ymin)
        { //Update min y value
            ymin = get<1>(*i);
            max_position[3] = (*i);
        }
        if (get<1>(*i) > ymax)
        { //Update max y value
            ymax = get<1>(*i);
            max_position[1] = (*i);
        }
    }
    std::cout <<"Max Positions are \n";
    for (I iter = max_position.begin(); iter != max_position.end(); iter ++){
        std::cout << get<0>(*iter) << " " << get<1>(*iter) << "\n";
    }
    // This part removes duplicates:get<0>(*iter)
    max_position.resize(std::distance(max_position.begin(), std::unique(max_position.begin(), max_position.end())));
    return max_position;
}

/*tuple<point, int> InCircle (vector<point> Extrema){
    //typedef CGAL::Exact_predicates_inexact_constructions_kernel K ;
    typedef K::Point_2                   altPoint ;
    //typedef CGAL::Polygon_2<K>           Polygon_2 ;
    //typedef CGAL::Straight_skeleton_2<K> Ss ;
    //typedef boost::shared_ptr<Ss> SsPtr ;
    std::vector< altPoint > points;
    for (I iter = Extrema.begin(); iter != Extrema.end(); iter++){
        points.push_back(altPoint((float)get<0>(*iter), (float)get<1>(*iter)));
    }
    //SsPtr iss = CGAL::create_interior_straight_skeleton_2(points.begin(), points.end());
    //std::cout << iss << std::endl;
    

}
*/


tuple<point, int> InCircle (vector<point> Extrema){
    typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
    typedef CGAL::Delaunay_triangulation_2<K>                   Delaunay;    
    typedef K::Point_2                                          altPoint;
    typedef K::Line_2                                           Line;
    std::vector< altPoint > points;
    //Extrema = Extrema_8(Extrema.begin(),Extrema.end());
    for (I iter = Extrema.begin(); iter != Extrema.end(); iter++){
        points.push_back(altPoint((float)get<0>(*iter), (float)get<1>(*iter)));
    }
    Delaunay dt;
    dt.insert(points.begin(), points.end());
    std::cout << dt.number_of_faces() << std::endl;
    Delaunay::Finite_faces_iterator it;
    altPoint Max;
    float maxRadius = -1;

    for (it = dt.finite_faces_begin(); it != dt.finite_faces_end(); it++) {   
        
        altPoint circum = CGAL::circumcenter(dt.triangle(it)[0],dt.triangle(it)[1],dt.triangle(it)[2]);
        float newRadius = CGAL::squared_distance(circum,dt.triangle(it)[2]); 
        if( newRadius > maxRadius){ 
            maxRadius = newRadius;
            Max = circum;}
    }
    


    Line line2(*prev(points.end()), *points.begin());
    float dist = CGAL::squared_distance(Max,line2);
    float minRadius = std::numeric_limits<float>::max();
    if (dist < minRadius){
            minRadius = dist;
        }
    for ( vector<altPoint>::iterator iter = points.begin(); iter != prev( points.end()); iter++){
        Line line(*iter, *next(iter));
        dist = CGAL::squared_distance(Max,line);
        if (dist < minRadius){
            minRadius = dist;
        }
    }
    
    point p = make_pair ((int)Max.x(),(int)Max.y());
    return {p,(int)(sqrt(minRadius))};      

}

// ############################################################################################
// ################################ Throwaways ################################################
// ############################################################################################
// returns true if the three points make a counter-clockwise turn

tuple<vector<point>, int, int> Prune(std::vector<point> Points, std::vector<point> Approx)
{
    vector<point> RESULT;
    bool dinmor;
    int comps = 0;
    int removed = 0;
    for (I iter = Points.begin(); iter != Points.end(); iter++)
    {
        dinmor = true;
        if (ccw(*prev((Approx.end())), *(Approx.begin()), (*iter)))
        {
            RESULT.push_back(*iter);
            dinmor = false;
        }
        comps++;
        for (I Extreme = Approx.begin(); Extreme != prev(Approx.end()); Extreme++)
        {
            if (dinmor)
            {
                comps++;
                if (ccw(*Extreme, *next(Extreme), *iter))
                {
                    RESULT.push_back(*iter);
                    dinmor = false;
                }
            }
        }
        if (dinmor)
        {
            removed++;
        }
    }
    return {PlaneSweep(RESULT), comps, removed};
}

tuple<vector<point>, int, int> CirclePrune(std::vector<point> Points, point Centrum, float Radius ){
    vector<point> RESULT;
    int comps = 0;
    int removed = 0;
    double Radius2  = pow(Radius,2);
    for(I iter = Points.begin(); iter != Points.end(); iter++){
        comps++;
        if (pow((double)(get<0>(*iter) - (double)get<0>(Centrum)),2) + 
            pow(((double)get<1>(*iter) - (double)get<1>(Centrum)),2) > Radius2  ){
            RESULT.push_back(*iter);            
        }        
        else {removed ++;}
    }
    return {PlaneSweep(RESULT), comps,removed};
}

tuple<vector<point>, int, int> ShoeLacePrune(std::vector<point> Points, std::vector<point> Approx ){
    vector<point> RESULT;
    bool dinmor;
    int comps = 0;
    int removed = 0;
    for(I iter = Points.begin(); iter != Points.end(); iter++){
        dinmor = true;
        if (ccw(*prev((Approx.end())), *(Approx.begin()),(*iter))){
            RESULT.push_back(*iter);
            dinmor =false;            
        }
        comps++;
        for( I Extreme = Approx.begin(); Extreme != prev(Approx.end()); Extreme++){        
            if (dinmor) {
                comps++;
                if ( ccw(*Extreme, *next(Extreme),*iter)){
                    RESULT.push_back(*iter);
                    dinmor = false;
                }
            }
        }
        if (dinmor) {removed ++;}
    }
    return {PlaneSweep(RESULT), comps,removed};
}

// ############################################################################
// ######################## SIMPLE THROW AWAY #################################
// ############################################################################

tuple<vector<point>, int, int> BasicThrowAway(std::vector<point> p)
{
    vector<point> max_position;
    tuple<vector<point>, int, int> RESULT;
    max_position = Extrema_8(p.begin(), p.end());
    RESULT = Prune(p, max_position);
    return {PlaneSweep(get<0>(RESULT)), get<1>(RESULT), get<2>(RESULT)};
}

// ############################################################################
// ######################## SQUARE THROW AWAY #################################
// ############################################################################

tuple<vector<point>, int, int> SquareThrowAway(std::vector<point> p)
{
    vector<point> max_position;
    tuple<vector<point>, int, int> RESULT;
    max_position = Extrema_4(p.begin(), p.end());
    RESULT = Prune(p, max_position);
    return {PlaneSweep(get<0>(RESULT)), get<1>(RESULT), get<2>(RESULT)};
}

// #############################################################################
// ######################### CIRCLE THROWAWAY ##################################
// #############################################################################


tuple<double, double, double> circleInTri(vector<point> p) {
    vector<int> xCoords, yCoords;
    for (std::size_t i = 0; i != p.size(); i++) {
        xCoords.push_back(get<0>(p[i]));
        yCoords.push_back(get<1>(p[i]));
    }
    double a = sqrt(pow((double) xCoords[1]-xCoords[2], 2.0) + pow((double) yCoords[1]-yCoords[2], 2.0));
    double b = sqrt(pow((double) xCoords[0]-xCoords[2], 2.0) + pow((double) yCoords[0]-yCoords[2], 2.0));
    double c = sqrt(pow((double) xCoords[0]-xCoords[1], 2.0) + pow((double) yCoords[0]-yCoords[1], 2.0));
    double s = (a+b+c)/2;
    cout << "abc is:" << a << " " << b << " " << c << " ";
    double r = sqrt(s*(s-a)*(s-b)*(s-c))/s;
    cout << "r is:\n" << r << "\n";
    double x_c = (a*xCoords[0]+b*xCoords[1]+c*xCoords[2])/(a+b+c);
    double y_c = (a*yCoords[0]+b*yCoords[1]+c*yCoords[2])/(a+b+c);
    tuple<double, double, double> RESULT = {x_c, y_c, r};
    cout << get<0>(RESULT) << "\n";
    cout << get<1>(RESULT) << "\n";
    cout << get<2>(RESULT) << "\n";
    return RESULT;
}

tuple<double, double, double> circleThing(std::vector<point> p) {
    vector<point> max_position = Extrema_4(p.begin(), p.end()); 
    vector<int> xCoords, yCoords;
    for (std::size_t i = 0; i < max_position.size(); i++) {
        xCoords.push_back(get<0>(max_position[i]));
        yCoords.push_back(get<1>(max_position[i]));
    }
    if (xCoords[0] == xCoords[1] && xCoords[2] == xCoords[3]) {
        double x_c = get<0>(max_position[3]) - (get<0>(max_position[0])/2);
        double y_c = get<1>(max_position[1]) - (get<1>(max_position[0])/2);
        double r = x_c - get<0>(max_position[0]);
        tuple<double, double, double> RESULT = {x_c, y_c, r};
        cout << "Perfect square with following incircle:\n";
        cout << "Center in:" << get<0>(RESULT) << get<1>(RESULT) << "Radius: " << get<2>(RESULT) << "\n";
        return RESULT;
    }
    vector<point> max_positionSort = max_position;
    std::sort(max_positionSort.begin(), max_positionSort.end());
    I ui = std::unique(max_positionSort.begin(), max_positionSort.end());
    max_positionSort.resize(std::distance(max_positionSort.begin(), ui));
    if (max_positionSort.size() == 3) {
        return circleInTri(max_positionSort);
    }
    

}

tuple<double, double, double> ApproxCircle(std::vector<point> p)
{
    vector<point> max_position = Extrema_8(p.begin(), p.end());
    vector<int> xCoords, yCoords;
    vector <double> uVals, vVals, u2Vals, v2Vals, u3Vals, v3Vals, uvVals, uvvVals, vuuVals;
    for (std::size_t i = 0; i < max_position.size(); i++) {
        xCoords.push_back(get<0>(max_position[i]));
        yCoords.push_back(get<1>(max_position[i]));
    }
    double xMean = (std::accumulate(xCoords.begin(), xCoords.end(), 0)) / max_position.size();
    double yMean = (std::accumulate(yCoords.begin(), yCoords.end(), 0)) / max_position.size();
    for (std::size_t i = 0; i < max_position.size(); i++) {
        uVals.push_back((double)(xCoords[i]-xMean));
        vVals.push_back((double)(yCoords[i]-yMean));
        u2Vals.push_back(pow(uVals[i], 2.0));
        v2Vals.push_back(pow(vVals[i], 2.0));
        u3Vals.push_back(pow(uVals[i], 3.0));
        v3Vals.push_back(pow(vVals[i], 3.0));
        uvVals.push_back(uVals[i]*vVals[i]);
        uvvVals.push_back(uVals[i]*vVals[i]*vVals[i]);
        vuuVals.push_back(vVals[i]*uVals[i]*uVals[i]);
    }
    double S_uu = (std::accumulate(u2Vals.begin(), u2Vals.end(), 0));
    double S_vv = (std::accumulate(v2Vals.begin(), v2Vals.end(), 0));
    double S_uv = (std::accumulate(uvVals.begin(), uvVals.end(), 0));
    double S_uuu = (std::accumulate(u3Vals.begin(), u3Vals.end(), 0));
    double S_vvv = (std::accumulate(v3Vals.begin(), v3Vals.end(), 0));
    double S_uvv = (std::accumulate(uvvVals.begin(), uvvVals.end(), 0));
    float S_vuu = (std::accumulate(vuuVals.begin(), vuuVals.end(), 0));
    double det = S_uu*S_vv - S_uv*S_uv;
    double a = 1/2*(S_uuu + S_uvv);
    double b = 1/2*(S_vvv + S_vuu);
    double u_c, v_c;
    //Solve equations u_c*S_uu + v_c*S_uv = a 
    // And u_c*S_uv + v_c*S_vv = b
    if (det != 0) {
        double u_c = (a*S_vv - S_uv*b)/det;
        double v_c = (S_uu*b - a*S_uv)/det;
    }    
    double x_c = u_c + xMean;
    double y_c = v_c + yMean;
    double alpha = pow(u_c, 2.0)+pow(v_c, 2.0)+(S_uu+S_vv)/max_position.size();
    double radius = sqrt(alpha);
    tuple<double, double, double> RESULT = {x_c, y_c, radius};
    cout << radius; 
    return RESULT;
}



tuple<vector<point>, int, int> CircleThrowAway(std::vector<point> p)
{
    vector<point> max_position = Extrema_4(p.begin(), p.end());
    tuple<vector<point>, int, int> RESULT;
    //RESULT = CirclePrune(p,get<0>(Circle),get<1>(Circle));
    return {PlaneSweep(get<0>(RESULT)), get<1>(RESULT),get<2>(RESULT)};
    
}

tuple<vector<point>, int, int> LauneyThrowAway(std::vector<point> p) {
    vector<point> max_position = Extrema_8(p.begin(), p.end());
    tuple<point, int> Circle = InCircle(max_position);
    tuple<vector<point>, int, int> RESULT;
    RESULT = CirclePrune(p,get<0>(Circle),get<1>(Circle));
    return {PlaneSweep(get<0>(RESULT)), get<1>(RESULT),get<2>(RESULT)};
}



// #############################################################################
// ######################### Shoelace Throwaway ################################
// #############################################################################
// Function to find eight extreme points
point calcCenter(I first, I last)
{
    float x = 0;
    float y = 0;
    float n = 0;
    for (I iter = first; iter != last; iter++)
    {
        x += get<0>(*iter);
        y += get<1>(*iter);
        n++;
    }
    x = x / n;
    y = y / n;
    return make_pair((int)x, (int)y);
}
tuple<vector<point>, int, int> ShoelaceThrowAwayHull(std::vector<point> p)
{
    vector<point> max_position;
    vector<point> RESULT;
    max_position = Extrema_4(p.begin(), p.end());
    bool dinmor;
    int comps = 0;
    int removed = 0;
    int Shoelace = 0;
    point center = calcCenter(max_position.begin(),max_position.end());
    return {PlaneSweep(RESULT), comps,removed};

}

//#############################################################################
//######################### STRUCTURES ########################################
//############################################################################

// ###########################################################################
// ####################### SIMPLE POINT PLANE ################################
// ###########################################################################

// Implement a struct, which contains the convex hull and number of comparisons.
struct PointPlane
{
    int NrComps;
    int removed;
    vector<point>::iterator ConveXHullstart;
    vector<point>::iterator ConveXHullEnd;
    vector<point> ConveXHull;
    vector<point> AllPoints;

    int GenerateSquarePoints (int RangeX,int RangeY, int n){
        unsigned seed = 0;
        std::default_random_engine gen (seed);
        std::uniform_int_distribution<int> uniform_dist_X(0, RangeX);
        std::uniform_int_distribution<int> uniform_dist_Y(0, RangeY);
        vector<point> points;
        for ( int i = 1; i <= n  ; i++){
            int x = uniform_dist_X(gen);
            int y = uniform_dist_Y(gen);
            point p = make_pair(x,y);   
            this->AllPoints.push_back(p);
        };
        return 0;
    };
    int GenerateCirclePoints (int Radius, int n){
        unsigned seed = 0;
        std::default_random_engine gen (seed);
        std::uniform_int_distribution<int> RanInt(-Radius, Radius);
        vector<point> points;
        for ( int i = 1; i <= n  ; i++){
            int x = RanInt(gen);
            int y = RanInt(gen);
            while ( pow(x,2) + pow(y,2) > Radius ){
                x = RanInt(gen);
                y = RanInt(gen);
            }
            point p = make_pair(x,y);   
            this->AllPoints.push_back(p);
        };
        return 0;
    };
    /*int GenerateCirclePoints2 (int Radius, int n){
        unsigned seed = 0;
        std::default_random_engine gen (seed);
        std::uniform_real_distribution<> uniform_dist_Radius(-(float)Radius, (float)Radius);
        std::uniform_real_distribution<> Angle(0.0, 1.0);
        vector<point> points;
        for ( int i = 1; i <= n  ; i++){
            float r = uniform_dist_Radius(gen);
            float theta = Angle(gen) * 2 * M_PI;
            int x = r * cos(theta);
            int y = r * sin(theta);
            point p = make_pair(x,y);   
            this->AllPoints.push_back(p);
        };
        return 0;
    };*/




    int GetPoints (){
        vector<point>::iterator iter;
        std::ofstream outfile;
        outfile.open("points.csv", std::ios_base::app);
        outfile << "x,y\n";
        for (iter = this->AllPoints.begin(); iter != this->AllPoints.end(); iter++)
        {

            outfile << get<0>(*iter) << " ,";
            outfile << get<1>(*iter) << " \n";
        }
        outfile.close();
        return 0;
    };
    int GetHull()
    {
        vector<point>::iterator iter;
        ofstream outfile;
        outfile.open("Hull.csv", std::ios_base::app);
        outfile << "x,y\n";
        for (vector<point>::iterator iter = this->ConveXHull.begin(); iter != this->ConveXHull.end(); iter++)
        {
            outfile << get<0>(*iter) << " ,";
            outfile << get<1>(*iter) << " \n";
        }
        outfile.close();
        return 0;
    };
    int GetCompsAndRemoved()
    {
        std::cout << "THIS is Nr of Comps " << this->NrComps << "\n";
        std::cout << "THIS is Nr of removed elements " << this->removed << "\n";
        return 0;
    }

    int solve()
    {
        this->ConveXHull = PlaneSweep(this->AllPoints);
        return 0;
    }
};

// #########################################################################
// ########################### BASIC THROWAWAY##############################
// #########################################################################
struct BaseHull : PointPlane
{
    int ThrowAwayHull()
    {
        tuple<vector<point>, int, int> res = BasicThrowAway(this->AllPoints);
        this->ConveXHull = get<0>(res);
        this->NrComps = get<1>(res);
        this->removed = get<2>(res);
        return 0;
    }
};

// ###############################################################
// ########################## SQUARE ThrowAWAY ####################
// ###############################################################
struct SquareHull : PointPlane
{
    int ThrowAwayHull()
    {
        tuple<vector<point>, int, int> res = SquareThrowAway(this->AllPoints);
        this->ConveXHull = get<0>(res);
        this->NrComps = get<1>(res);
        this->removed = get<2>(res);
        return 0;
    }
};


// ###############################################################
// ########################## CircleThrowAway ThrowAWAY ####################
// ###############################################################
struct CircleHull : PointPlane {
    int ThrowAwayHull(){
        tuple<vector<point>,int,int> res = LauneyThrowAway(this->AllPoints);
        this->ConveXHull = get<0>(res);
        this->NrComps = get<1>(res);
        this->removed =get<2>(res);
        return 0;
    }
};


int main() {   
    CircleHull plane;
    BaseHull plane2;
    BaseHull plane3;
    // plane.GenerateSquarePoints(100,100,100);
    plane.GenerateCirclePoints(1000,1000);
    plane.GetPoints();
    plane.ThrowAwayHull();
    plane.GetHull();
    plane.GetCompsAndRemoved();
    InCircle(plane.AllPoints);
    return 0;
}
