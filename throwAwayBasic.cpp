#include <iostream>
#include <vector> 
#include <stdlib.h>
#include <cassert>
#include <iterator>
#include <algorithm>    

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
};
constexpr void swap( vector<point>::iterator p1,
 vector<point>::iterator p2){
    iter_swap(p1,p2);
};

//Equality operator for our points
bool operator==(const point& p1, const point& p2) {
    return (p1.x == p2.x && p1.y == p2.y);
}
// NameSpace used for validation:
namespace validation {
    bool same_multiset(I p, I q, I r, I s) {
        using P = typename std::iterator_traits<I>::value_type;
        using S = std::vector<P>;
        std::size_t n = std::distance(p, q);
        std::size_t m = std::distance(r, s);
        if (m != n) {
            return false;
        }
        S backup; 
        backup.resize(n);
        std::copy(p, q, backup.begin());
        std::sort(backup.begin(), backup.end(),
        [](P const& a, P const& b) -> bool {
            return (a.x< b.x) or (a.x == b.x and a.y < b.y);
        });
        S other; 
        other.resize(n);
        std::copy(r, s, other.begin());
        std::sort(other.begin(), other.end(),
        [](P const& a, P const& b) -> bool {
            return (a.x < b.x or (a.x == b.x and a.y < b.y);
        });
        return std::equal(backup.begin(), backup.end(), other.begin(),
        [](P const& a, P const& b) -> bool {
            return (a.x == b.x) and (a.y == b.y);
        });
    }
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
    };
    int PlaneSweep () {
        return 0;
    }
};

// #################################################################
// ############# SIMPLE PLANESWEEP HULL ############################
//##################################################################



namespace planesweep {

bool comp(point a, point b){
    return (a. x < b. x) or (a. x == b. x and a. y < b. y);
}


template<typename I>
void swap_blocks(I source, I past, I target) {
    // retains the order in [ source , past)
    if (source == target || source == past) {
        return;
    }
    using P = typename std :: iterator_traits<I> :: value_type;
    I hole = target;
    P p = * target;
    I const last = std :: prev(past) ;
    while (true) {
        * hole = * source;
        ++hole;
        if (source == last) {
            break;
        }
        *source = *hole;
        ++source
    }
    *source = p;
}




template<typename I>
std :: pair<I, I> find_poles(I first, I past) {
using P = typename std : : iterator_traits<I> : : value_type;
auto pair = std :: minmax_element(first, past,comp)
return pair;
}

// remove all element on top of each oth
template<typename I, typename C>
    std :: pair<I, I> clean(I pole, I rest, C compare) {
    assert(pole != rest) ;
    I k = std :: next(pole) ;
    while (k !=rest && ( * k).x == ( * pole) . x) {
        ++k;
    }
    if (k == std :: next(pole)) {
        return std :: make_pair(pole, std :: next(pole)) ; 
    }
    I bottom = std :: min_element(pole, k, compare) ;
    std :: iter_swap(pole, bottom) ;
    I top = std :: max_element(std :: next(pole) , k, compare) ;
    if (( * top) . y == ( * pole) . y) {
    return std :: make_pair(pole, k) ;
    }
    std :: iter_swap(std :: next(pole) , top) ;
    return std :: make_pair(std :: next(pole) , k) ;
}
// using left turn test
template<typename I>
    I scan(I first, I top, I next, I past) {
    assert(first != past) ;
    for (I i = next; i != past; ++i) {
        while (top != first && not right_turn( * std :: prev(top) , * i) ) {
            -- top;
        }
        ++top;
        std :: iter_swap(i, top) ;
    }
    return ++top;
}
template<typename I>
I scan_one_more(I first, I top, I extra) {
    while (top  != first && not right_turn( * std :: prev(top) , * extra) ) {
        --top;
    }
    return ++top;
}

// Brute Force part

template<typename I, typename C>
I brute_force(I first, I past, C compare) {
    std :: size_t n = std : : distance(first, past) ;
    assert(n <= 4) ;
    if (n <= 1) {
        return past;
    }
    else if (n == 2) {
        if ( * first == * std : : next(first) ) {
            return std : : next(first) ;
        }
        if (not compare( * first, * std :: next(first) ) ) {
            std :: iter_swap(first, std :: next(first) ) ;
        }
        return past;
    }
    if (n == 3) {
        I west = std :: min_element(first, past, compare) ;
        std :: iter_swap(first, west) ;
        I east = std :: max_element(std :: next(first) , past, compare) ;
        if ( * first == * east) {
            return std : : next(first) ;
        }
        std :: iter_swap(std :: next(std :: next(first) ) , east) ;
        if (right_turn( * first, * std :: next(first) , * std : : next(std : : next(first) ) ) ) {
            return std : : next(std :: next(std :: next(first) ) ) ;
        }
        std :: iter_swap(std :: next(first) , std :: next(std :: next(first) ) ) ;
        return std :: next(std :: next(first) ) ;
    }
// n == 4
    I west = std :: min_element(first, past, compare) ;
    std :: iter_swap(first, west) ;
    I east = std :: max_element(std :: next(first) , past, compare) ;
    if ( * first == * east) {
        return std :: next(first) ;
    }
    std :: iter_swap(std :: next(std :: next(std :: next(first) ) ) , east) ;

    if (not compare( * std :: next(first) , * std :: next(std :: next(first) ) ) ) {
        std :: iter_swap(std :: next(first) , std :: next(std :: next(first) ) ) ;
    }
    I rest = scan(first, first, std :: next(first) , past) ;
    return rest;
}
template<typename I, typename C>
I conquer(I first, I past, C compare) {
    std :: size_t n = std :: distance(first, past) ;
    if (n <= 4) {
        return brute_force(first, past, compare) ;
    }
    I middle = first;
    std :: advance(middle, n / 2) ;
    I rest1 = conquer(first, middle, compare) ;
    I rest2 = conquer(middle, past, compare) ;
    std :: size_t h2 = std :: distance(middle, rest2 ) ;
    I rest = rest1 ;
    std :: advance(rest, h2 ) ;
    swap_blocks(middle, rest2 , rest1 ) ;
    std :: inplace_merge(first, rest1 , rest, compare) ;
    std :: pair<I, I> pair = clean(first, rest, compare) ;
    I top = std :: get<0> (pair) ;
    I next = std :: get<1> (pair) ;
    rest = scan(first, top, next, rest) ;
    return rest;
}





int orientation(point p, point q, point r) 
{ 
    int val = (q.y - p.y) * (r.x - q.x) - 
              (q.x - p.x) * (r.y - q.y); 
  
    if (val == 0) return 0;  // colinear 
    return (val > 0)? 1: 2; // clock or counterclock wise 
} 



vector<point>::iterator solve(vector<point>:: iterator first,
           vector<point>:: iterator past, int n) { 
    // There must be at least 3 points 
    if (n < 3) return; 
  
    // Initialize Result 
    vector<point> hull; 
  
    // Here we find the dividing line between the upper and lower part of the hull
    pair<vector<point>::iterator,
         vector<point>::iterator> POLES = find_poles(first,past);
    vector<point>::iterator left =  first;
    vector<point>::iterator right =  std::prev(past);
    iter_swap(std :: get<0> (POLES), left);
    iter_swap(std :: get<0> (POLES), right);
    int upper_limit = std :: max(( * left) . y, ( * right) . y) ;
    int lower_limit = std :: min(( * left) . y, ( * right) . y) ;
    // Now we partition the Hull in two part with partition:
    I i = std :: partition(std :: next(left) , right,
        [&](point const& q) {
            if (q. y >= upper_limit) {
                return true;
            }
            else if (q.y <= lower_limit) {
                return false ;
            }
            return not isLeftTurn( *left, q, *right) ;
        }) ;

    std :: iter_swap(i, right);

    // Now The Entire partition is complete and the lower partition begins with right.
    // Now we sort it and planesweep
    
    // Sort and sweep lower hull
    using P = typename std :: iterator_traits<I> :: value_type;
    vector<point>::iterator middle = std :: next( right);
    vector<point>::iterator rest1 = conquer(first, middle,
        [ ](P const& a, P const& b) {
            return a. x < b. x or (a. x == b. x and a. y < b. y) ;
        }) ;
    if (middle == past) {
        return rest1 ;
    }
    vector<point>::iterator west = first;
    std :: iter_swap(std :: prev(rest1 ) , std :: prev(middle) ) ;
    // Sort and sweep lower hull
    vector<point>::iterator east = std :: prev(middle) ;
    vector<point>::iterator rest2 = conquer(east, past,
    [ ](P const& a, P const& b) {
        return a. x > b. x or (a. x == b. x and a. y > b. y) ;
    });
    if (rest2 != east) {
        rest2 = scan_one_more(east, std :: prev(rest2 ) , west) ;
    }
    std :: size_t h2 = std :: distance(east, rest2 ) ;
    vector<point>::iterator rest = rest1 ;
    std :: advance(rest, h2 - 1) ;
    swap_blocks(east, rest2 , std :: prev(rest1 ) ) ;
    return rest;
      
}
// #################################################################
// ############# SIMPLE THROWAWAY HULL #############################
//##################################################################
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
        I mid = std::partition(first, past, [&](P const& r) -> bool {
            return not (r == *polygon);
        });
        return mid;
    }
    if (std::distance(polygon, rest) == 2) {
        I mid = std::partition(first, past, [&](P const& r) -> bool {
            return not (onLineSegment(*polygon, *std::next(polygon), r));
        });
        return mid;
    }
    assert(std::distance(polygon, rest) > 2);
    std::vector<P> apprHull;
    for (I k = polygon; k != rest; ++k) {
        apprHull.push_back(*k);
    }
    apprHull.push_back(*polygon);
    I upper = apprHull.begin();
    if ((*upper).x == (*std::next(upper)).x) {
        ++upper;
    }
    I lower = upper; 
    I last = apprHull.end();
    for (I k = std::next(upper); k != last; ++k) {
        if ((*lower).y < (*k).y) {
            lower = k;
        }
    }
    I uend = std::next(lower);
    if ((*lower).x== (*std::next(lower)).x) {
        ++lower;
    }
    I lend = apprHull.end();
    I mid = std::partition(first, past, [&](P const& r) -> bool {
        return not betweenPolynomialChains(r, upper, uend, lower, lend);
    });
    return mid; 
}

I prune(I first, I past) {
    if (std::distance(first, past) < 2) {
        return past; 
    }
    //Step 1
    std::vector<point> extreme = findExtrema(first, past);
    std::sort(extreme.begin(), extreme.end());
    removeDuplicates(extreme);
    I q = first;
    for (int i = 0; i != extreme.size(); ++i, ++q) {
        std::iter_swap(q, extreme[i]);
    }
    //Step 2
    I rest = plane_sweep::solve(first, q); 
    //Step 3
    I interior = eliminateInnerPoints(rest, past, first, rest);
    return interior;
}

I solve(I first, I past) {
    I rest = prune(first, past);
#if defined(PRUNING)
    prunningEfficiency += std::distance(first, rest);

#endif
    return plane_sweep::solve(first, rest);
}
// #############################################################
// ############## END Of Simple ThrowAway##################
// #############################################################
bool check(I first, I past) {
    using P = I; 
    using S = std::vector<I>;
    using J = typename S::iterator;
    S data;

    std::size_t n = std::distance(first, past);
    data.resize(n);
    std::copy(first, past, data.begin());
    J rest = solve(data.begin(), data.end());
    bool ok = (
        validation::same_multiset(data.begin(), data.end(), first, past) and
        validation::convex_polygon( data.begin(), rest) and 
        validation::all_inside(rest, data.end(), data.begin(), rest)
    );
    return ok;
}

// The general idea is, we create classes which inherit the baseclass 
// and have a additional function, which implements different throwaway tactics.
struct BasicThrowAway : PointPlane {
    int ThrowAwayHull(){
        std::cout  << " Dette skal implementeres ";
        return 0;
    }
};

int main() {   
    
    PointPlane plane;
    plane.GeneratePointList(20,100,10);
    plane.GetHullPoints();
    return 0;
};
