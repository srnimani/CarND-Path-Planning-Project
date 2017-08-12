#ifndef CAR_H
#define CAR_H
#include <iostream>
#include <math.h>
#include <vector>
#include <string>

using namespace std;

/* format for each car looks like this, [ id, x, y, vx, vy, s, d]. The id is a unique identifier for that car although that same car will always end up in the same position in the vector as well. The x, y values are in global map coordinates, and the vx, vy values are the velocity components also in reference to the global map. Finally s and d are the Frenet coordinates for that car.
*/

class Car
{
public:
    
    int id;
    double x ;
    double y ;
    double vx ;
    double vy ;
    double s ;
    double d ;
    
    Car(int id, double x, double y, double vx, double vy, double s, double d):
    id(id), x(x), y(y), vx(vx), vy(vy), s(s), d(d) { };

    /**
     * Constructor
     */
    Car() {};
    
    /**
     * Destructor
     */
    ~Car() {};
    /*
    const int get_car_lane()
    {
        int lane = 0;
        if (d >=0 && d < 4)  lane = 0;
        if (d >=4 && d < 8)  lane = 1;
        if (d >=8 && d < 12) lane = 2;
        return lane;
    }
    */
    
    const int get_car_lane(double d)
    {
        return int(d/4.0) ;
    }

    
    const double get_car_velocity()
    {
        return sqrt(vx * vx + vy * vy) ;
    }
    
};

#endif // CAR-H
