#include <fstream>
#include <math.h>
#include <uWS/uWS.h>
#include <chrono>
#include <iostream>
#include <thread>
#include <vector>
#include "Eigen-3.3/Eigen/Core"
#include "Eigen-3.3/Eigen/QR"
#include "json.hpp"

#include "spline.h"
#include "car.h"

using namespace std;

// for convenience
using json = nlohmann::json;

int current_lane            = 1;    // We are starting in lane 1

int next_lane               = 1;    // We are starting in lane 1

double ref_velocity         = 0;    // velocity we want to achieve, will change

const double max_velocity   = 49.0;  // Maximum velocty w/o violating speed limit with some margin

const double safe_distance  = 30 ;  // in meters

// Pathplanner states and initialized states
enum {
    _STARTING,      // Starting state
    _KL,            // Keep lane
    _PLCL,          // Plan Lane change left
    _CLL,           // Change lane left
    _PLCR,          // Plan Lane change right
    _CLR,           // Change lane right
};

unsigned int PP_current_state = _STARTING;
unsigned int PP_next_state    = _STARTING;

vector<Car> cars; // All cars in the vicinity

double track_length = 6945.554;



// For converting back and forth between radians and degrees.
constexpr double pi() { return M_PI; }
double deg2rad(double x) { return x * pi() / 180; }
double rad2deg(double x) { return x * 180 / pi(); }


// Checks if the SocketIO event has JSON data.
// If there is data the JSON object in string format will be returned,
// else the empty string "" will be returned.
string hasData(string s) {
  auto found_null = s.find("null");
  auto b1 = s.find_first_of("[");
  auto b2 = s.find_first_of("}");
  if (found_null != string::npos) {
    return "";
  } else if (b1 != string::npos && b2 != string::npos) {
    return s.substr(b1, b2 - b1 + 2);
  }
  return "";
}

double distance(double x1, double y1, double x2, double y2)
{
	return sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));
}

int ClosestWaypoint(double x, double y, vector<double> maps_x, vector<double> maps_y)
{

	double closestLen = 100000; //large number
	int closestWaypoint = 0;

	for(int i = 0; i < maps_x.size(); i++)
	{
		double map_x = maps_x[i];
		double map_y = maps_y[i];
		double dist = distance(x,y,map_x,map_y);
		if(dist < closestLen)
		{
			closestLen = dist;
			closestWaypoint = i;
		}

	}

	return closestWaypoint;

}

int NextWaypoint(double x, double y, double theta, vector<double> maps_x, vector<double> maps_y)
{

	int closestWaypoint = ClosestWaypoint(x,y,maps_x,maps_y);

	double map_x = maps_x[closestWaypoint];
	double map_y = maps_y[closestWaypoint];

	double heading = atan2( (map_y-y),(map_x-x) );

	double angle = abs(theta-heading);

	if(angle > pi()/4)
	{
		closestWaypoint++;
	}

	return closestWaypoint;

}

// Transform from Cartesian x,y coordinates to Frenet s,d coordinates
vector<double> getFrenet(double x, double y, double theta, vector<double> maps_x, vector<double> maps_y)
{
	int next_wp = NextWaypoint(x,y, theta, maps_x,maps_y);

	int prev_wp;
	prev_wp = next_wp-1;
	if(next_wp == 0)
	{
		prev_wp  = maps_x.size()-1;
	}

	double n_x = maps_x[next_wp]-maps_x[prev_wp];
	double n_y = maps_y[next_wp]-maps_y[prev_wp];
	double x_x = x - maps_x[prev_wp];
	double x_y = y - maps_y[prev_wp];

	// find the projection of x onto n
	double proj_norm = (x_x*n_x+x_y*n_y)/(n_x*n_x+n_y*n_y);
	double proj_x = proj_norm*n_x;
	double proj_y = proj_norm*n_y;

	double frenet_d = distance(x_x,x_y,proj_x,proj_y);

	//see if d value is positive or negative by comparing it to a center point

	double center_x = 1000-maps_x[prev_wp];
	double center_y = 2000-maps_y[prev_wp];
	double centerToPos = distance(center_x,center_y,x_x,x_y);
	double centerToRef = distance(center_x,center_y,proj_x,proj_y);

	if(centerToPos <= centerToRef)
	{
		frenet_d *= -1;
	}

	// calculate s value
	double frenet_s = 0;
	for(int i = 0; i < prev_wp; i++)
	{
		frenet_s += distance(maps_x[i],maps_y[i],maps_x[i+1],maps_y[i+1]);
	}

	frenet_s += distance(0,0,proj_x,proj_y);

	return {frenet_s,frenet_d};

}

// Transform from Frenet s,d coordinates to Cartesian x,y
vector<double> getXY(double s, double d, vector<double> maps_s, vector<double> maps_x, vector<double> maps_y)
{
	int prev_wp = -1;

	while(s > maps_s[prev_wp+1] && (prev_wp < (int)(maps_s.size()-1) ))
	{
		prev_wp++;
	}

	int wp2 = (prev_wp+1)%maps_x.size();

	double heading = atan2((maps_y[wp2]-maps_y[prev_wp]),(maps_x[wp2]-maps_x[prev_wp]));
	// the x,y,s along the segment
	double seg_s = (s-maps_s[prev_wp]);

	double seg_x = maps_x[prev_wp]+seg_s*cos(heading);
	double seg_y = maps_y[prev_wp]+seg_s*sin(heading);

	double perp_heading = heading-pi()/2;

	double x = seg_x + d*cos(perp_heading);
	double y = seg_y + d*sin(perp_heading);

	return {x,y};

}


int get_lane_number(double d)
{
    return int(d/4.0);
}


double check_space(double car_s, double lane_to_check)
{
    double freespace = 1000000 ; // assume free space to start with
    double current_gap = 0;
    
    double s1, s2;  // segment from s1 to s2 to check

    s1 = car_s - safe_distance;
    s2 = car_s + safe_distance;
    
    if ((s1 < 0) || (s2 > track_length)) // currently not allowing any passes around the zero crossing
    {
        return 0 ; // declare no free space or gap
    }
    
    for(int i = 0; i < cars.size(); ++i)
    {
        if (lane_to_check == get_lane_number(cars[i].d))
        {
            current_gap = abs(cars[i].s - car_s) ;
            if (current_gap < safe_distance) // lane is unsafe
            {
                freespace = current_gap ;
                cout << "Car_id:\t" << i << "\tlane:\t" << lane_to_check << "\tFree Sapce:\t " << freespace << endl;
                return freespace; // Exit as we have found a close car, no need to continue checking
            }
            else
            {
                if (freespace > current_gap) freespace = current_gap;
            }
        }
    }
    return freespace;
}



// Check available lane to switch to, return -1 if no lane available
int available_lane(double car_s, const int our_lane)
{
    int free_lane           = -1 ; // Means no lane free
    double left_free_space  =  0 ; // Assume no space to start with
    double right_free_space =  0 ; // Assume no space to start with
    
    cout << "Safe Distance: \t" << safe_distance << endl;

    switch (our_lane)
    
    {
        case 0 : // We are in the left most lane
        {
            right_free_space = check_space(car_s, our_lane + 1); // returns the distance to the closest car on right
            if (right_free_space > safe_distance) free_lane = our_lane + 1;
            break;
        }
            
        case 1 : // we are in the middle lane lane
        {
            right_free_space = check_space(car_s, our_lane + 1); // returns the distance to the closest car on right
            left_free_space = check_space(car_s, our_lane - 1) ; // returns the distance to the closest car on left
            if ((right_free_space > safe_distance * 1.5) && (right_free_space > left_free_space)) free_lane = our_lane + 1;
            if ((left_free_space > safe_distance * 1.5) && (left_free_space >= right_free_space)) free_lane = our_lane - 1;
            break;
        }

        case 2 : // We are in the right most lane
        {
            left_free_space = check_space(car_s, our_lane - 1); // returns the distance to the closest car on left
            if (left_free_space > safe_distance) free_lane = our_lane - 1;
            break;
        }
            
        default:
        {
            break;
        }
    }
    return free_lane;
}


// Car ahead of us, within 'collision distance'..

bool car_ahead_is_close(double our_car_s, int our_lane, int prev_size, double &velocity_of_car_ahead)
{
    int number_of_cars  = cars.size();
    
    bool too_close      = false;
    
    for(int i = 0; i < number_of_cars; ++i)
    {
        if(our_lane == get_lane_number(cars[i].d))
        {
            double check_speed  = cars[i].get_car_velocity(); // get the speed of the car ahead
            double check_car_s  = cars[i].s;
            
            check_car_s += ((double)prev_size * 0.02 * check_speed);
            
            // check s values greater than our car and s gap
            
            if ((check_car_s > our_car_s) && (check_car_s - our_car_s) < safe_distance)
            {
                too_close = true;
                velocity_of_car_ahead = check_speed;
                return too_close;
            }
        }
    }
    return too_close ;
}

 

//////////////////////////////////////////////////////////////////////////////////////////////////////////////


int main()
{
  uWS::Hub h;

  // Load up map values for waypoint's x,y,s and d normalized normal vectors
  vector<double> map_waypoints_x;
  vector<double> map_waypoints_y;
  vector<double> map_waypoints_s;
  vector<double> map_waypoints_dx;
  vector<double> map_waypoints_dy;

  // Waypoint map to read from
  string map_file_ = "../data/highway_map.csv";
  // The max s value before wrapping around the track back to 0
  double max_s = 6945.554;

  ifstream in_map_(map_file_.c_str(), ifstream::in);

  string line;
  while (getline(in_map_, line))
  {
  	istringstream iss(line);
  	double x;
  	double y;
  	float s;
  	float d_x;
  	float d_y;
  	iss >> x;
  	iss >> y;
  	iss >> s;
  	iss >> d_x;
  	iss >> d_y;
  	map_waypoints_x.push_back(x);
  	map_waypoints_y.push_back(y);
  	map_waypoints_s.push_back(s);
  	map_waypoints_dx.push_back(d_x);
  	map_waypoints_dy.push_back(d_y);
  }
    
    

  h.onMessage([&map_waypoints_x,&map_waypoints_y,&map_waypoints_s,&map_waypoints_dx,&map_waypoints_dy](uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length, uWS::OpCode opCode)
    {
    // "42" at the start of the message means there's a websocket message event.
    // The 4 signifies a websocket message
    // The 2 signifies a websocket event
    //auto sdata = string(data).substr(0, length);
    //cout << sdata << endl;
    if (length && length > 2 && data[0] == '4' && data[1] == '2') {

      auto s = hasData(data);

      if (s != "") {
        auto j = json::parse(s);
        
        string event = j[0].get<string>();
        
        if (event == "telemetry") {
          // j[1] is the data JSON object
          
        	// Main car's localization Data
          	double car_x = j[1]["x"];
          	double car_y = j[1]["y"];
          	double car_s = j[1]["s"];
          	double car_d = j[1]["d"];
          	double car_yaw = j[1]["yaw"];
          	double car_speed = j[1]["speed"];

          	// Previous path data given to the Planner
          	auto previous_path_x = j[1]["previous_path_x"];
          	auto previous_path_y = j[1]["previous_path_y"];
          	// Previous path's end s and d values 
          	double end_path_s = j[1]["end_path_s"];
          	double end_path_d = j[1]["end_path_d"];

          	// Sensor Fusion Data, a list of all other cars on the same side of the road.
          	vector<vector<double>> sensor_fusion = j[1]["sensor_fusion"];
            
            // Some local variables defined
            
            int prev_size       = previous_path_x.size();
            
            int current_lane    = get_lane_number(car_d);
            
            int last_lane       = get_lane_number(end_path_d); // to avoid double switching when we just changed lanes
            
            int free_lane       = current_lane ;
            
            if (prev_size > 0)
            {
                car_s = end_path_s;
            }
            
            double velocity_of_car_ahead = 0; // Used while maintaining the same lane behind a slow car
            
            // Get info. of all cars in the vicinity in to a class structure
            
            cars.clear();
            
            for(int i = 0; i < sensor_fusion.size(); ++i)
            {
                Car car(sensor_fusion[i][0], sensor_fusion[i][1], sensor_fusion[i][2], sensor_fusion[i][3], sensor_fusion[i][4], sensor_fusion[i][5], sensor_fusion[i][6]);
                cars.push_back(car);
            }
            
            //
            ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
            //
            // Path planner state machine .. transitions from start to other states based on conditions
            //
            ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
            //
            
            switch (PP_current_state)
            {
                case _STARTING: // Car just starting
                {
                    if (car_ahead_is_close(car_s, current_lane, prev_size, velocity_of_car_ahead))
                    {
                        ref_velocity -= 0.3;
                        PP_next_state = _STARTING ; // Maintain same state
                    }
                    else if (ref_velocity < max_velocity) // No car close, so ramp up speed
                    {
                        ref_velocity += 0.4 ;
                        PP_next_state = _STARTING ; // Maintain same lane
                    }
                    else
                    {
                        ref_velocity  = max_velocity ; // maximum speed
                        PP_next_state = _KL ; // switch to next state
                    }
                    next_lane = current_lane ;
                    cout << "Starting up.. new speed:\t " << ref_velocity << endl;
                    break;
                }
                    
                case _KL:
                {
                    if (!car_ahead_is_close(car_s, current_lane, prev_size, velocity_of_car_ahead))
                    // No need to change lane and keep maximum speed
                    {
                        if (ref_velocity < max_velocity)
                        {
                            ref_velocity += 0.3 ;
                        }
                        cout << "Staying in the same lane.. no need to switch.. current Lane :\t " << current_lane << endl;
                        PP_next_state = _KL ; // Maintain same lane
                    }
                    else
                    {
                        cout << "Car ahead is going slow and we are closing in fast.. check if we can swich lanes.. " << endl;
                        free_lane = available_lane(car_s, current_lane);
                        cout << "Free lane :\t" << free_lane << endl;
                        
                        switch (free_lane)
                        {
                            case 0: // Left lane is free
                            {
                                PP_next_state = _PLCL ;
                                cout << "Switching to left lane.. current_lane: \t " << current_lane << endl;
                                break;
                            }
                                
                            case 1:
                            {
                                if (last_lane == 0) // we are moving from left most lane
                                    PP_next_state = _PLCR ;
                                else PP_next_state = _PLCL ; // we are moving from right lane
                                break;
                            }
                                
                            case 2:
                            {
                                PP_next_state = _PLCR ;
                                cout << "Switching to right lane.. current_lane: \t " << current_lane << endl;
                                break;
                            }
                                
                            case -1: // No free lanes, maintain lane at reduced speed
                            {
                                if (ref_velocity > velocity_of_car_ahead)
                                {
                                    ref_velocity -= 0.3 ; // reduce speed
                                    if (ref_velocity < velocity_of_car_ahead)
                                        ref_velocity = velocity_of_car_ahead; // no need to reduce below the speed of the car ahead
                                    cout << "Can't switch lanes now.. staying in the same lane at the same speed of the car ahead\t" << "Lane :\t " << current_lane << endl;
                                }
                                PP_next_state = _KL ;
                                break;
                            }
                            default: // Maintain lane at reduced speed
                            {
                                if (ref_velocity > velocity_of_car_ahead)
                                {
                                    ref_velocity -= 0.3 ; // reduce speed
                                    if (ref_velocity < velocity_of_car_ahead)
                                        ref_velocity = velocity_of_car_ahead; // no need to reduce below the speed of the car ahead
                                    cout << "Default state .. how did we reach here???......" << endl;
                                }
                                PP_next_state = _KL ;
                                break;
                            }
                        }
                    }
                    break;
                }
                    
                case _PLCL:
                {
                    // Get the path adjusted to left lane
                    next_lane = current_lane - 1; // Move to left lane
                    PP_next_state = _CLL ;
                    cout << "Switching to left Lane :\t " << next_lane << endl;
                    
                    break;
                }
                    
                case _CLL:
                {
                    cout << "Switched to left lane :\t " << next_lane << endl;
                    
                    // Wait till the car fully moved to the lane before changing state
                    if ((next_lane == last_lane) && (current_lane == next_lane))
                    {
                        PP_next_state = _KL ; // swicthed to new lane fully, go to KL state
                    }
                    else
                    {
                        PP_next_state = _CLL ;
                    }
    
                    break;
                }
                    
                case _PLCR:
                {
                    // Get the path adjusted to right lane
                    next_lane = current_lane + 1; // Move to right lane
                    cout << "Switching to right Lane :\t " << next_lane << endl;
                    PP_next_state = _CLR ;
                    
                    break;
                }
                    
                case _CLR:
                {
                    cout << "Switched to right lane :\t " << next_lane << endl;
                    
                    // Wait till the car fully moved to the lane before changing state
                    if ((next_lane == last_lane) && (current_lane == next_lane))
                    {
                        PP_next_state = _KL ; // swicthed to new lane fully, go to KL state
                    }
                    else
                    {
                        PP_next_state = _CLR ;
                    }
            
                    break;
                }
            }
            //
            //
            ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
            //
            // Now generate s list of even and widely spaced x, y way points.. @ 30 metres apart
            //
            ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
            //
            
            vector<double> ptsx;
            vector<double> ptsy;
            
            // Refernce x, y, yaw states
            
            double ref_x    = car_x             ;
            double ref_y    = car_y             ;
            double ref_yaw  = deg2rad(car_yaw)  ;
            
            if (prev_size < 2)
            {
                // Use 2 points that make the path tangent to the car
                
                double prev_car_x = car_x - cos(car_yaw);
                double prev_car_y = car_y - sin(car_yaw);
                
                ptsx.push_back(prev_car_x);
                ptsx.push_back(car_x);
                
                ptsy.push_back(prev_car_y);
                ptsy.push_back(car_y);
                
            }
            // use the previous path's points as starting reference
            else
            {
                ref_x = previous_path_x[prev_size -1];
                ref_y = previous_path_y[prev_size -1];
                
                double ref_x_prev = previous_path_x[prev_size -2];
                double ref_y_prev = previous_path_y[prev_size -2];
                
                // Use 2 points that make the path tangent to the previous path's end points
                
                ptsx.push_back(ref_x_prev);
                ptsx.push_back(ref_x);
                
                ptsy.push_back(ref_y_prev);
                ptsy.push_back(ref_y);

            }
    
            vector<double> next_wp0 = getXY(car_s + 30, (2 + 4 * next_lane), map_waypoints_s, map_waypoints_x, map_waypoints_y);
            vector<double> next_wp1 = getXY(car_s + 60, (2 + 4 * next_lane), map_waypoints_s, map_waypoints_x, map_waypoints_y);
            vector<double> next_wp2 = getXY(car_s + 90, (2 + 4 * next_lane), map_waypoints_s, map_waypoints_x, map_waypoints_y);
            
            ptsx.push_back(next_wp0[0]);
            ptsx.push_back(next_wp1[0]);
            ptsx.push_back(next_wp2[0]);
            
            ptsy.push_back(next_wp0[1]);
            ptsy.push_back(next_wp1[1]);
            ptsy.push_back(next_wp2[1]);
            
            
            for (int i = 0; i < ptsx.size(); i++)
            {
                // Shift the car reference angle to 0 degrees.. simplifies the math..
                
                double shift_x = ptsx[i] - ref_x;
                double shift_y = ptsy[i] - ref_y;
                
                ptsx[i] = (shift_x * cos(0-ref_yaw) - shift_y * sin(0-ref_yaw));
                ptsy[i] = (shift_x * sin(0-ref_yaw) + shift_y * cos(0-ref_yaw));
                
            }
            
            // Create a spline path with the points generated
            
            tk::spline s ;
            
            s.set_points(ptsx, ptsy);
            
            
            // Define the actual points that will be used by the planner
            
            vector<double> next_x_vals;
            vector<double> next_y_vals;
            
            int previous_values = previous_path_x.size() ;
            
            // Start with all of the previous points from last time, not discarding..
            
            for (int i = 0; i < previous_path_x.size(); i++ )
            {
                next_x_vals.push_back(previous_path_x[i]);
                next_y_vals.push_back(previous_path_y[i]);

            }
            
            // Calculate how to break up spline points so that the car travels @ the desired speed
            
            double target_x     = 30;
            double target_y     = s(target_x);
            double target_dist  = sqrt((target_x * target_x) + (target_y * target_y));
            
            double x_add_on = 0;
            
            // Now fill up the rest of the points for our path planner after filling up with the previous points
            
            for (int i = 0; i <= 50 - previous_values; i++)
            {
                double N = (target_dist / (0.02 * ref_velocity/ 2.24)); // MPH to m/s
                double x_point = x_add_on + (target_x) / N;
                double y_point = s(x_point);
                
                x_add_on = x_point;
                
                double x_ref = x_point;
                double y_ref = y_point;
                
                // Rotate back to normal after rotating it earlier
                x_point = (x_ref * cos(ref_yaw) - y_ref*sin(ref_yaw));
                y_point = (x_ref * sin(ref_yaw) + y_ref*cos(ref_yaw));
                
                x_point += ref_x;
                y_point += ref_y;
                
                next_x_vals.push_back(x_point);
                next_y_vals.push_back(y_point);
            }
            
            //
            ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
            //
            // Load the simulator with the waypoints generated
            //
            ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
            //
            
          	json msgJson;
        
          	msgJson["next_x"] = next_x_vals;
          	msgJson["next_y"] = next_y_vals;

          	auto msg = "42[\"control\","+ msgJson.dump()+"]";

          	//this_thread::sleep_for(chrono::milliseconds(1000));
          	ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
            
            // Save the current state for the next run
            
            PP_current_state = PP_next_state;
          
        }
      } else {
        // Manual driving
        std::string msg = "42[\"manual\",{}]";
        ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
      }
    }
  });

  // We don't need this since we're not using HTTP but if it's removed the
  // program
  // doesn't compile :-(
  h.onHttpRequest([](uWS::HttpResponse *res, uWS::HttpRequest req, char *data,
                     size_t, size_t) {
    const std::string s = "<h1>Hello world!</h1>";
    if (req.getUrl().valueLength == 1) {
      res->end(s.data(), s.length());
    } else {
      // i guess this should be done more gracefully?
      res->end(nullptr, 0);
    }
  });

  h.onConnection([&h](uWS::WebSocket<uWS::SERVER> ws, uWS::HttpRequest req) {
    std::cout << "Connected!!!" << std::endl;
  });

  h.onDisconnection([&h](uWS::WebSocket<uWS::SERVER> ws, int code,
                         char *message, size_t length) {
    ws.close();
    std::cout << "Disconnected" << std::endl;
  });

  int port = 4567;
  if (h.listen(port)) {
    std::cout << "Listening to port " << port << std::endl;
  } else {
    std::cerr << "Failed to listen to port" << std::endl;
    return -1;
  }
  h.run();
}
















































































