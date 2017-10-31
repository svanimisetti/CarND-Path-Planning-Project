#include <fstream>
#include <math.h>
#include <uWS/uWS.h>
#include <chrono>
#include <iostream>
#include <thread>
#include <vector>
#include <cppad/cppad.hpp>
#include <cppad/ipopt/solve.hpp>
#include "Eigen-3.3/Eigen/Core"
#include "Eigen-3.3/Eigen/QR"
#include "json.hpp"
#include "spline.h"

using namespace std;

// for convenience
using json = nlohmann::json;

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

double distance(double x1, double y1, double x2, double y2) {
  return sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));
}

int ClosestWaypoint(double x, double y, vector<double> maps_x, vector<double> maps_y) {
  double closestLen = 100000; //large number
  int closestWaypoint = 0;
  for(int i = 0; i < maps_x.size(); i++) {
    double map_x = maps_x[i]; double map_y = maps_y[i];
    double dist = distance(x,y,map_x,map_y);
    if(dist < closestLen) { closestLen = dist; closestWaypoint = i; }
  }
  return closestWaypoint;
}

int NextWaypoint(double x, double y, double theta, vector<double> maps_x, vector<double> maps_y) {
  int closestWaypoint = ClosestWaypoint(x,y,maps_x,maps_y);
  double map_x = maps_x[closestWaypoint]; double map_y = maps_y[closestWaypoint];
  double heading = atan2( (map_y-y),(map_x-x) );
  double angle = abs(theta-heading);
  if(angle > pi()/4) {
    closestWaypoint++;
    if (closestWaypoint >= maps_x.size()) closestWaypoint = 0;
  }
  return closestWaypoint;
}

// Transform from Cartesian x,y coordinates to Frenet s,d coordinates
vector<double> getFrenet(double x, double y, double theta, vector<double> maps_x, vector<double> maps_y) {
  int next_wp = NextWaypoint(x,y, theta, maps_x,maps_y);

  int prev_wp;
  prev_wp = next_wp-1;
  if(next_wp == 0) {
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

  if(centerToPos <= centerToRef) {
    frenet_d *= -1;
  }

  // calculate s value
  double frenet_s = 0;
  for(int i = 0; i < prev_wp; i++) {
    frenet_s += distance(maps_x[i],maps_y[i],maps_x[i+1],maps_y[i+1]);
  }

  frenet_s += distance(0,0,proj_x,proj_y);

  return {frenet_s,frenet_d};
}

// !========== Helper Function ==========! //

// Coordinate transformation helper functions
// 1. global (world) to local (car)
vector<double> locXY(double cx, double cy, double ang, double wpx, double wpy) {
  vector<double> loc;
  loc.push_back( (wpx-cx)*cos(ang)+(wpy-cy)*sin(ang));
  loc.push_back(-(wpx-cx)*sin(ang)+(wpy-cy)*cos(ang));
  return loc;
}
// 2. local (car) to global (world)
vector<double> glbXY(double cx, double cy, double ang, double lx, double ly) {
  vector<double> glb;
  glb.push_back(lx*cos(ang)-ly*sin(ang)+cx);
  glb.push_back(lx*sin(ang)+ly*cos(ang)+cy);
  return glb;
}
// 3. Conversion of local WP segments
vector<vector<double>> locWPSeg(double cx, double cy, double yaw, double d,
							    vector<double> mx, vector<double> my,
							    vector<double> mdx, vector<double> mdy) {
  vector<double> wp_x, wp_y; vector<vector<double>> loc_seg;
  int prev_wp = ClosestWaypoint(cx, cy, mx, my)-6;
  if (prev_wp < 0) { prev_wp += mx.size(); }
  for (int i = 0; i < 25; i++) {
    int next_wp = (prev_wp+i)%mx.size();
    vector<double> lxy = locXY(cx, cy, deg2rad(yaw),
	                           (mx[next_wp]+d*mdx[next_wp]),
							   (my[next_wp]+d*mdy[next_wp]));
    wp_x.push_back(lxy[0]); wp_y.push_back(lxy[1]);
  }
  loc_seg.push_back(wp_x); loc_seg.push_back(wp_y);
  return loc_seg;
}
// 4. Convert set of global coords to local coords
vector<vector<double>> g2l_pts(double cx, double cy, double yaw,
                               vector<double> gx, vector<double> gy) {
  vector<double> lx, ly; vector<vector<double>> lpts;
  for (int i = 0; i < gx.size(); i++) {
    vector<double> lxy = locXY(cx, cy, deg2rad(yaw), gx[i], gy[i]);
    lx.push_back(lxy[0]); ly.push_back(lxy[1]);
  }
  lpts.push_back(lx); lpts.push_back(ly);
  return lpts;
}
// 5. Convert a set of local x,y vector coordinates to world x y vectors
vector<vector<double>> l2g_pts(double cx, double cy, double yaw,
                               vector<double> lx, vector<double> ly) {
  vector<double> gx, gy; vector<vector<double>> gpts;
  for (int i = 0; i < lx.size(); i++) {
    vector<double> gxy = glbXY(cx, cy, deg2rad(yaw), lx[i], ly[i]);
    gx.push_back(gxy[0]); gy.push_back(gxy[1]);
  }
  gpts.push_back(gx); gpts.push_back(gy);
  return gpts;
}

// !========== End of Helper Function ==========! //

int main() {
  uWS::Hub h;

  // Load up map values for waypoint's x,y and d normalized normal vectors
  vector<double> map_waypoints_x;
  vector<double> map_waypoints_y;
  vector<double> map_waypoints_dx;
  vector<double> map_waypoints_dy;

  // Waypoint map to read from
  string map_file_ = "../data/highway_map.csv";
  // The max s value before wrapping around the track back to 0
  double max_s = 6945.554;
  
  // !===== Var, Const and Param definitions =====! //
  vector<double> vx, vy, vd, xyd;
  const double inc_max = 0.445;
  double nextd = 6., dist_inc = inc_max;
  int timestep = 0, watchdog_timer = 0;
  bool lanechange = false;

  ifstream in_map_(map_file_.c_str(), ifstream::in);

  string line;
  while (getline(in_map_, line)) {
    istringstream iss(line);
    double x, y, s, d_x, d_y;
    iss >> x; iss >> y; iss >> s; iss >> d_x; iss >> d_y;
    map_waypoints_x.push_back(x); map_waypoints_y.push_back(y);
    map_waypoints_dx.push_back(d_x); map_waypoints_dy.push_back(d_y);
  }

  h.onMessage([&map_waypoints_x,&map_waypoints_y,&map_waypoints_dx,&map_waypoints_dy,&vx,&vy,&vd,&xyd,&nextd,&inc_max,&dist_inc,&timestep,&watchdog_timer,&lanechange]
              (uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length, uWS::OpCode opCode) {
    // "42" at the start of the message means there's a websocket message event.
    // The 4 signifies a websocket message
    // The 2 signifies a websocket event
    //auto sdata = string(data).substr(0, length);
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
          double angle = deg2rad(car_yaw);

          // Previous path data given to the Planner
          auto previous_path_x = j[1]["previous_path_x"];
          auto previous_path_y = j[1]["previous_path_y"];
          // Previous path's end s and d values 
          double end_path_s = j[1]["end_path_s"];
          double end_path_d = j[1]["end_path_d"];

          // Sensor Fusion Data, a list of all other cars on the same side of the road.
          auto sensor_fusion = j[1]["sensor_fusion"];

          json msgJson;

          vector<double> next_x_vals;
          vector<double> next_y_vals;
		  
		  // !===== Begin new path definition =====! //
		  
		  vector<double> lx, ly; // for lane and speed control
          int n_path = previous_path_x.size(), npts = 50;
          tk::spline lane_sm, vel_sm, local_sm;

          // get frenet and spline for lane points
          vector<double> frenet = getFrenet(car_x, car_y, deg2rad(car_yaw),
                                            map_waypoints_x, map_waypoints_y);
          vector<vector<double>> lwpxy = locWPSeg(car_x, car_y, car_yaw, nextd,
                                                  map_waypoints_x, map_waypoints_y,
                                                  map_waypoints_dx, map_waypoints_dy);
          // Perform basic checks
		  // 1. Wrong way
		  if (lwpxy[0][0] > 0.) {
            car_yaw += 180; angle = deg2rad(car_yaw);
            lwpxy = locWPSeg(car_x, car_y, car_yaw, nextd,
                             map_waypoints_x,  map_waypoints_y,
		                     map_waypoints_dx, map_waypoints_dy);
          }
		  lane_sm.set_points(lwpxy[0], lwpxy[1]);
		  // 2. If there is no path, process using following
		  if (n_path == 0) {
			vector<double> t, inc;
            t.insert(t.end(),{-1.0, 15.0, 25.0, double(npts), double(npts*2)});
			inc.insert(inc.end(), {dist_inc*0.01, dist_inc*0.10, dist_inc*0.25,
			                       dist_inc*0.65, dist_inc*1.0});
            vel_sm.set_points(t,inc);
			
            double lwpx_nxt = 0.0, lwpy_nxt = 0.0;
            for (int i = 0; i<npts; i++) {
              lwpx_nxt += vel_sm(double(i)); lwpy_nxt = lane_sm(lwpx_nxt);
              lx.push_back(lwpx_nxt); ly.push_back(lwpy_nxt);
              if (i > 0)
                vd.push_back(distance(lx[i-1], ly[i-1], lx[i], ly[i]));
              else
                vd.push_back(vel_sm(double(0)));
            }
			
            // Smoothen the path
            double localxx = 0., localxy = 0.;
            for(int i = 0; i < npts; i++) {
              ly[i] = lane_sm(lx[i]);
              double dist = distance(localxx, localxy, lx[i], ly[i]);
              double speed = vel_sm(double(i));
              if (dist > speed || dist < speed*0.8) {
                 double heading = atan2(ly[i]-localxy,lx[i]-localxx);
                 lx[i] = localxx + vel_sm(double(i))*cos(heading);
                 ly[i] = lane_sm(lx[i]);
                 dist = distance(localxx, localxy, lx[i], ly[i]);
              }
              localxx = lx[i]; localxy = ly[i];
            }

			// Convert smoothened path from local to global coordinates
            vector<vector<double>> worldxy = l2g_pts(car_x, car_y, car_yaw, lx, ly);
            for (int i=n_path; i<worldxy[0].size(); i++) {
              vx.push_back(worldxy[0][i]); vy.push_back(worldxy[1][i]);
              next_x_vals.push_back(worldxy[0][i]); next_y_vals.push_back(worldxy[1][i]);
            }

		  // If path already exists (n_path not 0), then car is moving and
		  // following code will attempt to update it using data from sensor_fusion
		  // and use path planner to determine new path
          } else {
			vector<vector<double>> previous_localxy = g2l_pts(car_x, car_y, car_yaw, previous_path_x, previous_path_y);
            lx = previous_localxy[0]; ly = previous_localxy[1];
			
            for (int i = 0; i < (npts-n_path); i++) {
              vector<double> frenet = getFrenet(vx[i], vy[i], deg2rad(car_yaw), map_waypoints_x, map_waypoints_y);
            }
            vx.erase(vx.begin(),vx.begin()+(npts-n_path));
            vy.erase(vy.begin(),vy.begin()+(npts-n_path));
            vd.erase(vd.begin(),vd.begin()+(npts-n_path));
            xyd.erase(xyd.begin(),xyd.begin()+(npts-n_path));

            if (lanechange && abs(lane_sm(0.)) < 0.01) {
              lanechange = false;
            }

            vector<vector<double>> newwxy = locWPSeg(car_x, car_y, car_yaw, nextd, map_waypoints_x, map_waypoints_y, map_waypoints_dx, map_waypoints_dy);
            if (newwxy[0][0] > 0.) {
              car_yaw += 180; angle = deg2rad(car_yaw);
              newwxy = locWPSeg(car_x, car_y, car_yaw, nextd, map_waypoints_x, map_waypoints_y, map_waypoints_dx, map_waypoints_dy);
            }
            tk::spline newlane;
			newlane.set_points(newwxy[0], newwxy[1]);
			vector<double> localwx; vector<double> localwy;
			for (int i=0; i<n_path; i++) {
              localwx.push_back(lx[i]); localwy.push_back(ly[i]);
            }
			double nextx = lx[n_path-1]+40;
			for (int i=0; i<n_path; i++) {
              localwx.push_back(nextx);
			  localwy.push_back(newlane(nextx));
              nextx += dist_inc;
            }
			lane_sm.set_points(localwx, localwy);
			
            for(int i = 0; i < n_path; i++) {
              next_x_vals.push_back(previous_path_x[i]);
              next_y_vals.push_back(previous_path_y[i]);
            }

			vector<double> t, inc;
            t.insert(t.end(),{0.0,125.0, 250.0});
			if (vd[0] < inc_max) {
			  inc.insert(inc.end(), {vd[0], 0.5*(vd[0]+dist_inc), dist_inc});
            } else {
              inc.insert(inc.end(), {vd[n_path-1], 0.5*(vd[n_path-1]+dist_inc), dist_inc});
            }
			vel_sm.set_points(t,inc);
			
            for(int i = n_path; i<npts; i++) {
              lx.push_back(lx[i-1]+vel_sm(double(i)));
			  ly.push_back(lane_sm(lx[i]));
              vx.push_back(0.0); vy.push_back(0.0);
              next_x_vals.push_back(0.0); next_y_vals.push_back(0.0);
            }
			
            double localxx = lx[0]; double localxy = ly[0];
            for(int i = 0; i < npts; i++) {
              ly[i] = lane_sm(lx[i]);
              double dist = distance(localxx, localxy, lx[i], ly[i]);
              if (dist > vel_sm(double(i))) {
                double heading = atan2(ly[i]-localxy,lx[i]-localxx);
                lx[i] = localxx + vel_sm(double(i))*cos(heading);
                ly[i] = lane_sm(lx[i]);
                dist = distance(localxx, localxy, lx[i], ly[i]);
              }
              if (i >= n_path) {
			    vd.push_back(dist);
			  } 
              localxx = lx[i]; localxy = ly[i];
            }

            vector<vector<double>> worldxy = l2g_pts(car_x, car_y, car_yaw, lx, ly);
            for (int i=n_path; i<worldxy[0].size(); i++) {
              vx[i] = worldxy[0][i]; vy[i] = worldxy[1][i];
              next_x_vals[i] = worldxy[0][i]; next_y_vals[i] = worldxy[1][i];
            }

          }

          // Following is the path planner
		  // Initialize individual lane data and append to "lanes" vector
          vector<double> l1, l2, l3;
          vector<vector<double>> lanes;
          int ourlane = round(round(nextd-2)/4);
          int bestlane = ourlane;
		  lanes.insert(lanes.end(),{l1,l2,l3});
          bool badsensor = false;
          int backvehicle_shift = 5;

		  // Process sensor fusion data to get location of cars in 3 lanes
		  // Retrieve data for each vehicles ID (vid) provided by sensor fusion
          for (int k = 0; k<sensor_fusion.size(); k++) {
            vector<double> vid = sensor_fusion[k];
            double vid_x = vid[1]+vid[3]*0.02;
            double vid_y = vid[2]+vid[4]*0.02;
            vector<double> vidlocal = locXY(car_x, car_y, deg2rad(car_yaw), vid_x, vid_y);
            double vid_dst = distance(car_x, car_y, vid[1], vid[2]);
            double vid_s = vidlocal[0] + backvehicle_shift;
            double vid_d = vid[6];
            sensor_fusion[k].push_back(vid_s);
            sensor_fusion[k].push_back(distance(0,0,vid[3],vid[4])*0.02);
            sensor_fusion[k].push_back(round(round(vid_d-2)/4));

            if (vel_sm(0.) < dist_inc) {
              if (vel_sm(0.) < vid[8]) vid_s += (vid[8]-vel_sm(0.))*120.;
            } else {
              if (dist_inc < vid[8]) vid_s += (vid[8]-dist_inc)*240.;
            }

			// Identify lane location for each vehicle ID
			// Try to identify if vehicle is in the middle of lane transition
            string lanestr = "error";
            if (vid_s > 0.) {
              if (vid_d < 12. && vid_d > 0.) {
				// inside lane 0
                if (vid_d <= 3.7) {
                  lanestr = "0"; lanes[0].push_back(vid_s);
                }
				// transition between lane 0 and 1
                if (vid_d > 3.7 && vid_d <= 4.3) {
                  lanestr = "0,1"; lanes[0].push_back(vid_s); lanes[1].push_back(vid_s);
                }
				// inside lane 1
                if (vid_d > 4.3 && vid_d <= 7.7) {
                  lanestr = "1"; lanes[1].push_back(vid_s);
                }
				// transition between lane 1 and 2
                if (vid_d > 7.7 && vid_d <= 8.3) {
                  lanestr = "1,2"; lanes[1].push_back(vid_s); lanes[2].push_back(vid_s);
                }
				// inside lane 2
                if (vid_d > 8.3) {
                  lanestr = "2"; lanes[2].push_back(vid_s);
                }
              } else {
                badsensor = true;
              }
            }
          }

          for (int lane = 0; lane<3; lane++) {
            if (lanes[lane].size() > 0) sort(lanes[lane].begin(),lanes[lane].end());
          }

		  // In case the path planner times out, restart by setting newlane
          if (watchdog_timer > 1000) {
            int newlane = -1;
            for (int lane=0; lane<lanes.size(); lane++) {
              if (newlane < 0 && ourlane != lane && abs(ourlane-lane)==1) {
                newlane = lane;
                for (int i=0; i<sensor_fusion.size(); i++) {
                  vector<double> oldvid = sensor_fusion[i];
                  vector<double> vid;
                  if (oldvid[7] == lanes[newlane][0] && (lanes[newlane][0] > 15 || lanes[newlane][0] < 5)) {
                    lanes[ourlane][0] = lanes[newlane][0]-1.;
                    for (int j=0; j<oldvid.size(); j++) {
                      vid.push_back(oldvid[j]);
                    }
                    vid[7] = lanes[ourlane][0];
                    vid[9] = double(ourlane);
                    sensor_fusion.push_back(vid);
                  } else {
                    watchdog_timer = 0;
          }}}}}

		  // Identify best lane to transition
          for (int lane = 0; lane<3; lane++) {
            if (lanes[lane].size() > 0) {
              if (lanes[bestlane].size() > 0 && (lanes[bestlane][0] < lanes[lane][0])) {
                if (lanes[ourlane].size() > 0 && lanes[ourlane][0] < 80. && abs(ourlane-lane)==1) {
                  if (abs(ourlane-lane) == 1) {
                    if (lanes[lane][0] > 20) {
                      bestlane = lane;
                    } else {
                      if (lanes[lane][0] > 5 && watchdog_timer > 1000) {
                        bestlane = lane;
              }}}}}
            } else {
              if (lanes[ourlane].size() > 0 && lanes[ourlane][0] < 80. && lanes[bestlane].size() > 0 && abs(ourlane-lane)==1) {
                if (abs(ourlane-lane) == 1) {
                  bestlane = lane;
                }
                if (dist_inc < inc_max) {
                  dist_inc = vd[n_path-1];
            }}}
          }
          int l0size = lanes[0].size();
          int l1size = lanes[1].size();
          int l2size = lanes[2].size();
          float l0closest = 0, l1closest = 0, l2closest = 0;
          if (l0size > 0) l0closest = lanes[0][0];
          if (l1size > 0) l1closest = lanes[1][0];
          if (l2size > 0) l2closest = lanes[2][0];

		  // If current lane is not best lane and sensor data is good,
		  // then plan for the lane change now
          if (timestep > 50 && ourlane != bestlane) {
            if ( not lanechange and not badsensor ) {
              nextd = bestlane*4+2;
              if (nextd > 7) nextd -= 0.3;
              lanechange = true; watchdog_timer = 0;
            }
          }

		  // If there are cars in current lane and best lane is current lane,
		  // then slow down to match speed of car ahead
          if (lanes[ourlane].size() > 0 && lanes[ourlane][0] < 50. && lanes[ourlane][0] > 5.) {
            if (ourlane != 1) watchdog_timer++;
            for (int i=0; i<sensor_fusion.size(); i++) {
              vector<double> vid = sensor_fusion[i];
              if (vid[7] == lanes[ourlane][0] && vid[9] == ourlane && (lanes[ourlane][0] < 40 || (lanes[ourlane][0] < 50 && vid[8] < 0.3))) {
                if (vid[8] > 0.25) {
                  if (dist_inc >= vid[8]) {
                    dist_inc = vid[8]*0.95;
                  } else {
                    if (vid[7] > 15) {
                      dist_inc = dist_inc + (vid[8] - dist_inc)/2.;
                    }
                  }
                } else {
                  dist_inc = 0.1;
          }}}} else {
            // slowly increase speed to avoid max acceleration...
            if (dist_inc < inc_max && not lanechange) {
              double interval = inc_max - vel_sm(0.);
              if (interval > 0.3) {
                double inc1 = (vel_sm(0.)*9. + inc_max)/10.;
                double inc2 = (dist_inc*9. + inc_max)/10.;
                if (inc1 > dist_inc && inc1 < inc2) {
                  dist_inc = inc1;
                } else {
                  dist_inc = inc2;
                }
              } else if (interval > 0.2) {
                double inc1 = (vel_sm(0.)*7. + inc_max)/8.;
                double inc2 = (dist_inc*7. + inc_max)/8.;
                if (inc1 > dist_inc && inc1 < inc2) {
                  dist_inc = inc1;
                } else {
                  dist_inc = inc2;
                }
              } else {
                double inc1 = (vel_sm(0.)*5. + inc_max)/6.;
                double inc2 = (dist_inc*5. + inc_max)/6.;
                if (inc1 > dist_inc && inc1 < inc2) {
                  dist_inc = inc1;
                } else {
                  dist_inc = inc2;
          }}}}
		  
		  // Append cleaned up path for sending back to the simulator
          vector<double> localx(next_x_vals.size());
          vector<double> localy(next_x_vals.size());
          for (int i=0; i < next_x_vals.size(); i++) {
            float next_x = (next_x_vals[i] - car_x);
            float next_y = (next_y_vals[i] - car_y);
            localx[i] = next_x*cos(angle) + next_y*sin(angle);
            localy[i] = -next_x*sin(angle) + next_y*cos(angle);
          }
		  local_sm.set_points(localx, localy);
		  double localxx = 0., localxy = 0.;
          for(int i = 0; i < npts; i++) {
            localy[i] = local_sm(localx[i]);
            double dist = distance(localxx, localxy, localx[i], localy[i]);
            if (dist > vel_sm(double(i))) {
               double heading = atan2(localy[i]-localxy,localx[i]-localxx);
               localx[i] = localxx + vel_sm(double(i))*cos(heading);
               localy[i] = local_sm(localx[i]);
               dist = distance(localxx, localxy, localx[i], localy[i]);
            }
            localxx = localx[i]; localxy = localy[i];
          }
		  
          for (int i=0; i<npts; i++) {
            next_x_vals[i] = localx[i]*cos(angle) - localy[i]*sin(angle) + car_x;
            next_y_vals[i] = localx[i]*sin(angle) + localy[i]*cos(angle) + car_y;
          }

          for (int i=n_path; i<npts; i++) {
            if (i > 0)
              xyd.push_back(distance(next_x_vals[i-1], next_y_vals[i-1], next_x_vals[i], next_y_vals[i]));
            else
              xyd.push_back(dist_inc);
          }
		  
		  // !===== End path definition =====! //
		  
		  msgJson["next_x"] = next_x_vals;
          msgJson["next_y"] = next_y_vals;

          auto msg = "42[\"control\","+ msgJson.dump()+"]";

          //this_thread::sleep_for(chrono::milliseconds(1000));
          //this_thread::sleep_for(chrono::milliseconds(500));
          ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
          timestep++;        
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

