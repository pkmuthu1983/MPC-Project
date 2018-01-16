# MPC controller Project

#### 1. Source code

'src' directory contains the implementation of the MPC controller (MPC.cpp), and
main.cpp. 'build' directory contains the executable (mpc). It can be run as
"./mpc".

#### 2. The Model:

The vehicle's state at time t, includes the (x,y) position, the heading
direction, and velocity along with the cross-track error and the heading
direction error. That is, state = (x(t), y(t), psi(t), v(t), cte(t), epsi(t)).
The vehicle is controlled via actuations (delta, a), where, 'delta' is the
steering control and 'a' is throttle.

Given the current state and current actuations, the vehicle's next state is
given by:

x(t+dt) = x(t) + v(t) * cos(psi(t)) * dt

y(t+dt) = y(t) + v(t) * sin(psi(t)) * dt

psi(t+dt) = psi(t) + v(t)/Lf * delta * dt

v(t+dt) = v(t) + a * dt

cte(t+dt) = f(x(t)) - y(t) + v(t) * sin(epsi(t)) * dt

epsi(t+dt) = (psi(t) - psi_des(t) + (v(t)/Lf * delta * dt

Here, f(.) is the polynomial which provides the reference trajectory.

The MPC planning horizon 'T' is divided into N timesteps of duration dt. Given a set
of waypoints near the vehicle's initial location, the actuations for each of the
N timesteps is obtained by solving an optimization problem, such that the
vehicle's path remains close to the trajectory given by waypoints, and the
vehicles velocity is close to the reference velocity (ref_v).

The objective function for the optimization problem is defined in src/MPC.cpp
'FG_eval' class. 

#### 3. Time Length and Duration (T, N, dt):

The simulator provides the locations of 6 nearby
waypoints.  First, the time horizon 'T', is chosen such that the MPC
controller does not look beyond (ie extrapolate) the given waypoints. This of
course depends on the vehicles veloctiy (ref_v). So, once a reference velocity
(60mph) is set, I experimented with few choices of T by plotting the waypoints,
and the MPC predicted points. When T was 0.75 seconds or less, the MPC
controller's planning horizon does not need to extrapolate the given waypoints.

After T is chosen, N was chosen based on experimentation as well. If N is small (<
8), then the last few MPC predicted points was not close to the fitted waypoint trajectory due to
discretization error, but the vehicle could steer. If N was > 12, then the
actuations were not optimal, and the vehicle went off the track.  However, a
choice of N = 10, worked just fine. Note that dt is automatically set to
T/N, once T and N are decided.

So, T = 0.75, N = 10, and dt = T/N worked fine. Refer src/MPC.cpp, lines 10-13.

I also tried the following choices:

1) T = 0.5 seconds, and N = 8, 10, 15: With these choices, the vehicle started
oscillating beyond a point.

2) T = 1 seconds, and N = 8, 10, 15: With this, vehicle was able to steer
succesfully. But, near curves, the controller looked beyond the given the
waypoints, therefore the MPC predicted points were little off.

3) T = 1.5 seconds, and N = 8, 10, 15: Again, this did not work very well due to
extrapolation.

#### 4. Polynomial Fitting and MPC Preprocessing:

The simulator provides the car's location (x,y), six waypoints, the car heading
direction (psi), current steering, and velocity. Note that the locations are in
map coordinate system.

The waypoints are first translated into car-coordinate system (see
src/main.cpp, line 139) centered at (x,y) and oriented with respect to the map
coordinate system by psi. Then a third order polynomial is fitted to these
translated waypoints. This polynomial is provided as input to the MPC solver for
calculating the optimal actuations.

Due to this transaltion, the initial location of the car (x(0), y(0)) is (0, 0).
The cte(0), and epsi(0) are calculated using the fitted trajectory.

cte(0) = f(x(0)) - y(0)

epsi(0) =  psi - atan(grad_f(x(0)))

#### 5. Model Predictive Control with Latency:

In order to incorporate actuation latency, the initial state of the car presented to MPC optimizer
is adjusted based on the latency. That is, instead of the above initial state,
the car's state after the latency, is estimated based on the update equations
in section 2, and then given to MPC optimizer.

Refer src/main.cpp, line 156, for the code which updates the initial state.
