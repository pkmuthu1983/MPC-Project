#include "MPC.h"
#include <cppad/cppad.hpp>
#include <cppad/ipopt/solve.hpp>
#include "Eigen-3.3/Eigen/Core"

using CppAD::AD;

typedef CPPAD_TESTVECTOR(double) Dvector;

double T = 0.75;
size_t N = 10;
double dt = T/(double) N;

// This value assumes the model presented in the classroom is used.
//
// It was obtained by measuring the radius formed by running the vehicle in the
// simulator around in a circle with a constant steering angle and velocity on a
// flat terrain.
//
// Lf was tuned until the the radius formed by the simulating the model
// presented in the classroom matched the previous radius.
//
// This is the length from front to CoG that has a similar radius.

int iter=0;
const double Lf = 2.67;
double ref_v = 60;
size_t x_start = 0;
size_t y_start = x_start + N;
size_t psi_start = y_start + N;
size_t v_start = psi_start + N;
size_t cte_start = v_start + N;
size_t epsi_start = cte_start + N;
size_t delta_start = epsi_start + N;
size_t a_start = delta_start + N - 1;

double polyeval(Eigen::VectorXd, double);

double gradx(Eigen::VectorXd coeffs, double x) {
	double res = 0;
	for (unsigned int i = 1; i < coeffs.size(); i++) {
		res += i * coeffs[i] * pow(x, i-1);
	}
	return (res);
}

void displaycosts(const Dvector& vars) {
    size_t t = 0;

    double cost_cte = 0, cost_epsi = 0;
    double cost_velocity = 0, cost_delta = 0, cost_a = 0, cost_delta_diff = 0, cost_a_diff = 0, cost_total = 0;
    double k_cte = 1500, k_epsi = 1500, k_v = 1, k_mag = 5, k_diff = 10;

    for (t = 0; t < N; t++) {
      cost_cte += k_cte * CppAD::pow(vars[cte_start + t], 2);
      cost_epsi += k_epsi * CppAD::pow(vars[epsi_start + t], 2);
      cost_velocity += k_v * CppAD::pow(vars[v_start + t] - ref_v, 2);
    }

    // Minimize the use of actuators.
    for (t = 0; t < N - 1; t++) {
      cost_delta += k_mag * CppAD::pow(vars[delta_start + t], 2);
      cost_a += k_mag * CppAD::pow(vars[a_start + t], 2);
    }

    // Minimize the value gap between sequential actuations.
    for (t = 0; t < N - 2; t++) {
      cost_delta_diff += k_diff * CppAD::pow(vars[delta_start + t + 1] - vars[delta_start + t], 2);
      cost_a_diff += k_diff * CppAD::pow(vars[a_start + t + 1] - vars[a_start + t], 2);
    }
    cost_total += cost_cte + cost_epsi + cost_velocity + cost_delta + cost_a + cost_delta_diff + cost_a_diff;

#ifdef DEBUG
    cout << "DEBUG costs " << cost_cte << "  " << cost_epsi << "  " << cost_velocity << "  " << cost_delta << "  " << cost_a << "  " << cost_delta_diff << "  " << cost_a_diff << " " << cost_total << endl;
#endif

}

class FG_eval {
 public:
  // Fitted polynomial coefficients
  Eigen::VectorXd coeffs;
  FG_eval(Eigen::VectorXd coeffs) { this->coeffs = coeffs; }

  typedef CPPAD_TESTVECTOR(AD<double>) ADvector;

  void operator()(ADvector& fg, const ADvector& vars) {
    size_t t;
    double k_cte = 1500, k_epsi = 1500, k_v = 1, k_mag = 5, k_diff = 10;

    fg[0] = 0;

    // The part of the cost based on the reference state.
    for (t = 0; t < N; t++) {
      fg[0] += k_cte * CppAD::pow(vars[cte_start + t], 2);
      fg[0] += k_epsi * CppAD::pow(vars[epsi_start + t], 2);
      fg[0] += k_v * CppAD::pow(vars[v_start + t] - ref_v, 2);
    }

    // Minimize the use of actuators.
    for (t = 0; t < N - 1; t++) {
      fg[0] += k_mag * CppAD::pow(vars[delta_start + t], 2);
      fg[0] += k_mag * CppAD::pow(vars[a_start + t], 2);
    }

    // Minimize the value gap between sequential actuations.
    for (t = 0; t < N - 2; t++) {
      fg[0] += k_diff * CppAD::pow(vars[delta_start + t + 1] - vars[delta_start + t], 2);
      fg[0] += k_diff * CppAD::pow(vars[a_start + t + 1] - vars[a_start + t], 2);
    }

    fg[1 + x_start] = vars[x_start];
    fg[1 + y_start] = vars[y_start];
    fg[1 + psi_start] = vars[psi_start];
    fg[1 + v_start] = vars[v_start];
    fg[1 + cte_start] = vars[cte_start];
    fg[1 + epsi_start] = vars[epsi_start];

    for (t = 1; t < N; t++) {
      AD<double> x1 = vars[x_start + t];
      AD<double> y1 = vars[y_start + t];
      AD<double> psi1 = vars[psi_start + t];
      AD<double> v1 = vars[v_start + t];
      AD<double> cte1 = vars[cte_start + t];
      AD<double> epsi1 = vars[epsi_start + t];

      AD<double> x0 = vars[x_start + t - 1];
      AD<double> y0 = vars[y_start + t - 1];
      AD<double> psi0 = vars[psi_start + t - 1];
      AD<double> v0 = vars[v_start + t - 1];
      AD<double> epsi0 = vars[epsi_start + t - 1];

      AD<double> delta0 = vars[delta_start + t - 1];
      AD<double> a0 = vars[a_start + t - 1];

      AD<double> f0 = coeffs[0] + coeffs[1] * x0 + coeffs[2] * x0 * x0 + coeffs[3] * x0 * x0 * x0;
      AD<double> gradx0 = coeffs[1] + 2 * coeffs[2] * x0 + 3 * coeffs[3] * x0 * x0;
      AD<double> psides0 = CppAD::atan(gradx0);

      fg[1 + x_start + t] = x1 - (x0 + v0 * CppAD::cos(psi0) * dt);
      fg[1 + y_start + t] = y1 - (y0 + v0 * CppAD::sin(psi0) * dt);
      fg[1 + psi_start + t] = psi1 - (psi0 - (v0/Lf) * delta0 * dt);
      fg[1 + v_start + t] = v1 - (v0 + a0 * dt);
      fg[1 + cte_start + t] =  cte1 - (f0 - y0 + v0 * CppAD::sin(epsi0) * dt);
      fg[1 + epsi_start + t] = epsi1 - ((psi0 - psides0) - (v0/Lf) * delta0 * dt);
  }
  }

};

//
// MPC class definition implementation.
//
MPC::MPC() {}
MPC::~MPC() {}

vector<double> MPC::Solve(Eigen::VectorXd state, Eigen::VectorXd coeffs,
    vector<double> &mpc_x_vals, vector<double> &mpc_y_vals) {
  bool ok = true;
  size_t i;

  // TODO: Set the number of model variables (includes both states and inputs).
  // For example: If the state is a 4 element vector, the actuators is a 2
  // element vector and there are 10 timesteps. The number of variables is:
  //
  // 4 * 10 + 2 * 9

  size_t n_vars = N * 6 + (N-1) * 2;
  size_t n_constraints = N * 6;

  // Initial value of the independent variables.
  // Should be 0 except for the initial values.
  Dvector vars(n_vars);
  for (i = 0; i < n_vars; i++) {
    vars[i] = 0.0;
  }
  // Set the initial variable values
  double x, y, psi, v, cte, epsi;
  x = state[0]; y = state[1]; psi = state[2]; v = state[3]; cte = state[4]; epsi = state[5];

  vars[x_start] = x;
  vars[y_start] = y;
  vars[psi_start] = psi;
  vars[v_start] = v;
  vars[cte_start] = cte;
  vars[epsi_start] = epsi;

  // Lower and upper limits for x
  Dvector vars_lowerbound(n_vars);
  Dvector vars_upperbound(n_vars);

  // Set all non-actuators upper and lowerlimits
  // to the max negative and positive values.
  for (i = 0; i < delta_start; i++) {
    vars_lowerbound[i] = -1.0e19;
    vars_upperbound[i] = 1.0e19;
  }

  // The upper and lower limits of delta are set to -25 and 25
  // degrees (values in radians).
  // NOTE: Feel free to change this to something else.
  for (i = delta_start; i < a_start; i++) {
    vars_lowerbound[i] = -0.436332 * Lf;
    vars_upperbound[i] = 0.436332 * Lf;
  }

  // Acceleration/decceleration upper and lower limits.
  // NOTE: Feel free to change this to something else.
  for (i = a_start; i < n_vars; i++) {
    vars_lowerbound[i] = -1.0;
    vars_upperbound[i] = 1.0;
  }

  // Lower and upper limits for the constraints
  // Should be 0 besides initial state.
  Dvector constraints_lowerbound(n_constraints);
  Dvector constraints_upperbound(n_constraints);

  for (i = 0; i < n_constraints; i++) {
    constraints_lowerbound[i] = 0;
    constraints_upperbound[i] = 0;
  }
  constraints_lowerbound[x_start] = x;
  constraints_lowerbound[y_start] = y;
  constraints_lowerbound[psi_start] = psi;
  constraints_lowerbound[v_start] = v;
  constraints_lowerbound[cte_start] = cte;
  constraints_lowerbound[epsi_start] = epsi;

  constraints_upperbound[x_start] = x;
  constraints_upperbound[y_start] = y;
  constraints_upperbound[psi_start] = psi;
  constraints_upperbound[v_start] = v;
  constraints_upperbound[cte_start] = cte;
  constraints_upperbound[epsi_start] = epsi;

  // object that computes objective and constraints
  FG_eval fg_eval(coeffs);

  //
  // NOTE: You don't have to worry about these options
  //
  // options for IPOPT solver
  std::string options;
  // Uncomment this if you'd like more print information
  options += "Integer print_level  0\n";
  // NOTE: Setting sparse to true allows the solver to take advantage
  // of sparse routines, this makes the computation MUCH FASTER. If you
  // can uncomment 1 of these and see if it makes a difference or not but
  // if you uncomment both the computation time should go up in orders of
  // magnitude.
  options += "Sparse  true        forward\n";
  options += "Sparse  true        reverse\n";
  // NOTE: Currently the solver has a maximum time limit of 0.5 seconds.
  // Change this as you see fit.
  options += "Numeric max_cpu_time          0.5\n";

  // place to return solution
  CppAD::ipopt::solve_result<Dvector> solution;

  // solve the problem
  CppAD::ipopt::solve<Dvector, FG_eval>(
      options, vars, vars_lowerbound, vars_upperbound, constraints_lowerbound,
      constraints_upperbound, fg_eval, solution);

  // Check some of the solution values
  ok &= solution.status == CppAD::ipopt::solve_result<Dvector>::success;

  // Cost
  auto cost = solution.obj_value;
  std::cout << "Cost " << cost << std::endl;

  // TODO: Return the first actuator values. The variables can be accessed with
  // `solution.x[i]`.
  //
  // {...} is shorthand for creating a vector, so auto x1 = {1.0,2.0}
  // creates a 2 element double vector.
  for (i = 1; i < N; i++){
	mpc_x_vals.push_back(solution.x[x_start+i]);
	mpc_y_vals.push_back(solution.x[y_start+i]);
  }

#ifdef DEBUG
    std::cout << "DEBUG Iteration "<< iter++ << std::endl;
    std::cout << "DEBUG x = " << solution.x[x_start] << std::endl;
    std::cout << "DEBUG y = " << solution.x[y_start] << std::endl;
    std::cout << "DEBUG psi = " << solution.x[psi_start] << std::endl;
    std::cout << "DEBUG v = " << solution.x[v_start] << std::endl;
    std::cout << "DEBUG cte = " << solution.x[cte_start] << std::endl;
    std::cout << "DEBUG epsi = " << solution.x[epsi_start] << std::endl;
    std::cout << "DEBUG delta = " << solution.x[delta_start] << std::endl;
    std::cout << "DEBUG a = " << solution.x[a_start] << std::endl;
    displaycosts(solution.x);
#endif

  return {solution.x[delta_start], solution.x[a_start]};
}
