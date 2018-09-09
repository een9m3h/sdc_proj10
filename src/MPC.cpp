#include "MPC.h"
#include <cppad/cppad.hpp>
#include <cppad/ipopt/solve.hpp>
#include "Eigen-3.3/Eigen/Core"

using CppAD::AD;

// TODO: Set the timestep length and duration
size_t N = 10;
double dt = 0.2;  // sampling period increased for act latency

size_t x_start 		= 0;
size_t y_start 		= x_start + N;
size_t psi_start 	= y_start + N;
size_t v_start 		= psi_start + N;
size_t cte_start 	= v_start + N;
size_t epsi_start 	= cte_start + N;
size_t delta_start  = epsi_start + N;
size_t a_start 		= delta_start + N - 1;

double cte_weight = 100; // weight for not having a low cross track error
double epsi_weight = 100; // for having an angle error
double speed_weight = 1; // not following the speed limit
double steer_use_weight = 1; // for steering the car
double a_use_weight = 1; // using the throttle
double steer_change_weight = 10; // having sharp / large steer angles between steps
double a_change_weight = 10; // accelerating or braking fast
double speed_vs_steer_weight = 90;

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
const double Lf = 2.67;

// Both the reference cross track and orientation errors are 0.
// The reference velocity is set to 40 mph.
double ref_v = 40.0;

class FG_eval {
 public:
  // Fitted polynomial coefficients
  Eigen::VectorXd coeffs;
  FG_eval(Eigen::VectorXd coeffs) { this->coeffs = coeffs; }

  typedef CPPAD_TESTVECTOR(AD<double>) ADvector;
  void operator()(ADvector& fg, const ADvector& vars) {
    // TODO: implement MPC
    // `fg` a vector of the cost constraints, `vars` is a vector of variable values (state & actuators)
    // NOTE: You'll probably go back and forth between this function and
    // the Solver function below.
	
	fg[0] = 0;
	
	// The part of the cost based on the reference state.
    for (int t = 0; t < N; t++) {
		//cross-track error
		fg[0] += cte_weight*CppAD::pow(vars[cte_start + t], 2);
		//heading error
		fg[0] += epsi_weight*CppAD::pow(vars[epsi_start + t], 2);
		//velocity error
		fg[0] += speed_weight*CppAD::pow(vars[v_start + t] - ref_v, 2);
    }

    // Minimize the use of actuators - constrain erratic control inputs.
    for (int t = 0; t < N - 1; t++) {
		//minimize wheel turn
		fg[0] += steer_use_weight*CppAD::pow(vars[delta_start + t], 2);
		// minimize acceleration
		fg[0] += a_use_weight*CppAD::pow(vars[a_start + t], 2);
		// tradeoff speed vs steer
		fg[0] += speed_vs_steer_weight*CppAD::pow(vars[delta_start + t] * vars[v_start+t], 2);
    }

    /* Minimize the value gap between sequential actuations.
	The goal of this final loop is to make control decisions more consistent, or smoother.
	The next control input should be similar to the current one.
	*/
    for (int t = 0; t < N - 2; t++) {
		//minimize wheel turn deltas
		fg[0] += steer_change_weight*CppAD::pow(vars[delta_start + t + 1] - vars[delta_start + t], 2);
		// minimize acceleration deltas
		fg[0] += a_change_weight*CppAD::pow(vars[a_start + t + 1] - vars[a_start + t], 2);
    }
	
	//initialize the model to the initial state (offset by 1 for cost)
	fg[1 + x_start] 	= vars[x_start];
	fg[1 + y_start] 	= vars[y_start];
	fg[1 + psi_start] 	= vars[psi_start];
	fg[1 + v_start] 	= vars[v_start];
	fg[1 + cte_start] 	= vars[cte_start];
	fg[1 + epsi_start] 	= vars[epsi_start];
	
	for (int t = 1; t < N; t++) {
		// The state at time t+1 .
		AD<double> x1 = vars[x_start + t];
		AD<double> y1 = vars[y_start + t];
		AD<double> psi1 = vars[psi_start + t];
		AD<double> v1 = vars[v_start + t];
		AD<double> cte1 = vars[cte_start + t];
		AD<double> epsi1 = vars[epsi_start + t];

		// The state at time t.
		AD<double> x0 = vars[x_start + t - 1];
		AD<double> y0 = vars[y_start + t - 1];
		AD<double> psi0 = vars[psi_start + t - 1];
		AD<double> v0 = vars[v_start + t - 1];
		AD<double> cte0 = vars[cte_start + t - 1];
		AD<double> epsi0 = vars[epsi_start + t - 1];

		// Only consider the actuation at time t.
		AD<double> delta0 = vars[delta_start + t - 1];
		AD<double> a0 = vars[a_start + t - 1];
		
		AD<double> f0 = coeffs[0] + coeffs[1] * x0 + coeffs[2] * CppAD::pow(x0, 2) + coeffs[3] * CppAD::pow(x0, 3);
		AD<double> psides0 = CppAD::atan(coeffs[1] + 2 * coeffs[2] * x0 + 3 * coeffs[3] * CppAD::pow(x0, 2));

		// Here's `x` to get you started.
		// The idea here is to constraint this value to be 0.
		//
		// Recall the equations for the model:
		// x_[t] = x[t-1] + v[t-1] * cos(psi[t-1]) * dt
		// y_[t] = y[t-1] + v[t-1] * sin(psi[t-1]) * dt
		// psi_[t] = psi[t-1] + v[t-1] / Lf * delta[t-1] * dt
		// v_[t] = v[t-1] + a[t-1] * dt
		// cte[t] = f(x[t-1]) - y[t-1] + v[t-1] * sin(epsi[t-1]) * dt
		// epsi[t] = psi[t] - psides[t-1] + v[t-1] * delta[t-1] / Lf * dt
		fg[1 + x_start + t] = x1 - (x0 + v0 * CppAD::cos(psi0) * dt);
		fg[1 + y_start + t] = y1 - (y0 + v0 * CppAD::sin(psi0) * dt);
		fg[1 + psi_start + t] = psi1 - (psi0 - v0 * delta0 / Lf * dt);
		fg[1 + v_start + t] = v1 - (v0 + a0 * dt);
		fg[1 + cte_start + t] =
		  cte1 - ((f0 - y0) + (v0 * CppAD::sin(epsi0) * dt));
		fg[1 + epsi_start + t] =
		  epsi1 - ((psi0 - psides0) - ((v0 * delta0) / Lf) * dt);
	}
	
  }
};

//
// MPC class definition implementation.
//
MPC::MPC() {}
MPC::~MPC() {}

vector<double> MPC::Solve(Eigen::VectorXd x0, Eigen::VectorXd coeffs) {
	bool ok = true;
	size_t i;
	typedef CPPAD_TESTVECTOR(double) Dvector;

	// TODO: Set the number of model variables (includes both states and inputs).
	// For example: If the state is a 4 element vector, the actuators is a 2
	// element vector and there are 10 timesteps. The number of variables is:
	//
	// 4 * 10 + 2 * 9
	size_t n_vars = 6 * N + 2 * (N-1);
	// Number of constraints
	size_t n_constraints = N * 6;

	// Initial value of the independent variables.
	// SHOULD BE 0 besides initial state.
	Dvector vars(n_vars);
	for (int i = 0; i < n_vars; i++) {
		vars[i] = 0;
	}
  
	double x 	= x0[0];
	double y 	= x0[1];
	double psi 	= x0[2];
	double v 	= x0[3];
	double cte 	= x0[4];
	double epsi = x0[5];
	
	std::cout << "Set the initial variable values" << std::endl;

	// Set the initial variable values
	vars[x_start] 		= x;
	vars[y_start] 		= y;
	vars[psi_start] 	= psi;
	vars[v_start] 		= v;
	vars[cte_start] 	= cte;
	vars[epsi_start] 	= epsi;

	Dvector vars_lowerbound(n_vars);
	Dvector vars_upperbound(n_vars);
	// Set all non-actuators upper and lowerlimits
  
	std::cout << "Svars_lowerbound[i] = -1.0e19;" << std::endl;
  
	// to the max negative and positive values.
	for (int i = 0; i < delta_start; i++) {
		vars_lowerbound[i] = -1.0e19;
		vars_upperbound[i] = 1.0e19;
	}

	std::cout << "vars_lowerbound[i] = -0.436332;" << std::endl;
	
	// The upper and lower limits of delta are set to -25 and 25
	// degrees (values in radians).
	// NOTE: Feel free to change this to something else.
	for (int i = delta_start; i < a_start; i++) {
		vars_lowerbound[i] = -0.436332;
		vars_upperbound[i] = 0.436332;
	}

	std::cout << "vars_lowerbound[i] = -1.0;" << std::endl;
	
	// Acceleration/decceleration upper and lower limits.
	// NOTE: Feel free to change this to something else.
	for (int i = a_start; i < n_vars; i++) {
		vars_lowerbound[i] = -1.0;
		vars_upperbound[i] = 1.0;
	}
	
	std::cout << "Lower and upper limits for constraints" << std::endl;

	// Lower and upper limits for constraints
	// All of these should be 0 except the initial
	// state indices.
	Dvector constraints_lowerbound(n_constraints);
	Dvector constraints_upperbound(n_constraints);
	for (int i = 0; i < n_constraints; i++) {
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

	std::cout << "solve the problem" << std::endl;
	
	// solve the problem
	CppAD::ipopt::solve<Dvector, FG_eval>(
	  options, vars, vars_lowerbound, vars_upperbound, constraints_lowerbound,
	  constraints_upperbound, fg_eval, solution);

	// Check some of the solution values
	ok &= solution.status == CppAD::ipopt::solve_result<Dvector>::success;

	// Cost
	auto cost = solution.obj_value;
	std::cout << "Cost " << cost << std::endl;
	
	std::cout << "return vector " << std::endl;

	// TODO: Return the first actuator values. The variables can be accessed with
	// `solution.x[i]`.
	//
	// {...} is shorthand for creating a vector, so auto x1 = {1.0,2.0}
	// creates a 2 element double vector.
	
	vector<double> result;

	result.push_back(solution.x[delta_start]);
	result.push_back(solution.x[a_start]);

	for (int i = 0; i < N-1; i++) {
	result.push_back(solution.x[x_start + i + 1]);
	result.push_back(solution.x[y_start + i + 1]);
	}

	return result;
	
	/*return {solution.x[delta_start], solution.x[a_start], solution.x[x_start], 
	solution.x[x_start+1] , solution.x[x_start+2] , solution.x[x_start+3] , solution.x[x_start+4], solution.x[y_start], 
	solution.x[y_start+1] , solution.x[y_start+2] , solution.x[y_start+3] , solution.x[y_start+4]};*/
}
