// Not copyrighted - public domain.
//
// This sample parser implementation was generated by CodeSynthesis XSD,
// an XML Schema to C++ data binding compiler. You may use it in your
// programs without any restrictions.
//

#include "MolSimImpl.h"
#include <iostream>

// benchmark_pimpl
//

void benchmark_pimpl::
pre ()
{
}

library::benchmark benchmark_pimpl::
post_benchmark ()
{
    const ::std::string& v (post_string ());

    // TODO
    return v;
    //
}

// time_pimpl
//

void time_pimpl::
pre ()
{
}

library::time time_pimpl::
post_time ()
{
    double v (post_double ());

    return v;
}

// epsilon_pimpl
//

void epsilon_pimpl::
pre ()
{
}

library::epsilon epsilon_pimpl::
post_epsilon ()
{
    double v (post_double ());

    return v;
}

// sigma_pimpl
//

void sigma_pimpl::
pre ()
{
}

library::sigma sigma_pimpl::
post_sigma ()
{
    double v (post_double ());

    return v;
}

// dimension_pimpl
//

void dimension_pimpl::
pre ()
{
}

void dimension_pimpl::
x (int x)
{
  // TODO
  dimension_.x(x);
  //
}

void dimension_pimpl::
y (int y)
{
  // TODO
  dimension_.y(y);
  //
}

void dimension_pimpl::
z (int z)
{
  // TODO
  dimension_.z(z);
  //
}

library::dimension dimension_pimpl::
post_dimension ()
{
    return dimension_;
}

// point_pimpl
//

void point_pimpl::
pre ()
{
}

void point_pimpl::
x (double x)
{
  // TODO
  point_.x(x);
  //
}

void point_pimpl::
y (double y)
{
  // TODO
  point_.y(y);
  //
}

void point_pimpl::
z (double z)
{
  // TODO
  point_.z(z);
  //
}

library::point point_pimpl::
post_point ()
{
    return point_;
}

// velocity_pimpl
//

void velocity_pimpl::
pre ()
{
}

void velocity_pimpl::
x (double x)
{
  // TODO
  velocity_.x(x);
  //
}

void velocity_pimpl::
y (double y)
{
  // TODO
  velocity_.y(y);
  //
}

void velocity_pimpl::
z (double z)
{
  // TODO
  velocity_.z(z);
  //
}

library::velocity velocity_pimpl::
post_velocity ()
{
    return velocity_;
}

// algorithm_pimpl
//

void algorithm_pimpl::
pre ()
{
}

const std::string& algorithm_pimpl::
post_algorithm ()
{
  // TODO
    const ::std::string& v (post_string ());
    algorithm_.assign(v);
    return algorithm_;
  //
}

// containerAlgorithm_pimpl
//

void containerAlgorithm_pimpl::
pre ()
{
}

library::containerAlgorithm containerAlgorithm_pimpl::
post_containerAlgorithm ()
{

  return post_string();
  //
}

// simulationContainer_pimpl
//

void simulationContainer_pimpl::
pre ()
{
}

void simulationContainer_pimpl::
boundaryConditions (library::boundaryConditions b)
{
    simulationContainer_.boundaryConditions(b);
}

void simulationContainer_pimpl::
dimension (library::dimension dimension)
{
    simulationContainer_.dimension(dimension);
}

void simulationContainer_pimpl::
mesh (double mesh)
{
  // TODO
  simulationContainer_.mesh(mesh);
  //
}

void simulationContainer_pimpl::
cutOff (double cutOff)
{
  // TODO
  simulationContainer_.cutOff(cutOff);
  //
}

void simulationContainer_pimpl::
containerAlgorithm (library::containerAlgorithm c)
{
    simulationContainer_.containerAlgorithm(c);
}

library::simulationContainer simulationContainer_pimpl::
post_simulationContainer ()
{
    return simulationContainer_;
}

// boundaryConditions_pimpl
//

void boundaryConditions_pimpl::
pre ()
{
}

library::boundaryConditions boundaryConditions_pimpl::
post_boundaryConditions ()
{
  // TODO
  return post_string();
  //
}



// Cube_pimpl
//

void Cube_pimpl::
pre ()
{
}

void Cube_pimpl::
dimension (library::dimension d)
{
    cube_.dimension(d);
}

void Cube_pimpl::
startPoint (library::point p)
{
    cube_.startPoint(p);
}

void Cube_pimpl::
h (double h)
{
  // TODO
  cube_.h(h);
  //
}

void Cube_pimpl::
mass (double mass)
{
  // TODO
  cube_.m(mass);
  //
}

void Cube_pimpl::
velocity (library::velocity v)
{
    cube_.velocity(v);
}

void Cube_pimpl::
epsilon (library::epsilon e)
{
    cube_.epsilon(e);
}

void Cube_pimpl::
sigma (library::sigma s)
{
    cube_.sigma(s);
}

library::Cube Cube_pimpl::
post_Cube ()
{
    return cube_;
}

// Sphere_pimpl
//

void Sphere_pimpl::
pre ()
{
}

void Sphere_pimpl::
center (library::point p)
{
    sphere_.center(p);
}

void Sphere_pimpl::
radius (double radius)
{
  // TODO
    sphere_.radius(radius);
  //
}

void Sphere_pimpl::
h (double h)
{
  // TODO
  sphere_.h(h);
  //
}

void Sphere_pimpl::
mass (double mass)
{
  // TODO
  sphere_.m(mass);
  //
}

void Sphere_pimpl::
velocity (library::velocity v)
{
    sphere_.velocity(v);
}

void Sphere_pimpl::
epsilon (library::epsilon e)
{
    sphere_.epsilon(e);
}

void Sphere_pimpl::
sigma (library::sigma s)
{
    sphere_.sigma(s);
}

library::Sphere Sphere_pimpl::
post_Sphere ()
{
    return sphere_;
}

// particles_pimpl
//

void particles_pimpl::
pre ()
{
}

void particles_pimpl::
Cube (library::Cube c)
{
    particles_.cube().push_back(c);
}

void particles_pimpl::
Sphere (library::Sphere s)
{
    particles_.sphere().push_back(s);
}

library::particles particles_pimpl::
post_particles ()
{
    return particles_;
}

// thermostats_pimpl
//

void thermostats_pimpl::
pre ()
{
}

void thermostats_pimpl::
initialTemperature (double initialTemperature)
{
  // TODO
  thermostats_.initialTemperature(initialTemperature);
  //
}

void thermostats_pimpl::
targetTemperature (double targetTemperature)
{
  // TODO
  thermostats_.targetTemperature(targetTemperature);
    std::cout << "ini di impl thermo, tt = " << thermostats_.targetTemperature() << std::endl;
  //
}

void thermostats_pimpl::
maxDelta (double maxDelta)
{
  // TODO
  thermostats_.maxDelta(maxDelta);
  //
}

void thermostats_pimpl::
stepSize (int stepSize)
{
  // TODO
  thermostats_.stepSize(stepSize);
  //
}

library::thermostats thermostats_pimpl::
post_thermostats ()
{
    return thermostats_;
}

// molsim_pimpl
//

void molsim_pimpl::
pre ()
{
}

void molsim_pimpl::
input_file (const ::std::string& input_file)
{
  // TODO
  molsim_.input(input_file);
  //
}

void molsim_pimpl::
delta_t (library::time t)
{
    molsim_.delta_t(t);
}

void molsim_pimpl::
end_time (library::time t)
{
    molsim_.endtime(t);
}

void molsim_pimpl::
output_step (int output_step)
{
  // TODO
  molsim_.outputStep(output_step);
  //
}

void molsim_pimpl::
epsilon (library::epsilon e)
{
    molsim_.epsilon(e);
}

void molsim_pimpl::
sigma (library::sigma s)
{
    molsim_.sigma(s);
}

void molsim_pimpl::
gravity (double gravity)
{
  // TODO
  molsim_.gravity(gravity);
  //
}

void molsim_pimpl::
averageV (double averageV)
{
  // TODO
  molsim_.averageV(averageV);
  //
}

void molsim_pimpl::
algorithm (const std::string& a)
{
    molsim_.algorithm(a);
}

void molsim_pimpl::
simulationContainer (library::simulationContainer s)
{
    molsim_.simu(s);
}

void molsim_pimpl::
particles (library::particles p)
{
    molsim_.particles(p);
}

void molsim_pimpl::
thermostats (library::thermostats t)
{
    molsim_.thermostats(t);
}

void molsim_pimpl::
benchmark (library::benchmark b)
{
    molsim_.benchmark(b);
}

library::molsim molsim_pimpl::
post_molsim ()
{
    return molsim_;
}

