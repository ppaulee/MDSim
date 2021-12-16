// Copyright (c) 2005-2014 Code Synthesis Tools CC
//
// This program was generated by CodeSynthesis XSD, an XML Schema to
// C++ data binding compiler.
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License version 2 as
// published by the Free Software Foundation.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA
//
// In addition, as a special exception, Code Synthesis Tools CC gives
// permission to link this program with the Xerces-C++ library (or with
// modified versions of Xerces-C++ that use the same license as Xerces-C++),
// and distribute linked combinations including the two. You must obey
// the GNU General Public License version 2 in all respects for all of
// the code used other than Xerces-C++. If you modify this copy of the
// program, you may extend this exception to your version of the program,
// but you are not obligated to do so. If you do not wish to do so, delete
// this exception statement from your version.
//
// Furthermore, Code Synthesis Tools CC makes a special exception for
// the Free/Libre and Open Source Software (FLOSS) which is described
// in the accompanying FLOSSE file.
//

// Begin prologue.
//
//
// End prologue.

#include <xsd/cxx/pre.hxx>

#include "MolSimSkl.h"

// dimension_pskel
//

void dimension_pskel::
x_parser (::xml_schema::int_pskel& p)
{
  this->x_parser_ = &p;
}

void dimension_pskel::
y_parser (::xml_schema::int_pskel& p)
{
  this->y_parser_ = &p;
}

void dimension_pskel::
z_parser (::xml_schema::int_pskel& p)
{
  this->z_parser_ = &p;
}

void dimension_pskel::
parsers (::xml_schema::int_pskel& x,
         ::xml_schema::int_pskel& y,
         ::xml_schema::int_pskel& z)
{
  this->x_parser_ = &x;
  this->y_parser_ = &y;
  this->z_parser_ = &z;
}

dimension_pskel::
dimension_pskel ()
: x_parser_ (0),
  y_parser_ (0),
  z_parser_ (0)
{
}

// point_pskel
//

void point_pskel::
x_parser (::xml_schema::double_pskel& p)
{
  this->x_parser_ = &p;
}

void point_pskel::
y_parser (::xml_schema::double_pskel& p)
{
  this->y_parser_ = &p;
}

void point_pskel::
z_parser (::xml_schema::double_pskel& p)
{
  this->z_parser_ = &p;
}

void point_pskel::
parsers (::xml_schema::double_pskel& x,
         ::xml_schema::double_pskel& y,
         ::xml_schema::double_pskel& z)
{
  this->x_parser_ = &x;
  this->y_parser_ = &y;
  this->z_parser_ = &z;
}

point_pskel::
point_pskel ()
: x_parser_ (0),
  y_parser_ (0),
  z_parser_ (0)
{
}

// velocity_pskel
//

void velocity_pskel::
x_parser (::xml_schema::double_pskel& p)
{
  this->x_parser_ = &p;
}

void velocity_pskel::
y_parser (::xml_schema::double_pskel& p)
{
  this->y_parser_ = &p;
}

void velocity_pskel::
z_parser (::xml_schema::double_pskel& p)
{
  this->z_parser_ = &p;
}

void velocity_pskel::
parsers (::xml_schema::double_pskel& x,
         ::xml_schema::double_pskel& y,
         ::xml_schema::double_pskel& z)
{
  this->x_parser_ = &x;
  this->y_parser_ = &y;
  this->z_parser_ = &z;
}

velocity_pskel::
velocity_pskel ()
: x_parser_ (0),
  y_parser_ (0),
  z_parser_ (0)
{
}

// simulationContainer_pskel
//

void simulationContainer_pskel::
boundaryConditions_parser (::boundaryConditions_pskel& p)
{
  this->boundaryConditions_parser_ = &p;
}

void simulationContainer_pskel::
dimension_parser (::dimension_pskel& p)
{
  this->dimension_parser_ = &p;
}

void simulationContainer_pskel::
mesh_parser (::xml_schema::double_pskel& p)
{
  this->mesh_parser_ = &p;
}

void simulationContainer_pskel::
cutOff_parser (::xml_schema::double_pskel& p)
{
  this->cutOff_parser_ = &p;
}

void simulationContainer_pskel::
containerAlgorithm_parser (::containerAlgorithm_pskel& p)
{
  this->containerAlgorithm_parser_ = &p;
}

void simulationContainer_pskel::
parsers (::boundaryConditions_pskel& boundaryConditions,
         ::dimension_pskel& dimension,
         ::xml_schema::double_pskel& mesh,
         ::xml_schema::double_pskel& cutOff,
         ::containerAlgorithm_pskel& containerAlgorithm)
{
  this->boundaryConditions_parser_ = &boundaryConditions;
  this->dimension_parser_ = &dimension;
  this->mesh_parser_ = &mesh;
  this->cutOff_parser_ = &cutOff;
  this->containerAlgorithm_parser_ = &containerAlgorithm;
}

simulationContainer_pskel::
simulationContainer_pskel ()
: boundaryConditions_parser_ (0),
  dimension_parser_ (0),
  mesh_parser_ (0),
  cutOff_parser_ (0),
  containerAlgorithm_parser_ (0)
{
}

// Cube_pskel
//

void Cube_pskel::
dimension_parser (::dimension_pskel& p)
{
  this->dimension_parser_ = &p;
}

void Cube_pskel::
startPoint_parser (::point_pskel& p)
{
  this->startPoint_parser_ = &p;
}

void Cube_pskel::
h_parser (::xml_schema::double_pskel& p)
{
  this->h_parser_ = &p;
}

void Cube_pskel::
mass_parser (::xml_schema::double_pskel& p)
{
  this->mass_parser_ = &p;
}

void Cube_pskel::
velocity_parser (::velocity_pskel& p)
{
  this->velocity_parser_ = &p;
}

void Cube_pskel::
epsilon_parser (::epsilon_pskel& p)
{
  this->epsilon_parser_ = &p;
}

void Cube_pskel::
sigma_parser (::sigma_pskel& p)
{
  this->sigma_parser_ = &p;
}

void Cube_pskel::
parsers (::dimension_pskel& dimension,
         ::point_pskel& startPoint,
         ::xml_schema::double_pskel& h,
         ::xml_schema::double_pskel& mass,
         ::velocity_pskel& velocity,
         ::epsilon_pskel& epsilon,
         ::sigma_pskel& sigma)
{
  this->dimension_parser_ = &dimension;
  this->startPoint_parser_ = &startPoint;
  this->h_parser_ = &h;
  this->mass_parser_ = &mass;
  this->velocity_parser_ = &velocity;
  this->epsilon_parser_ = &epsilon;
  this->sigma_parser_ = &sigma;
}

Cube_pskel::
Cube_pskel ()
: dimension_parser_ (0),
  startPoint_parser_ (0),
  h_parser_ (0),
  mass_parser_ (0),
  velocity_parser_ (0),
  epsilon_parser_ (0),
  sigma_parser_ (0)
{
}

// Sphere_pskel
//

void Sphere_pskel::
center_parser (::point_pskel& p)
{
  this->center_parser_ = &p;
}

void Sphere_pskel::
radius_parser (::xml_schema::double_pskel& p)
{
  this->radius_parser_ = &p;
}

void Sphere_pskel::
h_parser (::xml_schema::double_pskel& p)
{
  this->h_parser_ = &p;
}

void Sphere_pskel::
mass_parser (::xml_schema::double_pskel& p)
{
  this->mass_parser_ = &p;
}

void Sphere_pskel::
velocity_parser (::velocity_pskel& p)
{
  this->velocity_parser_ = &p;
}

void Sphere_pskel::
epsilon_parser (::epsilon_pskel& p)
{
  this->epsilon_parser_ = &p;
}

void Sphere_pskel::
sigma_parser (::sigma_pskel& p)
{
  this->sigma_parser_ = &p;
}

void Sphere_pskel::
parsers (::point_pskel& center,
         ::xml_schema::double_pskel& radius,
         ::xml_schema::double_pskel& h,
         ::xml_schema::double_pskel& mass,
         ::velocity_pskel& velocity,
         ::epsilon_pskel& epsilon,
         ::sigma_pskel& sigma)
{
  this->center_parser_ = &center;
  this->radius_parser_ = &radius;
  this->h_parser_ = &h;
  this->mass_parser_ = &mass;
  this->velocity_parser_ = &velocity;
  this->epsilon_parser_ = &epsilon;
  this->sigma_parser_ = &sigma;
}

Sphere_pskel::
Sphere_pskel ()
: center_parser_ (0),
  radius_parser_ (0),
  h_parser_ (0),
  mass_parser_ (0),
  velocity_parser_ (0),
  epsilon_parser_ (0),
  sigma_parser_ (0)
{
}

// particles_pskel
//

void particles_pskel::
Cube_parser (::Cube_pskel& p)
{
  this->Cube_parser_ = &p;
}

void particles_pskel::
Sphere_parser (::Sphere_pskel& p)
{
  this->Sphere_parser_ = &p;
}

void particles_pskel::
parsers (::Cube_pskel& Cube,
         ::Sphere_pskel& Sphere)
{
  this->Cube_parser_ = &Cube;
  this->Sphere_parser_ = &Sphere;
}

particles_pskel::
particles_pskel ()
: Cube_parser_ (0),
  Sphere_parser_ (0)
{
}

// thermostats_pskel
//

void thermostats_pskel::
initialTemperature_parser (::xml_schema::double_pskel& p)
{
  this->initialTemperature_parser_ = &p;
}

void thermostats_pskel::
targetTemperature_parser (::xml_schema::double_pskel& p)
{
  this->targetTemperature_parser_ = &p;
}

void thermostats_pskel::
maxDelta_parser (::xml_schema::double_pskel& p)
{
  this->maxDelta_parser_ = &p;
}

void thermostats_pskel::
stepSize_parser (::xml_schema::int_pskel& p)
{
  this->stepSize_parser_ = &p;
}

void thermostats_pskel::
parsers (::xml_schema::double_pskel& initialTemperature,
         ::xml_schema::double_pskel& targetTemperature,
         ::xml_schema::double_pskel& maxDelta,
         ::xml_schema::int_pskel& stepSize)
{
  this->initialTemperature_parser_ = &initialTemperature;
  this->targetTemperature_parser_ = &targetTemperature;
  this->maxDelta_parser_ = &maxDelta;
  this->stepSize_parser_ = &stepSize;
}

thermostats_pskel::
thermostats_pskel ()
: initialTemperature_parser_ (0),
  targetTemperature_parser_ (0),
  maxDelta_parser_ (0),
  stepSize_parser_ (0)
{
}

// molsim_pskel
//

void molsim_pskel::
input_file_parser (::xml_schema::string_pskel& p)
{
  this->input_file_parser_ = &p;
}

void molsim_pskel::
delta_t_parser (::time_pskel& p)
{
  this->delta_t_parser_ = &p;
}

void molsim_pskel::
end_time_parser (::time_pskel& p)
{
  this->end_time_parser_ = &p;
}

void molsim_pskel::
output_step_parser (::xml_schema::int_pskel& p)
{
  this->output_step_parser_ = &p;
}

void molsim_pskel::
epsilon_parser (::epsilon_pskel& p)
{
  this->epsilon_parser_ = &p;
}

void molsim_pskel::
sigma_parser (::sigma_pskel& p)
{
  this->sigma_parser_ = &p;
}

void molsim_pskel::
gravity_parser (::xml_schema::double_pskel& p)
{
  this->gravity_parser_ = &p;
}

void molsim_pskel::
averageV_parser (::xml_schema::double_pskel& p)
{
  this->averageV_parser_ = &p;
}

void molsim_pskel::
algorithm_parser (::algorithm_pskel& p)
{
  this->algorithm_parser_ = &p;
}

void molsim_pskel::
simulationContainer_parser (::simulationContainer_pskel& p)
{
  this->simulationContainer_parser_ = &p;
}

void molsim_pskel::
particles_parser (::particles_pskel& p)
{
  this->particles_parser_ = &p;
}

void molsim_pskel::
thermostats_parser (::thermostats_pskel& p)
{
  this->thermostats_parser_ = &p;
}

void molsim_pskel::
benchmark_parser (::benchmark_pskel& p)
{
  this->benchmark_parser_ = &p;
}

void molsim_pskel::
parsers (::xml_schema::string_pskel& input_file,
         ::time_pskel& delta_t,
         ::time_pskel& end_time,
         ::xml_schema::int_pskel& output_step,
         ::epsilon_pskel& epsilon,
         ::sigma_pskel& sigma,
         ::xml_schema::double_pskel& gravity,
         ::xml_schema::double_pskel& averageV,
         ::algorithm_pskel& algorithm,
         ::simulationContainer_pskel& simulationContainer,
         ::particles_pskel& particles,
         ::thermostats_pskel& thermostats,
         ::benchmark_pskel& benchmark)
{
  this->input_file_parser_ = &input_file;
  this->delta_t_parser_ = &delta_t;
  this->end_time_parser_ = &end_time;
  this->output_step_parser_ = &output_step;
  this->epsilon_parser_ = &epsilon;
  this->sigma_parser_ = &sigma;
  this->gravity_parser_ = &gravity;
  this->averageV_parser_ = &averageV;
  this->algorithm_parser_ = &algorithm;
  this->simulationContainer_parser_ = &simulationContainer;
  this->particles_parser_ = &particles;
  this->thermostats_parser_ = &thermostats;
  this->benchmark_parser_ = &benchmark;
}

molsim_pskel::
molsim_pskel ()
: input_file_parser_ (0),
  delta_t_parser_ (0),
  end_time_parser_ (0),
  output_step_parser_ (0),
  epsilon_parser_ (0),
  sigma_parser_ (0),
  gravity_parser_ (0),
  averageV_parser_ (0),
  algorithm_parser_ (0),
  simulationContainer_parser_ (0),
  particles_parser_ (0),
  thermostats_parser_ (0),
  benchmark_parser_ (0)
{
}

// time_pskel
//

// epsilon_pskel
//

// sigma_pskel
//

// dimension_pskel
//

void dimension_pskel::
x (int)
{
}

void dimension_pskel::
y (int)
{
}

void dimension_pskel::
z (int)
{
}

bool dimension_pskel::
_start_element_impl (const ::xml_schema::ro_string& ns,
                     const ::xml_schema::ro_string& n,
                     const ::xml_schema::ro_string* t)
{
  XSD_UNUSED (t);

  if (this->::xml_schema::complex_content::_start_element_impl (ns, n, t))
    return true;

  if (n == "x" && ns.empty ())
  {
    this->::xml_schema::complex_content::context_.top ().parser_ = this->x_parser_;

    if (this->x_parser_)
      this->x_parser_->pre ();

    return true;
  }

  if (n == "y" && ns.empty ())
  {
    this->::xml_schema::complex_content::context_.top ().parser_ = this->y_parser_;

    if (this->y_parser_)
      this->y_parser_->pre ();

    return true;
  }

  if (n == "z" && ns.empty ())
  {
    this->::xml_schema::complex_content::context_.top ().parser_ = this->z_parser_;

    if (this->z_parser_)
      this->z_parser_->pre ();

    return true;
  }

  return false;
}

bool dimension_pskel::
_end_element_impl (const ::xml_schema::ro_string& ns,
                   const ::xml_schema::ro_string& n)
{
  if (this->::xml_schema::complex_content::_end_element_impl (ns, n))
    return true;

  if (n == "x" && ns.empty ())
  {
    if (this->x_parser_)
      this->x (this->x_parser_->post_int ());

    return true;
  }

  if (n == "y" && ns.empty ())
  {
    if (this->y_parser_)
      this->y (this->y_parser_->post_int ());

    return true;
  }

  if (n == "z" && ns.empty ())
  {
    if (this->z_parser_)
      this->z (this->z_parser_->post_int ());

    return true;
  }

  return false;
}

// point_pskel
//

void point_pskel::
x (double)
{
}

void point_pskel::
y (double)
{
}

void point_pskel::
z (double)
{
}

bool point_pskel::
_start_element_impl (const ::xml_schema::ro_string& ns,
                     const ::xml_schema::ro_string& n,
                     const ::xml_schema::ro_string* t)
{
  XSD_UNUSED (t);

  if (this->::xml_schema::complex_content::_start_element_impl (ns, n, t))
    return true;

  if (n == "x" && ns.empty ())
  {
    this->::xml_schema::complex_content::context_.top ().parser_ = this->x_parser_;

    if (this->x_parser_)
      this->x_parser_->pre ();

    return true;
  }

  if (n == "y" && ns.empty ())
  {
    this->::xml_schema::complex_content::context_.top ().parser_ = this->y_parser_;

    if (this->y_parser_)
      this->y_parser_->pre ();

    return true;
  }

  if (n == "z" && ns.empty ())
  {
    this->::xml_schema::complex_content::context_.top ().parser_ = this->z_parser_;

    if (this->z_parser_)
      this->z_parser_->pre ();

    return true;
  }

  return false;
}

bool point_pskel::
_end_element_impl (const ::xml_schema::ro_string& ns,
                   const ::xml_schema::ro_string& n)
{
  if (this->::xml_schema::complex_content::_end_element_impl (ns, n))
    return true;

  if (n == "x" && ns.empty ())
  {
    if (this->x_parser_)
      this->x (this->x_parser_->post_double ());

    return true;
  }

  if (n == "y" && ns.empty ())
  {
    if (this->y_parser_)
      this->y (this->y_parser_->post_double ());

    return true;
  }

  if (n == "z" && ns.empty ())
  {
    if (this->z_parser_)
      this->z (this->z_parser_->post_double ());

    return true;
  }

  return false;
}

// velocity_pskel
//

void velocity_pskel::
x (double)
{
}

void velocity_pskel::
y (double)
{
}

void velocity_pskel::
z (double)
{
}


bool velocity_pskel::
_start_element_impl (const ::xml_schema::ro_string& ns,
                     const ::xml_schema::ro_string& n,
                     const ::xml_schema::ro_string* t)
{
  XSD_UNUSED (t);

  if (this->::xml_schema::complex_content::_start_element_impl (ns, n, t))
    return true;

  if (n == "x" && ns.empty ())
  {
    this->::xml_schema::complex_content::context_.top ().parser_ = this->x_parser_;

    if (this->x_parser_)
      this->x_parser_->pre ();

    return true;
  }

  if (n == "y" && ns.empty ())
  {
    this->::xml_schema::complex_content::context_.top ().parser_ = this->y_parser_;

    if (this->y_parser_)
      this->y_parser_->pre ();

    return true;
  }

  if (n == "z" && ns.empty ())
  {
    this->::xml_schema::complex_content::context_.top ().parser_ = this->z_parser_;

    if (this->z_parser_)
      this->z_parser_->pre ();

    return true;
  }

  return false;
}

bool velocity_pskel::
_end_element_impl (const ::xml_schema::ro_string& ns,
                   const ::xml_schema::ro_string& n)
{
  if (this->::xml_schema::complex_content::_end_element_impl (ns, n))
    return true;

  if (n == "x" && ns.empty ())
  {
    if (this->x_parser_)
      this->x (this->x_parser_->post_double ());

    return true;
  }

  if (n == "y" && ns.empty ())
  {
    if (this->y_parser_)
      this->y (this->y_parser_->post_double ());

    return true;
  }

  if (n == "z" && ns.empty ())
  {
    if (this->z_parser_)
      this->z (this->z_parser_->post_double ());

    return true;
  }

  return false;
}

// algorithm_pskel
//

// containerAlgorithm_pskel
//

// simulationContainer_pskel
//

void simulationContainer_pskel::
boundaryConditions (const ::library::boundaryConditions)
{
}

void simulationContainer_pskel::
dimension (const ::library::dimension)
{
}

void simulationContainer_pskel::
mesh (double)
{
}

void simulationContainer_pskel::
cutOff (double)
{
}

void simulationContainer_pskel::
containerAlgorithm (const ::library::containerAlgorithm)
{
}

bool simulationContainer_pskel::
_start_element_impl (const ::xml_schema::ro_string& ns,
                     const ::xml_schema::ro_string& n,
                     const ::xml_schema::ro_string* t)
{
  XSD_UNUSED (t);

  if (this->::xml_schema::complex_content::_start_element_impl (ns, n, t))
    return true;

  if (n == "boundaryConditions" && ns.empty ())
  {
    this->::xml_schema::complex_content::context_.top ().parser_ = this->boundaryConditions_parser_;

    if (this->boundaryConditions_parser_)
      this->boundaryConditions_parser_->pre ();

    return true;
  }

  if (n == "dimension" && ns.empty ())
  {
    this->::xml_schema::complex_content::context_.top ().parser_ = this->dimension_parser_;

    if (this->dimension_parser_)
      this->dimension_parser_->pre ();

    return true;
  }

  if (n == "mesh" && ns.empty ())
  {
    this->::xml_schema::complex_content::context_.top ().parser_ = this->mesh_parser_;

    if (this->mesh_parser_)
      this->mesh_parser_->pre ();

    return true;
  }

  if (n == "cutOff" && ns.empty ())
  {
    this->::xml_schema::complex_content::context_.top ().parser_ = this->cutOff_parser_;

    if (this->cutOff_parser_)
      this->cutOff_parser_->pre ();

    return true;
  }

  return false;
}

bool simulationContainer_pskel::
_end_element_impl (const ::xml_schema::ro_string& ns,
                   const ::xml_schema::ro_string& n)
{
  if (this->::xml_schema::complex_content::_end_element_impl (ns, n))
    return true;

  if (n == "boundaryConditions" && ns.empty ())
  {
    if (this->boundaryConditions_parser_)
    {
      this->boundaryConditions (this->boundaryConditions_parser_->post_boundaryConditions ());
    }

    return true;
  }

  if (n == "dimension" && ns.empty ())
  {
    if (this->dimension_parser_)
    {
      this->dimension (this->dimension_parser_->post_dimension ());
    }

    return true;
  }

  if (n == "mesh" && ns.empty ())
  {
    if (this->mesh_parser_)
      this->mesh (this->mesh_parser_->post_double ());

    return true;
  }

  if (n == "cutOff" && ns.empty ())
  {
    if (this->cutOff_parser_)
      this->cutOff (this->cutOff_parser_->post_double ());

    return true;
  }

  return false;
}

bool simulationContainer_pskel::
_attribute_impl (const ::xml_schema::ro_string& ns,
                 const ::xml_schema::ro_string& n,
                 const ::xml_schema::ro_string& v)
{
  if (this->::xml_schema::complex_content::_attribute_impl (ns, n, v))
    return true;

  if (n == "containerAlgorithm" && ns.empty ())
  {
    if (this->containerAlgorithm_parser_)
    {
      this->containerAlgorithm_parser_->pre ();
      this->containerAlgorithm_parser_->_pre_impl ();
      this->containerAlgorithm_parser_->_characters (v);
      this->containerAlgorithm_parser_->_post_impl ();
      this->containerAlgorithm (this->containerAlgorithm_parser_->post_containerAlgorithm ());
    }

    return true;
  }

  return false;
}

// boundaryConditions_pskel
//

// benchmark_pskel
//

// Cube_pskel
//

void Cube_pskel::
dimension (const ::library::dimension)
{
}

void Cube_pskel::
startPoint (const ::library::point)
{
}

void Cube_pskel::
h (double)
{
}

void Cube_pskel::
mass (double)
{
}

void Cube_pskel::
velocity (const ::library::velocity)
{
}

void Cube_pskel::
epsilon (const ::library::epsilon)
{
}

void Cube_pskel::
sigma (const ::library::sigma)
{
}

bool Cube_pskel::
_start_element_impl (const ::xml_schema::ro_string& ns,
                     const ::xml_schema::ro_string& n,
                     const ::xml_schema::ro_string* t)
{
  XSD_UNUSED (t);

  if (this->::xml_schema::complex_content::_start_element_impl (ns, n, t))
    return true;

  if (n == "dimension" && ns.empty ())
  {
    this->::xml_schema::complex_content::context_.top ().parser_ = this->dimension_parser_;

    if (this->dimension_parser_)
      this->dimension_parser_->pre ();

    return true;
  }

  if (n == "startPoint" && ns.empty ())
  {
    this->::xml_schema::complex_content::context_.top ().parser_ = this->startPoint_parser_;

    if (this->startPoint_parser_)
      this->startPoint_parser_->pre ();

    return true;
  }

  if (n == "h" && ns.empty ())
  {
    this->::xml_schema::complex_content::context_.top ().parser_ = this->h_parser_;

    if (this->h_parser_)
      this->h_parser_->pre ();

    return true;
  }

  if (n == "mass" && ns.empty ())
  {
    this->::xml_schema::complex_content::context_.top ().parser_ = this->mass_parser_;

    if (this->mass_parser_)
      this->mass_parser_->pre ();

    return true;
  }

  if (n == "velocity" && ns.empty ())
  {
    this->::xml_schema::complex_content::context_.top ().parser_ = this->velocity_parser_;

    if (this->velocity_parser_)
      this->velocity_parser_->pre ();

    return true;
  }

  if (n == "epsilon" && ns.empty ())
  {
    this->::xml_schema::complex_content::context_.top ().parser_ = this->epsilon_parser_;

    if (this->epsilon_parser_)
      this->epsilon_parser_->pre ();

    return true;
  }

  if (n == "sigma" && ns.empty ())
  {
    this->::xml_schema::complex_content::context_.top ().parser_ = this->sigma_parser_;

    if (this->sigma_parser_)
      this->sigma_parser_->pre ();

    return true;
  }

  return false;
}

bool Cube_pskel::
_end_element_impl (const ::xml_schema::ro_string& ns,
                   const ::xml_schema::ro_string& n)
{
  if (this->::xml_schema::complex_content::_end_element_impl (ns, n))
    return true;

  if (n == "dimension" && ns.empty ())
  {
    if (this->dimension_parser_)
    {
      this->dimension (this->dimension_parser_->post_dimension ());
    }

    return true;
  }

  if (n == "startPoint" && ns.empty ())
  {
    if (this->startPoint_parser_)
    {
      this->startPoint (this->startPoint_parser_->post_point ());
    }

    return true;
  }

  if (n == "h" && ns.empty ())
  {
    if (this->h_parser_)
      this->h (this->h_parser_->post_double ());

    return true;
  }

  if (n == "mass" && ns.empty ())
  {
    if (this->mass_parser_)
      this->mass (this->mass_parser_->post_double ());

    return true;
  }

  if (n == "velocity" && ns.empty ())
  {
    if (this->velocity_parser_)
    {
      this->velocity (this->velocity_parser_->post_velocity ());
    }

    return true;
  }

  if (n == "epsilon" && ns.empty ())
  {
    if (this->epsilon_parser_)
    {
      this->epsilon (this->epsilon_parser_->post_epsilon ());
    }

    return true;
  }

  if (n == "sigma" && ns.empty ())
  {
    if (this->sigma_parser_)
    {
      this->sigma (this->sigma_parser_->post_sigma ());
    }

    return true;
  }

  return false;
}

// Sphere_pskel
//

void Sphere_pskel::
center (const ::library::point)
{
}

void Sphere_pskel::
radius (double)
{
}

void Sphere_pskel::
h (double)
{
}

void Sphere_pskel::
mass (double)
{
}

void Sphere_pskel::
velocity (const ::library::velocity)
{
}

void Sphere_pskel::
epsilon (::library::epsilon)
{
}

void Sphere_pskel::
sigma (::library::epsilon)
{
}

bool Sphere_pskel::
_start_element_impl (const ::xml_schema::ro_string& ns,
                     const ::xml_schema::ro_string& n,
                     const ::xml_schema::ro_string* t)
{
  XSD_UNUSED (t);

  if (this->::xml_schema::complex_content::_start_element_impl (ns, n, t))
    return true;

  if (n == "center" && ns.empty ())
  {
    this->::xml_schema::complex_content::context_.top ().parser_ = this->center_parser_;

    if (this->center_parser_)
      this->center_parser_->pre ();

    return true;
  }

  if (n == "radius" && ns.empty ())
  {
    this->::xml_schema::complex_content::context_.top ().parser_ = this->radius_parser_;

    if (this->radius_parser_)
      this->radius_parser_->pre ();

    return true;
  }

  if (n == "h" && ns.empty ())
  {
    this->::xml_schema::complex_content::context_.top ().parser_ = this->h_parser_;

    if (this->h_parser_)
      this->h_parser_->pre ();

    return true;
  }

  if (n == "mass" && ns.empty ())
  {
    this->::xml_schema::complex_content::context_.top ().parser_ = this->mass_parser_;

    if (this->mass_parser_)
      this->mass_parser_->pre ();

    return true;
  }

  if (n == "velocity" && ns.empty ())
  {
    this->::xml_schema::complex_content::context_.top ().parser_ = this->velocity_parser_;

    if (this->velocity_parser_)
      this->velocity_parser_->pre ();

    return true;
  }

  if (n == "epsilon" && ns.empty ())
  {
    this->::xml_schema::complex_content::context_.top ().parser_ = this->epsilon_parser_;

    if (this->epsilon_parser_)
      this->epsilon_parser_->pre ();

    return true;
  }

  if (n == "sigma" && ns.empty ())
  {
    this->::xml_schema::complex_content::context_.top ().parser_ = this->sigma_parser_;

    if (this->sigma_parser_)
      this->sigma_parser_->pre ();

    return true;
  }

  return false;
}

bool Sphere_pskel::
_end_element_impl (const ::xml_schema::ro_string& ns,
                   const ::xml_schema::ro_string& n)
{
  if (this->::xml_schema::complex_content::_end_element_impl (ns, n))
    return true;

  if (n == "center" && ns.empty ())
  {
    if (this->center_parser_)
    {
      this->center (this->center_parser_->post_point ());
    }

    return true;
  }

  if (n == "radius" && ns.empty ())
  {
    if (this->radius_parser_)
      this->radius (this->radius_parser_->post_double ());

    return true;
  }

  if (n == "h" && ns.empty ())
  {
    if (this->h_parser_)
      this->h (this->h_parser_->post_double ());

    return true;
  }

  if (n == "mass" && ns.empty ())
  {
    if (this->mass_parser_)
      this->mass (this->mass_parser_->post_double ());

    return true;
  }

  if (n == "velocity" && ns.empty ())
  {
    if (this->velocity_parser_)
    {
      this->velocity (this->velocity_parser_->post_velocity ());
    }

    return true;
  }

  if (n == "epsilon" && ns.empty ())
  {
    if (this->epsilon_parser_)
    {
      this->epsilon (this->epsilon_parser_->post_epsilon ());
    }

    return true;
  }

  if (n == "sigma" && ns.empty ())
  {
    if (this->sigma_parser_)
    {
      this->sigma (this->sigma_parser_->post_sigma ());
    }

    return true;
  }

  return false;
}

// particles_pskel
//

void particles_pskel::
Cube (const ::library::Cube)
{
}

void particles_pskel::
Sphere (const ::library::Sphere)
{
}

bool particles_pskel::
_start_element_impl (const ::xml_schema::ro_string& ns,
                     const ::xml_schema::ro_string& n,
                     const ::xml_schema::ro_string* t)
{
  XSD_UNUSED (t);

  if (this->::xml_schema::complex_content::_start_element_impl (ns, n, t))
    return true;

  if (n == "Cube" && ns.empty ())
  {
    this->::xml_schema::complex_content::context_.top ().parser_ = this->Cube_parser_;

    if (this->Cube_parser_)
      this->Cube_parser_->pre ();

    return true;
  }

  if (n == "Sphere" && ns.empty ())
  {
    this->::xml_schema::complex_content::context_.top ().parser_ = this->Sphere_parser_;

    if (this->Sphere_parser_)
      this->Sphere_parser_->pre ();

    return true;
  }

  return false;
}

bool particles_pskel::
_end_element_impl (const ::xml_schema::ro_string& ns,
                   const ::xml_schema::ro_string& n)
{
  if (this->::xml_schema::complex_content::_end_element_impl (ns, n))
    return true;

  if (n == "Cube" && ns.empty ())
  {
    if (this->Cube_parser_)
    {
      this->Cube (this->Cube_parser_->post_Cube ());
    }

    return true;
  }

  if (n == "Sphere" && ns.empty ())
  {
    if (this->Sphere_parser_)
    {
      this->Sphere (this->Sphere_parser_->post_Sphere ());
    }

    return true;
  }

  return false;
}

// thermostats_pskel
//

void thermostats_pskel::
initialTemperature (double)
{
}

void thermostats_pskel::
targetTemperature (double)
{
}

void thermostats_pskel::
maxDelta (double)
{
}

void thermostats_pskel::
stepSize (int)
{
}

bool thermostats_pskel::
_start_element_impl (const ::xml_schema::ro_string& ns,
                     const ::xml_schema::ro_string& n,
                     const ::xml_schema::ro_string* t)
{
  XSD_UNUSED (t);

  if (this->::xml_schema::complex_content::_start_element_impl (ns, n, t))
    return true;

  if (n == "initialTemperature" && ns.empty ())
  {
    this->::xml_schema::complex_content::context_.top ().parser_ = this->initialTemperature_parser_;

    if (this->initialTemperature_parser_)
      this->initialTemperature_parser_->pre ();

    return true;
  }

  if (n == "targetTemperature" && ns.empty ())
  {
    this->::xml_schema::complex_content::context_.top ().parser_ = this->targetTemperature_parser_;

    if (this->targetTemperature_parser_)
      this->targetTemperature_parser_->pre ();

    return true;
  }

  if (n == "maxDelta" && ns.empty ())
  {
    this->::xml_schema::complex_content::context_.top ().parser_ = this->maxDelta_parser_;

    if (this->maxDelta_parser_)
      this->maxDelta_parser_->pre ();

    return true;
  }

  if (n == "stepSize" && ns.empty ())
  {
    this->::xml_schema::complex_content::context_.top ().parser_ = this->stepSize_parser_;

    if (this->stepSize_parser_)
      this->stepSize_parser_->pre ();

    return true;
  }

  return false;
}

bool thermostats_pskel::
_end_element_impl (const ::xml_schema::ro_string& ns,
                   const ::xml_schema::ro_string& n)
{
  if (this->::xml_schema::complex_content::_end_element_impl (ns, n))
    return true;

  if (n == "initialTemperature" && ns.empty ())
  {
    if (this->initialTemperature_parser_)
      this->initialTemperature (this->initialTemperature_parser_->post_double ());

    return true;
  }

  if (n == "targetTemperature" && ns.empty ())
  {
    if (this->targetTemperature_parser_)
      this->targetTemperature (this->targetTemperature_parser_->post_double ());

    return true;
  }

  if (n == "maxDelta" && ns.empty ())
  {
    if (this->maxDelta_parser_)
      this->maxDelta (this->maxDelta_parser_->post_double ());

    return true;
  }

  if (n == "stepSize" && ns.empty ())
  {
    if (this->stepSize_parser_)
      this->stepSize (this->stepSize_parser_->post_int ());

    return true;
  }

  return false;
}

// molsim_pskel
//

void molsim_pskel::
input_file (const ::std::string&)
{
}

void molsim_pskel::
delta_t (::library::time)
{
}

void molsim_pskel::
end_time (::library::time)
{
}

void molsim_pskel::
output_step (int)
{
}

void molsim_pskel::
epsilon (::library::epsilon)
{
}

void molsim_pskel::
sigma (::library::sigma)
{
}

void molsim_pskel::
gravity (double)
{
}

void molsim_pskel::
averageV (double)
{
}

void molsim_pskel::
algorithm (const std::string&)
{
}

void molsim_pskel::
simulationContainer (::library::simulationContainer)
{
}

void molsim_pskel::
particles (::library::particles)
{
}

void molsim_pskel::
thermostats (::library::thermostats)
{
}

void molsim_pskel::
benchmark (::library::benchmark)
{
}

bool molsim_pskel::
_start_element_impl (const ::xml_schema::ro_string& ns,
                     const ::xml_schema::ro_string& n,
                     const ::xml_schema::ro_string* t)
{
  XSD_UNUSED (t);

  if (this->::xml_schema::complex_content::_start_element_impl (ns, n, t))
    return true;

  if (n == "input_file" && ns.empty ())
  {
    this->::xml_schema::complex_content::context_.top ().parser_ = this->input_file_parser_;

    if (this->input_file_parser_)
      this->input_file_parser_->pre ();

    return true;
  }

  if (n == "delta_t" && ns.empty ())
  {
    this->::xml_schema::complex_content::context_.top ().parser_ = this->delta_t_parser_;

    if (this->delta_t_parser_)
      this->delta_t_parser_->pre ();

    return true;
  }

  if (n == "end_time" && ns.empty ())
  {
    this->::xml_schema::complex_content::context_.top ().parser_ = this->end_time_parser_;

    if (this->end_time_parser_)
      this->end_time_parser_->pre ();

    return true;
  }

  if (n == "output_step" && ns.empty ())
  {
    this->::xml_schema::complex_content::context_.top ().parser_ = this->output_step_parser_;

    if (this->output_step_parser_)
      this->output_step_parser_->pre ();

    return true;
  }

  if (n == "epsilon" && ns.empty ())
  {
    this->::xml_schema::complex_content::context_.top ().parser_ = this->epsilon_parser_;

    if (this->epsilon_parser_)
      this->epsilon_parser_->pre ();

    return true;
  }

  if (n == "sigma" && ns.empty ())
  {
    this->::xml_schema::complex_content::context_.top ().parser_ = this->sigma_parser_;

    if (this->sigma_parser_)
      this->sigma_parser_->pre ();

    return true;
  }

  if (n == "gravity" && ns.empty ())
  {
    this->::xml_schema::complex_content::context_.top ().parser_ = this->gravity_parser_;

    if (this->gravity_parser_)
      this->gravity_parser_->pre ();

    return true;
  }

  if (n == "averageV" && ns.empty ())
  {
    this->::xml_schema::complex_content::context_.top ().parser_ = this->averageV_parser_;

    if (this->averageV_parser_)
      this->averageV_parser_->pre ();

    return true;
  }

  if (n == "algorithm" && ns.empty ())
  {
    this->::xml_schema::complex_content::context_.top ().parser_ = this->algorithm_parser_;

    if (this->algorithm_parser_)
      this->algorithm_parser_->pre ();

    return true;
  }

  if (n == "simulationContainer" && ns.empty ())
  {
    this->::xml_schema::complex_content::context_.top ().parser_ = this->simulationContainer_parser_;

    if (this->simulationContainer_parser_)
      this->simulationContainer_parser_->pre ();

    return true;
  }

  if (n == "particles" && ns.empty ())
  {
    this->::xml_schema::complex_content::context_.top ().parser_ = this->particles_parser_;

    if (this->particles_parser_)
      this->particles_parser_->pre ();

    return true;
  }

  if (n == "thermostats" && ns.empty ())
  {
    this->::xml_schema::complex_content::context_.top ().parser_ = this->thermostats_parser_;

    if (this->thermostats_parser_)
      this->thermostats_parser_->pre ();

    return true;
  }

  if (n == "benchmark" && ns.empty ())
  {
    this->::xml_schema::complex_content::context_.top ().parser_ = this->benchmark_parser_;

    if (this->benchmark_parser_)
      this->benchmark_parser_->pre ();

    return true;
  }

  return false;
}

bool molsim_pskel::
_end_element_impl (const ::xml_schema::ro_string& ns,
                   const ::xml_schema::ro_string& n)
{
  if (this->::xml_schema::complex_content::_end_element_impl (ns, n))
    return true;

  if (n == "input_file" && ns.empty ())
  {
    if (this->input_file_parser_)
      this->input_file (this->input_file_parser_->post_string ());

    return true;
  }

  if (n == "delta_t" && ns.empty ())
  {
    if (this->delta_t_parser_)
    {
      this->delta_t (this->delta_t_parser_->post_time ());
    }

    return true;
  }

  if (n == "end_time" && ns.empty ())
  {
    if (this->end_time_parser_)
    {
      this->end_time (this->end_time_parser_->post_time ());
    }

    return true;
  }

  if (n == "output_step" && ns.empty ())
  {
    if (this->output_step_parser_)
      this->output_step (this->output_step_parser_->post_int ());

    return true;
  }

  if (n == "epsilon" && ns.empty ())
  {
    if (this->epsilon_parser_)
    {
      this->epsilon (this->epsilon_parser_->post_epsilon ());
    }

    return true;
  }

  if (n == "sigma" && ns.empty ())
  {
    if (this->sigma_parser_)
    {
      this->sigma (this->sigma_parser_->post_sigma ());
    }

    return true;
  }

  if (n == "gravity" && ns.empty ())
  {
    if (this->gravity_parser_)
      this->gravity (this->gravity_parser_->post_double ());

    return true;
  }

  if (n == "averageV" && ns.empty ())
  {
    if (this->averageV_parser_)
      this->averageV (this->averageV_parser_->post_double ());

    return true;
  }

  if (n == "algorithm" && ns.empty ())
  {
    if (this->algorithm_parser_)
    {
      this->algorithm (this->algorithm_parser_->post_algorithm ());
    }

    return true;
  }

  if (n == "simulationContainer" && ns.empty ())
  {
    if (this->simulationContainer_parser_)
    {
      this->simulationContainer (this->simulationContainer_parser_->post_simulationContainer ());
    }

    return true;
  }

  if (n == "particles" && ns.empty ())
  {
    if (this->particles_parser_)
    {
      this->particles (this->particles_parser_->post_particles ());
    }

    return true;
  }

  if (n == "thermostats" && ns.empty ())
  {
    if (this->thermostats_parser_)
    {
      this->thermostats ( this->thermostats_parser_->post_thermostats ());
    }

    return true;
  }

  if (n == "benchmark" && ns.empty ())
  {
    if (this->benchmark_parser_)
    {
      this->benchmark (this->benchmark_parser_->post_benchmark ());
    }

    return true;
  }

  return false;
}

#include <xsd/cxx/post.hxx>

// Begin epilogue.
//
//
// End epilogue.

