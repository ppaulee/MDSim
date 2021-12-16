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

#ifndef MOL_SIM_SKL_H
#define MOL_SIM_SKL_H

#ifndef XSD_CXX11
#define XSD_CXX11
#endif

// Begin prologue.
//
//
// End prologue.

#include <xsd/cxx/config.hxx>

#if (XSD_INT_VERSION != 4000000L)
#error XSD runtime version mismatch
#endif

#include <xsd/cxx/pre.hxx>

// Forward declarations
//
class time_pskel;
class epsilon_pskel;
class sigma_pskel;
class dimension_pskel;
class point_pskel;
class velocity_pskel;
class algorithm_pskel;
class containerAlgorithm_pskel;
class simulationContainer_pskel;
class boundaryConditions_pskel;
class benchmark_pskel;
class Cube_pskel;
class Sphere_pskel;
class particles_pskel;
class thermostats_pskel;
class molsim_pskel;

#ifndef XSD_USE_CHAR
#define XSD_USE_CHAR
#endif

#ifndef XSD_CXX_PARSER_USE_CHAR
#define XSD_CXX_PARSER_USE_CHAR
#endif

#include <xsd/cxx/xml/char-utf8.hxx>
#include <xsd/cxx/xml/error-handler.hxx>
#include <xsd/cxx/parser/exceptions.hxx>
#include <xsd/cxx/parser/elements.hxx>
#include <xsd/cxx/parser/xml-schema.hxx>
#include <xsd/cxx/parser/non-validating/parser.hxx>
#include <xsd/cxx/parser/non-validating/xml-schema-pskel.hxx>
#include <xsd/cxx/parser/non-validating/xml-schema-pimpl.hxx>
#include <xsd/cxx/parser/xerces/elements.hxx>

#include "library.h"

namespace xml_schema
{
  // Built-in XML Schema types mapping.
  //
  typedef ::xsd::cxx::parser::string_sequence< char > string_sequence;
  typedef ::xsd::cxx::parser::qname< char > qname;
  typedef ::xsd::cxx::parser::buffer buffer;
  typedef ::xsd::cxx::parser::time_zone time_zone;
  typedef ::xsd::cxx::parser::gday gday;
  typedef ::xsd::cxx::parser::gmonth gmonth;
  typedef ::xsd::cxx::parser::gyear gyear;
  typedef ::xsd::cxx::parser::gmonth_day gmonth_day;
  typedef ::xsd::cxx::parser::gyear_month gyear_month;
  typedef ::xsd::cxx::parser::date date;
  typedef ::xsd::cxx::parser::time time;
  typedef ::xsd::cxx::parser::date_time date_time;
  typedef ::xsd::cxx::parser::duration duration;

  // Base parser skeletons.
  //
  typedef ::xsd::cxx::parser::parser_base< char > parser_base;
  typedef ::xsd::cxx::parser::non_validating::empty_content< char > empty_content;
  typedef ::xsd::cxx::parser::non_validating::simple_content< char > simple_content;
  typedef ::xsd::cxx::parser::non_validating::complex_content< char > complex_content;
  typedef ::xsd::cxx::parser::non_validating::list_base< char > list_base;

  // Parser skeletons and implementations for the XML Schema
  // built-in types.
  //
  typedef ::xsd::cxx::parser::non_validating::any_type_pskel< char > any_type_pskel;
  typedef ::xsd::cxx::parser::non_validating::any_type_pimpl< char > any_type_pimpl;

  typedef ::xsd::cxx::parser::non_validating::any_simple_type_pskel< char > any_simple_type_pskel;
  typedef ::xsd::cxx::parser::non_validating::any_simple_type_pimpl< char > any_simple_type_pimpl;

  typedef ::xsd::cxx::parser::non_validating::byte_pskel< char > byte_pskel;
  typedef ::xsd::cxx::parser::non_validating::byte_pimpl< char > byte_pimpl;

  typedef ::xsd::cxx::parser::non_validating::unsigned_byte_pskel< char > unsigned_byte_pskel;
  typedef ::xsd::cxx::parser::non_validating::unsigned_byte_pimpl< char > unsigned_byte_pimpl;

  typedef ::xsd::cxx::parser::non_validating::short_pskel< char > short_pskel;
  typedef ::xsd::cxx::parser::non_validating::short_pimpl< char > short_pimpl;

  typedef ::xsd::cxx::parser::non_validating::unsigned_short_pskel< char > unsigned_short_pskel;
  typedef ::xsd::cxx::parser::non_validating::unsigned_short_pimpl< char > unsigned_short_pimpl;

  typedef ::xsd::cxx::parser::non_validating::int_pskel< char > int_pskel;
  typedef ::xsd::cxx::parser::non_validating::int_pimpl< char > int_pimpl;

  typedef ::xsd::cxx::parser::non_validating::unsigned_int_pskel< char > unsigned_int_pskel;
  typedef ::xsd::cxx::parser::non_validating::unsigned_int_pimpl< char > unsigned_int_pimpl;

  typedef ::xsd::cxx::parser::non_validating::long_pskel< char > long_pskel;
  typedef ::xsd::cxx::parser::non_validating::long_pimpl< char > long_pimpl;

  typedef ::xsd::cxx::parser::non_validating::unsigned_long_pskel< char > unsigned_long_pskel;
  typedef ::xsd::cxx::parser::non_validating::unsigned_long_pimpl< char > unsigned_long_pimpl;

  typedef ::xsd::cxx::parser::non_validating::integer_pskel< char > integer_pskel;
  typedef ::xsd::cxx::parser::non_validating::integer_pimpl< char > integer_pimpl;

  typedef ::xsd::cxx::parser::non_validating::non_positive_integer_pskel< char > non_positive_integer_pskel;
  typedef ::xsd::cxx::parser::non_validating::non_positive_integer_pimpl< char > non_positive_integer_pimpl;

  typedef ::xsd::cxx::parser::non_validating::non_negative_integer_pskel< char > non_negative_integer_pskel;
  typedef ::xsd::cxx::parser::non_validating::non_negative_integer_pimpl< char > non_negative_integer_pimpl;

  typedef ::xsd::cxx::parser::non_validating::positive_integer_pskel< char > positive_integer_pskel;
  typedef ::xsd::cxx::parser::non_validating::positive_integer_pimpl< char > positive_integer_pimpl;

  typedef ::xsd::cxx::parser::non_validating::negative_integer_pskel< char > negative_integer_pskel;
  typedef ::xsd::cxx::parser::non_validating::negative_integer_pimpl< char > negative_integer_pimpl;

  typedef ::xsd::cxx::parser::non_validating::boolean_pskel< char > boolean_pskel;
  typedef ::xsd::cxx::parser::non_validating::boolean_pimpl< char > boolean_pimpl;

  typedef ::xsd::cxx::parser::non_validating::float_pskel< char > float_pskel;
  typedef ::xsd::cxx::parser::non_validating::float_pimpl< char > float_pimpl;

  typedef ::xsd::cxx::parser::non_validating::double_pskel< char > double_pskel;
  typedef ::xsd::cxx::parser::non_validating::double_pimpl< char > double_pimpl;

  typedef ::xsd::cxx::parser::non_validating::decimal_pskel< char > decimal_pskel;
  typedef ::xsd::cxx::parser::non_validating::decimal_pimpl< char > decimal_pimpl;

  typedef ::xsd::cxx::parser::non_validating::string_pskel< char > string_pskel;
  typedef ::xsd::cxx::parser::non_validating::string_pimpl< char > string_pimpl;

  typedef ::xsd::cxx::parser::non_validating::normalized_string_pskel< char > normalized_string_pskel;
  typedef ::xsd::cxx::parser::non_validating::normalized_string_pimpl< char > normalized_string_pimpl;

  typedef ::xsd::cxx::parser::non_validating::token_pskel< char > token_pskel;
  typedef ::xsd::cxx::parser::non_validating::token_pimpl< char > token_pimpl;

  typedef ::xsd::cxx::parser::non_validating::name_pskel< char > name_pskel;
  typedef ::xsd::cxx::parser::non_validating::name_pimpl< char > name_pimpl;

  typedef ::xsd::cxx::parser::non_validating::nmtoken_pskel< char > nmtoken_pskel;
  typedef ::xsd::cxx::parser::non_validating::nmtoken_pimpl< char > nmtoken_pimpl;

  typedef ::xsd::cxx::parser::non_validating::nmtokens_pskel< char > nmtokens_pskel;
  typedef ::xsd::cxx::parser::non_validating::nmtokens_pimpl< char > nmtokens_pimpl;

  typedef ::xsd::cxx::parser::non_validating::ncname_pskel< char > ncname_pskel;
  typedef ::xsd::cxx::parser::non_validating::ncname_pimpl< char > ncname_pimpl;

  typedef ::xsd::cxx::parser::non_validating::language_pskel< char > language_pskel;
  typedef ::xsd::cxx::parser::non_validating::language_pimpl< char > language_pimpl;

  typedef ::xsd::cxx::parser::non_validating::id_pskel< char > id_pskel;
  typedef ::xsd::cxx::parser::non_validating::id_pimpl< char > id_pimpl;

  typedef ::xsd::cxx::parser::non_validating::idref_pskel< char > idref_pskel;
  typedef ::xsd::cxx::parser::non_validating::idref_pimpl< char > idref_pimpl;

  typedef ::xsd::cxx::parser::non_validating::idrefs_pskel< char > idrefs_pskel;
  typedef ::xsd::cxx::parser::non_validating::idrefs_pimpl< char > idrefs_pimpl;

  typedef ::xsd::cxx::parser::non_validating::uri_pskel< char > uri_pskel;
  typedef ::xsd::cxx::parser::non_validating::uri_pimpl< char > uri_pimpl;

  typedef ::xsd::cxx::parser::non_validating::qname_pskel< char > qname_pskel;
  typedef ::xsd::cxx::parser::non_validating::qname_pimpl< char > qname_pimpl;

  typedef ::xsd::cxx::parser::non_validating::base64_binary_pskel< char > base64_binary_pskel;
  typedef ::xsd::cxx::parser::non_validating::base64_binary_pimpl< char > base64_binary_pimpl;

  typedef ::xsd::cxx::parser::non_validating::hex_binary_pskel< char > hex_binary_pskel;
  typedef ::xsd::cxx::parser::non_validating::hex_binary_pimpl< char > hex_binary_pimpl;

  typedef ::xsd::cxx::parser::non_validating::date_pskel< char > date_pskel;
  typedef ::xsd::cxx::parser::non_validating::date_pimpl< char > date_pimpl;

  typedef ::xsd::cxx::parser::non_validating::date_time_pskel< char > date_time_pskel;
  typedef ::xsd::cxx::parser::non_validating::date_time_pimpl< char > date_time_pimpl;

  typedef ::xsd::cxx::parser::non_validating::duration_pskel< char > duration_pskel;
  typedef ::xsd::cxx::parser::non_validating::duration_pimpl< char > duration_pimpl;

  typedef ::xsd::cxx::parser::non_validating::gday_pskel< char > gday_pskel;
  typedef ::xsd::cxx::parser::non_validating::gday_pimpl< char > gday_pimpl;

  typedef ::xsd::cxx::parser::non_validating::gmonth_pskel< char > gmonth_pskel;
  typedef ::xsd::cxx::parser::non_validating::gmonth_pimpl< char > gmonth_pimpl;

  typedef ::xsd::cxx::parser::non_validating::gmonth_day_pskel< char > gmonth_day_pskel;
  typedef ::xsd::cxx::parser::non_validating::gmonth_day_pimpl< char > gmonth_day_pimpl;

  typedef ::xsd::cxx::parser::non_validating::gyear_pskel< char > gyear_pskel;
  typedef ::xsd::cxx::parser::non_validating::gyear_pimpl< char > gyear_pimpl;

  typedef ::xsd::cxx::parser::non_validating::gyear_month_pskel< char > gyear_month_pskel;
  typedef ::xsd::cxx::parser::non_validating::gyear_month_pimpl< char > gyear_month_pimpl;

  typedef ::xsd::cxx::parser::non_validating::time_pskel< char > time_pskel;
  typedef ::xsd::cxx::parser::non_validating::time_pimpl< char > time_pimpl;

  // Exceptions. See libxsd/xsd/cxx/parser/exceptions.hxx for details.
  //
  typedef ::xsd::cxx::parser::exception< char > exception;

  // Parsing diagnostics.
  //
  typedef ::xsd::cxx::parser::severity severity;
  typedef ::xsd::cxx::parser::error< char > error;
  typedef ::xsd::cxx::parser::diagnostics< char > diagnostics;
  typedef ::xsd::cxx::parser::parsing< char > parsing;

  // Error handler. See libxsd/xsd/cxx/xml/error-handler.hxx for details.
  //
  typedef ::xsd::cxx::xml::error_handler< char > error_handler;

  // Read-only string.
  //
  typedef ::xsd::cxx::ro_string< char > ro_string;

  // Parsing flags. See libxsd/xsd/cxx/parser/xerces/elements.hxx
  // for details.
  //
  typedef ::xsd::cxx::parser::xerces::flags flags;

  // Parsing properties. See libxsd/xsd/cxx/parser/xerces/elements.hxx
  // for details.
  //
  typedef ::xsd::cxx::parser::xerces::properties< char > properties;

  // Document type. See libxsd/xsd/cxx/parser/xerces/elements.hxx
  // for details.
  //
  typedef ::xsd::cxx::parser::xerces::document< char > document;
}

class time_pskel: public virtual ::xml_schema::double_pskel
{
  public:
  // Parser callbacks. Override them in your implementation.
  //
  // virtual void
  // pre ();

  virtual library::time
  post_time ()=0;
};

class epsilon_pskel: public virtual ::xml_schema::double_pskel
{
  public:
  // Parser callbacks. Override them in your implementation.
  //
  // virtual void
  // pre ();

  virtual library::epsilon
  post_epsilon ()=0;
};

class sigma_pskel: public virtual ::xml_schema::double_pskel
{
  public:
  // Parser callbacks. Override them in your implementation.
  //
  // virtual void
  // pre ();

  virtual library::sigma
  post_sigma ()=0;
};

class dimension_pskel: public ::xml_schema::complex_content
{
  public:
  // Parser callbacks. Override them in your implementation.
  //
  // virtual void
  // pre ();

  virtual void
  x (int);

  virtual void
  y (int);

  virtual void
  z (int);

  virtual library::dimension
  post_dimension ()=0;

  // Parser construction API.
  //
  void
  x_parser (::xml_schema::int_pskel&);

  void
  y_parser (::xml_schema::int_pskel&);

  void
  z_parser (::xml_schema::int_pskel&);

  void
  parsers (::xml_schema::int_pskel& /* x */,
           ::xml_schema::int_pskel& /* y */,
           ::xml_schema::int_pskel& /* z */);

  // Constructor.
  //
  dimension_pskel ();

  // Implementation.
  //
  protected:
  virtual bool
  _start_element_impl (const ::xml_schema::ro_string&,
                       const ::xml_schema::ro_string&,
                       const ::xml_schema::ro_string*);

  virtual bool
  _end_element_impl (const ::xml_schema::ro_string&,
                     const ::xml_schema::ro_string&);

  protected:
  ::xml_schema::int_pskel* x_parser_;
  ::xml_schema::int_pskel* y_parser_;
  ::xml_schema::int_pskel* z_parser_;
};

class point_pskel: public ::xml_schema::complex_content
{
  public:
  // Parser callbacks. Override them in your implementation.
  //
  // virtual void
  // pre ();

  virtual void
  x (double);

  virtual void
  y (double);

  virtual void
  z (double);

  virtual library::point
  post_point ()=0;

  // Parser construction API.
  //
  void
  x_parser (::xml_schema::double_pskel&);

  void
  y_parser (::xml_schema::double_pskel&);

  void
  z_parser (::xml_schema::double_pskel&);

  void
  parsers (::xml_schema::double_pskel& /* x */,
           ::xml_schema::double_pskel& /* y */,
           ::xml_schema::double_pskel& /* z */);

  // Constructor.
  //
  point_pskel ();

  // Implementation.
  //
  protected:
  virtual bool
  _start_element_impl (const ::xml_schema::ro_string&,
                       const ::xml_schema::ro_string&,
                       const ::xml_schema::ro_string*);

  virtual bool
  _end_element_impl (const ::xml_schema::ro_string&,
                     const ::xml_schema::ro_string&);

  protected:
  ::xml_schema::double_pskel* x_parser_;
  ::xml_schema::double_pskel* y_parser_;
  ::xml_schema::double_pskel* z_parser_;
};

class velocity_pskel: public ::xml_schema::complex_content
{
  public:
  // Parser callbacks. Override them in your implementation.
  //
  // virtual void
  // pre ();

  virtual void
  x (double);

  virtual void
  y (double);

  virtual void
  z (double);

  virtual library::velocity
  post_velocity ()=0;

  // Parser construction API.
  //
  void
  x_parser (::xml_schema::double_pskel&);

  void
  y_parser (::xml_schema::double_pskel&);

  void
  z_parser (::xml_schema::double_pskel&);

  void
  parsers (::xml_schema::double_pskel& /* x */,
           ::xml_schema::double_pskel& /* y */,
           ::xml_schema::double_pskel& /* z */);

  // Constructor.
  //
  velocity_pskel ();

  // Implementation.
  //
  protected:
  virtual bool
  _start_element_impl (const ::xml_schema::ro_string&,
                       const ::xml_schema::ro_string&,
                       const ::xml_schema::ro_string*);

  virtual bool
  _end_element_impl (const ::xml_schema::ro_string&,
                     const ::xml_schema::ro_string&);

  protected:
  ::xml_schema::double_pskel* x_parser_;
  ::xml_schema::double_pskel* y_parser_;
  ::xml_schema::double_pskel* z_parser_;
};

class algorithm_pskel: public virtual ::xml_schema::string_pskel
{
  public:
  // Parser callbacks. Override them in your implementation.
  //
  // virtual void
  // pre ();

  virtual const std::string&
  post_algorithm ()=0;
};

class containerAlgorithm_pskel: public virtual ::xml_schema::string_pskel
{
  public:
  // Parser callbacks. Override them in your implementation.
  //
  // virtual void
  // pre ();

  virtual const std::string&
  post_containerAlgorithm ()=0;
};

class simulationContainer_pskel: public ::xml_schema::complex_content
{
  public:
  // Parser callbacks. Override them in your implementation.
  //
  // virtual void
  // pre ();

  virtual void
  boundaryConditions (const ::library::boundaryConditions);

  virtual void
  dimension (const ::library::dimension);

  virtual void
  mesh (double);

  virtual void
  cutOff (double);

  virtual void
  containerAlgorithm (const ::library::containerAlgorithm);

  virtual library::simulationContainer
  post_simulationContainer ()=0;

  // Parser construction API.
  //
  void
  boundaryConditions_parser (::boundaryConditions_pskel&);

  void
  dimension_parser (::dimension_pskel&);

  void
  mesh_parser (::xml_schema::double_pskel&);

  void
  cutOff_parser (::xml_schema::double_pskel&);

  void
  containerAlgorithm_parser (::containerAlgorithm_pskel&);

  void
  parsers (::boundaryConditions_pskel& /* boundaryConditions */,
           ::dimension_pskel& /* dimension */,
           ::xml_schema::double_pskel& /* mesh */,
           ::xml_schema::double_pskel& /* cutOff */,
           ::containerAlgorithm_pskel& /* containerAlgorithm */);

  // Constructor.
  //
  simulationContainer_pskel ();

  // Implementation.
  //
  protected:
  virtual bool
  _start_element_impl (const ::xml_schema::ro_string&,
                       const ::xml_schema::ro_string&,
                       const ::xml_schema::ro_string*);

  virtual bool
  _end_element_impl (const ::xml_schema::ro_string&,
                     const ::xml_schema::ro_string&);

  virtual bool
  _attribute_impl (const ::xml_schema::ro_string&,
                   const ::xml_schema::ro_string&,
                   const ::xml_schema::ro_string&);

  protected:
  ::boundaryConditions_pskel* boundaryConditions_parser_;
  ::dimension_pskel* dimension_parser_;
  ::xml_schema::double_pskel* mesh_parser_;
  ::xml_schema::double_pskel* cutOff_parser_;
  ::containerAlgorithm_pskel* containerAlgorithm_parser_;
};

class boundaryConditions_pskel: public virtual ::xml_schema::string_pskel
{
  public:
  // Parser callbacks. Override them in your implementation.
  //
  // virtual void
  // pre ();

  virtual library::boundaryConditions
  post_boundaryConditions ()=0;
};

class benchmark_pskel: public virtual ::xml_schema::string_pskel
{
  public:
  // Parser callbacks. Override them in your implementation.
  //
  // virtual void
  // pre ();

  virtual library::benchmark
  post_benchmark ()=0;
};

class Cube_pskel: public ::xml_schema::complex_content
{
  public:
  // Parser callbacks. Override them in your implementation.
  //
  // virtual void
  // pre ();

  virtual void
  dimension (const ::library::dimension);

  virtual void
  startPoint (const ::library::point);

  virtual void
  h (double);

  virtual void
  mass (double);

  virtual void
  velocity (const ::library::velocity);

  virtual void
  epsilon (const ::library::epsilon);

  virtual void
  sigma (const ::library::sigma);

  virtual library::Cube
  post_Cube ()=0;

  // Parser construction API.
  //
  void
  dimension_parser (::dimension_pskel&);

  void
  startPoint_parser (::point_pskel&);

  void
  h_parser (::xml_schema::double_pskel&);

  void
  mass_parser (::xml_schema::double_pskel&);

  void
  velocity_parser (::velocity_pskel&);

  void
  epsilon_parser (::epsilon_pskel&);

  void
  sigma_parser (::sigma_pskel&);

  void
  parsers (::dimension_pskel& /* dimension */,
           ::point_pskel& /* startPoint */,
           ::xml_schema::double_pskel& /* h */,
           ::xml_schema::double_pskel& /* mass */,
           ::velocity_pskel& /* velocity */,
           ::epsilon_pskel& /* epsilon */,
           ::sigma_pskel& /* sigma */);

  // Constructor.
  //
  Cube_pskel ();

  // Implementation.
  //
  protected:
  virtual bool
  _start_element_impl (const ::xml_schema::ro_string&,
                       const ::xml_schema::ro_string&,
                       const ::xml_schema::ro_string*);

  virtual bool
  _end_element_impl (const ::xml_schema::ro_string&,
                     const ::xml_schema::ro_string&);

  protected:
  ::dimension_pskel* dimension_parser_;
  ::point_pskel* startPoint_parser_;
  ::xml_schema::double_pskel* h_parser_;
  ::xml_schema::double_pskel* mass_parser_;
  ::velocity_pskel* velocity_parser_;
  ::epsilon_pskel* epsilon_parser_;
  ::sigma_pskel* sigma_parser_;
};

class Sphere_pskel: public ::xml_schema::complex_content
{
  public:
  // Parser callbacks. Override them in your implementation.
  //
  // virtual void
  // pre ();

  virtual void
  center (const ::library::point);

  virtual void
  radius (double);

  virtual void
  h (double);

  virtual void
  mass (double);

  virtual void
  velocity (const ::library::velocity);

  virtual void
  epsilon (::library::epsilon);

  virtual void
  sigma (::library::sigma);

  virtual library::Sphere
  post_Sphere ()=0;

  // Parser construction API.
  //
  void
  center_parser (::point_pskel&);

  void
  radius_parser (::xml_schema::double_pskel&);

  void
  h_parser (::xml_schema::double_pskel&);

  void
  mass_parser (::xml_schema::double_pskel&);

  void
  velocity_parser (::velocity_pskel&);

  void
  epsilon_parser (::epsilon_pskel&);

  void
  sigma_parser (::sigma_pskel&);

  void
  parsers (::point_pskel& /* center */,
           ::xml_schema::double_pskel& /* radius */,
           ::xml_schema::double_pskel& /* h */,
           ::xml_schema::double_pskel& /* mass */,
           ::velocity_pskel& /* velocity */,
           ::epsilon_pskel& /* epsilon */,
           ::sigma_pskel& /* sigma */);

  // Constructor.
  //
  Sphere_pskel ();

  // Implementation.
  //
  protected:
  virtual bool
  _start_element_impl (const ::xml_schema::ro_string&,
                       const ::xml_schema::ro_string&,
                       const ::xml_schema::ro_string*);

  virtual bool
  _end_element_impl (const ::xml_schema::ro_string&,
                     const ::xml_schema::ro_string&);

  protected:
  ::point_pskel* center_parser_;
  ::xml_schema::double_pskel* radius_parser_;
  ::xml_schema::double_pskel* h_parser_;
  ::xml_schema::double_pskel* mass_parser_;
  ::velocity_pskel* velocity_parser_;
  ::epsilon_pskel* epsilon_parser_;
  ::sigma_pskel* sigma_parser_;
};

class particles_pskel: public ::xml_schema::complex_content
{
  public:
  // Parser callbacks. Override them in your implementation.
  //
  // virtual void
  // pre ();

  virtual void
  Cube (const ::library::Cube);

  virtual void
  Sphere (::library::Sphere);

  virtual library::particles
  post_particles ()=0;

  // Parser construction API.
  //
  void
  Cube_parser (::Cube_pskel&);

  void
  Sphere_parser (::Sphere_pskel&);

  void
  parsers (::Cube_pskel& /* Cube */,
           ::Sphere_pskel& /* Sphere */);

  // Constructor.
  //
  particles_pskel ();

  // Implementation.
  //
  protected:
  virtual bool
  _start_element_impl (const ::xml_schema::ro_string&,
                       const ::xml_schema::ro_string&,
                       const ::xml_schema::ro_string*);

  virtual bool
  _end_element_impl (const ::xml_schema::ro_string&,
                     const ::xml_schema::ro_string&);

  protected:
  ::Cube_pskel* Cube_parser_;
  ::Sphere_pskel* Sphere_parser_;
};

class thermostats_pskel: public ::xml_schema::complex_content
{
  public:
  // Parser callbacks. Override them in your implementation.
  //
  // virtual void
  // pre ();

  virtual void
  initialTemperature (double);

  virtual void
  targetTemperature (double);

  virtual void
  maxDelta (double);

  virtual void
  stepSize (int);

  virtual library::thermostats
  post_thermostats ()=0;

  // Parser construction API.
  //
  void
  initialTemperature_parser (::xml_schema::double_pskel&);

  void
  targetTemperature_parser (::xml_schema::double_pskel&);

  void
  maxDelta_parser (::xml_schema::double_pskel&);

  void
  stepSize_parser (::xml_schema::int_pskel&);

  void
  parsers (::xml_schema::double_pskel& /* initialTemperature */,
           ::xml_schema::double_pskel& /* targetTemperature */,
           ::xml_schema::double_pskel& /* maxDelta */,
           ::xml_schema::int_pskel& /* stepSize */);

  // Constructor.
  //
  thermostats_pskel ();

  // Implementation.
  //
  protected:
  virtual bool
  _start_element_impl (const ::xml_schema::ro_string&,
                       const ::xml_schema::ro_string&,
                       const ::xml_schema::ro_string*);

  virtual bool
  _end_element_impl (const ::xml_schema::ro_string&,
                     const ::xml_schema::ro_string&);

  protected:
  ::xml_schema::double_pskel* initialTemperature_parser_;
  ::xml_schema::double_pskel* targetTemperature_parser_;
  ::xml_schema::double_pskel* maxDelta_parser_;
  ::xml_schema::int_pskel* stepSize_parser_;
};

class molsim_pskel: public ::xml_schema::complex_content
{
  public:
  // Parser callbacks. Override them in your implementation.
  //
  // virtual void
  // pre ();

  virtual void
  input_file (const ::std::string&);

  virtual void
  delta_t (::library::time);

  virtual void
  end_time (::library::time);

  virtual void
  output_step (int);

  virtual void
  epsilon (::library::epsilon);

  virtual void
  sigma (::library::epsilon);

  virtual void
  gravity (double);

  virtual void
  averageV (double);

  virtual void
  algorithm (const std::string&);

  virtual void
  simulationContainer (::library::simulationContainer);

  virtual void
  particles (::library::particles);

  virtual void
  thermostats (::library::thermostats);

  virtual void
  benchmark (::library::benchmark);

  virtual library::molsim
  post_molsim ()=0;

  // Parser construction API.
  //
  void
  input_file_parser (::xml_schema::string_pskel&);

  void
  delta_t_parser (::time_pskel&);

  void
  end_time_parser (::time_pskel&);

  void
  output_step_parser (::xml_schema::int_pskel&);

  void
  epsilon_parser (::epsilon_pskel&);

  void
  sigma_parser (::sigma_pskel&);

  void
  gravity_parser (::xml_schema::double_pskel&);

  void
  averageV_parser (::xml_schema::double_pskel&);

  void
  algorithm_parser (::algorithm_pskel&);

  void
  simulationContainer_parser (::simulationContainer_pskel&);

  void
  particles_parser (::particles_pskel&);

  void
  thermostats_parser (::thermostats_pskel&);

  void
  benchmark_parser (::benchmark_pskel&);

  void
  parsers (::xml_schema::string_pskel& /* input_file */,
           ::time_pskel& /* delta_t */,
           ::time_pskel& /* end_time */,
           ::xml_schema::int_pskel& /* output_step */,
           ::epsilon_pskel& /* epsilon */,
           ::sigma_pskel& /* sigma */,
           ::xml_schema::double_pskel& /* gravity */,
           ::xml_schema::double_pskel& /* averageV */,
           ::algorithm_pskel& /* algorithm */,
           ::simulationContainer_pskel& /* simulationContainer */,
           ::particles_pskel& /* particles */,
           ::thermostats_pskel& /* thermostats */,
           ::benchmark_pskel& /* benchmark */);

  // Constructor.
  //
  molsim_pskel ();

  // Implementation.
  //
  protected:
  virtual bool
  _start_element_impl (const ::xml_schema::ro_string&,
                       const ::xml_schema::ro_string&,
                       const ::xml_schema::ro_string*);

  virtual bool
  _end_element_impl (const ::xml_schema::ro_string&,
                     const ::xml_schema::ro_string&);

  protected:
  ::xml_schema::string_pskel* input_file_parser_;
  ::time_pskel* delta_t_parser_;
  ::time_pskel* end_time_parser_;
  ::xml_schema::int_pskel* output_step_parser_;
  ::epsilon_pskel* epsilon_parser_;
  ::sigma_pskel* sigma_parser_;
  ::xml_schema::double_pskel* gravity_parser_;
  ::xml_schema::double_pskel* averageV_parser_;
  ::algorithm_pskel* algorithm_parser_;
  ::simulationContainer_pskel* simulationContainer_parser_;
  ::particles_pskel* particles_parser_;
  ::thermostats_pskel* thermostats_parser_;
  ::benchmark_pskel* benchmark_parser_;
};

#include <xsd/cxx/post.hxx>

// Begin epilogue.
//
//
// End epilogue.

#endif // MOL_SIM_SKL_H