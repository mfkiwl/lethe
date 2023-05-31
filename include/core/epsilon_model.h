/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2019 -  by the Lethe authors
 *
 * This file is part of the Lethe library
 *
 * The Lethe library is free software; you can use it, redistribute
 * it, and/or modify it under the terms of the GNU Lesser General
 * Public License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * The full text of the license can be found in the file LICENSE at
 * the top level of the Lethe distribution.
 *
 * ---------------------------------------------------------------------
 */

#ifndef lethe_epsilon_model_h
#define lethe_epsilon_model_h

#include <core/physical_property_model.h>

/**
 * @brief EpsilonModel. Abstract class that allows to calculate epsilon.
 */
class EpsilonModel : public PhysicalPropertyModel
{
public:
  /**
   * @brief Instantiates and returns a pointer to a EpsilonModel object by casting it to
   * the proper child class
   *
   * @param material_properties Parameters for a single fluid
   */
  static std::shared_ptr<EpsilonModel>
  model_cast(const Parameters::Material &material_properties);
};


/**
 * @brief Constant epsilon, the only case considered
 */
class ConstantEpsilon : public EpsilonModel
{
public:
  /**
   * @brief Default constructor
   */
  ConstantEpsilon(const double epsilon_constant)
    : epsilon(epsilon_constant)
  {}

  /**
   * @brief value Calculates epsilon
   * @param fields_value Value of the various field on which the property may depend.
   * @return value of the physical property calculated with the fields_value
   */
  double
  value(const std::map<field, double> & /*fields_value*/) override
  {
    return epsilon;
  }

  /**
   * @brief vector_value Calculates the vector of epsilon
   * @param field_vectors Vectors of the fields on which epsilon may depend.
   * @param property_vector Vectors of epsilon values
   */
  void
  vector_value(const std::map<field, std::vector<double>> & /*field_vectors*/,
               std::vector<double> &property_vector) override
  {
    std::fill(property_vector.begin(), property_vector.end(), epsilon);
  }

  /**
   * @brief jacobian Calculates the jacobian (the partial derivative) of epsilon with respect to a field
   * @param field_values Value of the various fields on which the property may depend.
   * @param id Indicator of the field with respect to which the jacobian
   * should be calculated.
   * @return value of the partial derivative of epsilon with respect to the field.
   */

  double
  jacobian(const std::map<field, double> & /*field_values*/,
           field /*id*/) override
  {
    return 0;
  }

  /**
   * @brief vector_jacobian Calculates the derivative of epsilon with respect to a field.
   * @param field_vectors Vector for the values of the fields used to evaluate the property.
   * @param id Identifier of the field with respect to which a derivative should be calculated.
   * @param jacobian vector of the value of the derivative of epsilon with respect to the field id.
   */

  void
  vector_jacobian(
    const std::map<field, std::vector<double>> & /*field_vectors*/,
    const field /*id*/,
    std::vector<double> &jacobian_vector) override
  {
    std::fill(jacobian_vector.begin(), jacobian_vector.end(), 0);
  }

private:
  const double epsilon;
};

#endif
