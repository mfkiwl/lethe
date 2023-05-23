#include <core/bdf.h>
#include <core/time_integration_utilities.h>

#include <solvers/copy_data.h>
#include <solvers/cahn_hilliard_assemblers.h>


template <int dim>
void
CahnHilliardAssemblerCore<dim>::assemble_matrix(CahnHilliardScratchData<dim> &scratch_data,
                                          StabilizedMethodsCopyData &copy_data)
{
  // Scheme and physical properties
//  const std::vector<double> &diffusivity_vector =
//    scratch_data.tracer_diffusivity;
  const auto method = this->simulation_control->get_assembly_method();

  // Loop and quadrature informations
  const auto &       JxW_vec    = scratch_data.JxW;
  const unsigned int n_q_points = scratch_data.n_q_points;
  const unsigned int n_dofs     = scratch_data.n_dofs;
  const double       h          = scratch_data.cell_size;

  // Time steps and inverse time steps which is used for stabilization constant
  std::vector<double> time_steps_vector =
    this->simulation_control->get_time_steps_vector();
  const double dt  = time_steps_vector[0];
  const double sdt = 1. / dt;


  // Copy data elements
  auto &strong_jacobian_vec = copy_data.strong_jacobian;
  auto &local_matrix        = copy_data.local_matrix;

  // assembling local matrix and right hand side
  for (unsigned int q = 0; q < n_q_points; ++q)
    {
      // Gather into local variables the relevant fields
      //const double         diffusivity     = diffusivity_vector[q];
      const Tensor<1, dim> phase_order_gradient = scratch_data.phase_order_gradients[q];
      //const Tensor<1, dim> velocity        = scratch_data.velocity_values[q];
      const Tensor<1, dim> velocity;
      // Store JxW in local variable for faster access;
      const double JxW = JxW_vec[q];

      //Diffusivity parameter for the GLS stabilization
      //NEEDS TO BE CHANGED
      const double D = 1.0;

      // Shock capturing viscosity term
      const double order = scratch_data.fe_values_ch.get_fe().degree;

      // Calculation of the magnitude of the velocity for the
      // stabilization parameter
      const double u_mag = std::max(velocity.norm(), 1e-12);

      // Calculation of the GLS stabilization parameter. The
      // stabilization parameter used is different if the simulation is
      // steady or unsteady. In the unsteady case it includes the value
      // of the time-step
      const double tau =
        is_steady(method) ?
          1. / std::sqrt(std::pow(2. * u_mag / h, 2) +
                         9 * std::pow(4 * D / (h * h), 2)) :
          1. / std::sqrt(std::pow(sdt, 2) + std::pow(2. * u_mag / h, 2) +
                         9 * std::pow(4 * D/ (h * h), 2));

      for (unsigned int j = 0; j < n_dofs; ++j)
        {
          strong_jacobian_vec[q][j] += 0;

        }

      for (unsigned int i = 0; i < n_dofs; ++i)
        {


          for (unsigned int j = 0; j < n_dofs; ++j)
            {


              // Weak form : - D * laplacian T +  u * gradT - f=0
              local_matrix(i, j) += 0;



              local_matrix(i, j) += 0;
            }
        }
    }
} // end loop on quadrature points



template <int dim>
void
CahnHilliardAssemblerCore<dim>::assemble_rhs(CahnHilliardScratchData<dim> &   scratch_data,
                                       StabilizedMethodsCopyData &copy_data)
{
  // Scheme and physical properties
//  const std::vector<double> &diffusivity_vector =
//    scratch_data.tracer_diffusivity;
  const auto method = this->simulation_control->get_assembly_method();

  // Loop and quadrature informations
  const auto &       JxW_vec    = scratch_data.JxW;
  const unsigned int n_q_points = scratch_data.n_q_points;
  const unsigned int n_dofs     = scratch_data.n_dofs;
  const double       h          = scratch_data.cell_size;

  // Time steps and inverse time steps which is used for stabilization constant
  std::vector<double> time_steps_vector =
    this->simulation_control->get_time_steps_vector();
  const double dt  = time_steps_vector[0];
  const double sdt = 1. / dt;

  // Copy data elements
  auto &strong_residual_vec = copy_data.strong_residual;
  auto &local_rhs           = copy_data.local_rhs;

  // assembling local matrix and right hand side
  for (unsigned int q = 0; q < n_q_points; ++q)
    {
      // Gather into local variables the relevant fields
      const Tensor<1, dim> velocity;

      // Store JxW in local variable for faster access;
      const double JxW = JxW_vec[q];


      // Calculation of the magnitude of the velocity for the
      // stabilization parameter
      const double u_mag = std::max(velocity.norm(), 1e-12);
      const double D = 1.0;
      // Calculation of the GLS stabilization parameter. The
      // stabilization parameter used is different if the simulation is
      // steady or unsteady. In the unsteady case it includes the value
      // of the time-step
      const double tau =
        is_steady(method) ?
          1. / std::sqrt(std::pow(2. * u_mag / h, 2) +
                         9 * std::pow(4 * D / (h * h), 2)) :
          1. / std::sqrt(std::pow(sdt, 2) + std::pow(2. * u_mag / h, 2) +
                         9 * std::pow(4 * D / (h * h), 2));

      // Calculate the strong residual for GLS stabilization

      for (unsigned int i = 0; i < n_dofs; ++i)
        {

          // rhs for : - D * laplacian T +  u * grad T - f=0
          local_rhs(i) -= 0;

          local_rhs(i) -=
            0;

        }
    } // end loop on quadrature points
}

template class CahnHilliardAssemblerCore<2>;
template class CahnHilliardAssemblerCore<3>;

template <int dim>
void
CahnHilliardAssemblerBDF<dim>::assemble_matrix(CahnHilliardScratchData<dim> &scratch_data,
                                         StabilizedMethodsCopyData &copy_data)
{
  // Loop and quadrature informations
  const auto &       JxW        = scratch_data.JxW;
  const unsigned int n_q_points = scratch_data.n_q_points;
  const unsigned int n_dofs     = scratch_data.n_dofs;

  // Copy data elements
  auto &strong_residual = copy_data.strong_residual;
  auto &strong_jacobian = copy_data.strong_jacobian;
  auto &local_matrix    = copy_data.local_matrix;

  // Time stepping information
  const auto          method = this->simulation_control->get_assembly_method();
  std::vector<double> time_steps_vector =
    this->simulation_control->get_time_steps_vector();

  // Vector for the BDF coefficients
  Vector<double>      bdf_coefs = bdf_coefficients(method, time_steps_vector);
  std::vector<double> phase_order(1 + number_of_previous_solutions(method));

  // Loop over the quadrature points
  for (unsigned int q = 0; q < n_q_points; ++q)
    {
      phase_order[0] = scratch_data.phase_order_values[q];
      for (unsigned int p = 0; p < number_of_previous_solutions(method); ++p)
        phase_order[p + 1] = scratch_data.previous_phase_order_values[p][q];

      for (unsigned int p = 0; p < number_of_previous_solutions(method) + 1;
           ++p)
        {
          strong_residual[q] += bdf_coefs[p] * phase_order[p];
        }

      for (unsigned int j = 0; j < n_dofs; ++j)
        {
          strong_jacobian[q][j] += bdf_coefs[0] * scratch_data.phi_phase[q][j];
        }


      for (unsigned int i = 0; i < n_dofs; ++i)
        {
          const double phi_u_i = scratch_data.phi_phase[q][i];
          for (unsigned int j = 0; j < n_dofs; ++j)
            {
              const double phi_u_j = scratch_data.phi_phase[q][j];

              local_matrix(i, j) += phi_u_j * phi_u_i * bdf_coefs[0] * JxW[q];
            }
        }
    }
}

template <int dim>
void
CahnHilliardAssemblerBDF<dim>::assemble_rhs(CahnHilliardScratchData<dim> &   scratch_data,
                                      StabilizedMethodsCopyData &copy_data)
{
  // Loop and quadrature informations
  const auto &       JxW        = scratch_data.JxW;
  const unsigned int n_q_points = scratch_data.n_q_points;
  const unsigned int n_dofs     = scratch_data.n_dofs;

  // Copy data elements
  auto &strong_residual = copy_data.strong_residual;
  auto &local_rhs       = copy_data.local_rhs;

  // Time stepping information
  const auto          method = this->simulation_control->get_assembly_method();
  std::vector<double> time_steps_vector =
    this->simulation_control->get_time_steps_vector();

  // Vector for the BDF coefficients
  Vector<double>      bdf_coefs = bdf_coefficients(method, time_steps_vector);
  std::vector<double> phase_order(1 + number_of_previous_solutions(method));

  // Loop over the quadrature points
  for (unsigned int q = 0; q < n_q_points; ++q)
    {
      phase_order[0] = scratch_data.phase_order_values[q];
      for (unsigned int p = 0; p < number_of_previous_solutions(method); ++p)
        phase_order[p + 1] = scratch_data.previous_phase_order_values[p][q];

      for (unsigned int p = 0; p < number_of_previous_solutions(method) + 1;
           ++p)
        {
          strong_residual[q] += (bdf_coefs[p] * phase_order[p]);
        }

      for (unsigned int i = 0; i < n_dofs; ++i)
        {
          const double phi_u_i     = scratch_data.phi_phase[q][i];
          double       local_rhs_i = 0;
          for (unsigned int p = 0; p < number_of_previous_solutions(method) + 1;
               ++p)
            {
              local_rhs_i -= bdf_coefs[p] * (phase_order[p] * phi_u_i);
            }
          local_rhs(i) += local_rhs_i * JxW[q];
        }
    }
}

template class CahnHilliardAssemblerBDF<2>;
template class CahnHilliardAssemblerBDF<3>;
