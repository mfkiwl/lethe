#include <core/bdf.h>
#include <core/sdirk.h>

#include <solvers/cahn_hilliard_scratch_data.h>
template <int dim>
void
CahnHilliardScratchData<dim>::allocate()
{
  // Initialize size of arrays
  this->n_q_points = fe_values_ch.get_quadrature().size();
  this->n_dofs     = fe_values_ch.get_fe().n_dofs_per_cell();

  // Initialize arrays related to quadrature
  this->JxW = std::vector<double>(n_q_points);

  // Forcing term array
  this->source = std::vector<double>(n_q_points);

  // Initialize arrays related to velocity and pressure
  this->velocities.first_vector_component = 0;
  // Velocity
  this->velocity_values = std::vector<Tensor<1, dim>>(n_q_points);
  // Tracer
  this->phase_order_values        = std::vector<double>(n_q_points);
  this->phase_order_gradients     = std::vector<Tensor<1, dim>>(n_q_points);
  this->phase_order_laplacians    = std::vector<double>(n_q_points);
  this->chemical_potential_values        = std::vector<double>(n_q_points);
  this->chemical_potential_gradients     = std::vector<Tensor<1, dim>>(n_q_points);
  this->chemical_potential_laplacians    = std::vector<double>(n_q_points);

//  this->tracer_diffusivity   = std::vector<double>(n_q_points);
//  this->tracer_diffusivity_0 = std::vector<double>(n_q_points);
//  this->tracer_diffusivity_1 = std::vector<double>(n_q_points);


  // Velocity for BDF schemes
  this->previous_phase_order_values =
    std::vector<std::vector<double>>(maximum_number_of_previous_solutions(),
                                     std::vector<double>(n_q_points));

  this->previous_chemical_potential_values =
    std::vector<std::vector<double>>(maximum_number_of_previous_solutions(),
                                     std::vector<double>(n_q_points));
  // Velocity for SDIRK schemes
  this->stages_phase_order_values =
    std::vector<std::vector<double>>(max_number_of_intermediary_stages(),
                                     std::vector<double>(n_q_points));

  this->stages_chemical_potential_values =
    std::vector<std::vector<double>>(max_number_of_intermediary_stages(),
                                     std::vector<double>(n_q_points));

  // Initialize arrays related to shape functions
  // Phase-order shape functions
  this->phi_phase =
    std::vector<std::vector<double>>(n_q_points, std::vector<double>(n_dofs));
  this->grad_phi_phase = std::vector<std::vector<Tensor<1, dim>>>(
    n_q_points, std::vector<Tensor<1, dim>>(n_dofs));
  this->hess_phi_phase = std::vector<std::vector<Tensor<2, dim>>>(
    n_q_points, std::vector<Tensor<2, dim>>(n_dofs));
  this->laplacian_phi_phase =
    std::vector<std::vector<double>>(n_q_points, std::vector<double>(n_dofs));

  // Chemical potential shape functions
  this->phi_potential =
    std::vector<std::vector<double>>(n_q_points, std::vector<double>(n_dofs));
  this->grad_phi_potential = std::vector<std::vector<Tensor<1, dim>>>(
    n_q_points, std::vector<Tensor<1, dim>>(n_dofs));
  this->hess_phi_potential = std::vector<std::vector<Tensor<2, dim>>>(
    n_q_points, std::vector<Tensor<2, dim>>(n_dofs));
  this->laplacian_phi_potential =
    std::vector<std::vector<double>>(n_q_points, std::vector<double>(n_dofs));
}

//Must be completed to compute the physical properties everywhere depending on the phase order.
//template <int dim>
//void
//CahnHilliardScratchData<dim>::calculate_physical_properties()
//{
//  // Case where you have one fluid
//  switch (properties_manager.get_number_of_fluids())
//    {
//      case 1:
//        {
//          // In this case, only viscosity is the required property
//          const auto diffusivity_model =
//            properties_manager.get_tracer_diffusivity();
//          diffusivity_model->vector_value(fields, tracer_diffusivity);
//          break;
//        }
//      case 2:
//        {
//          // In this case,  we need both density and viscosity
//          const auto diffusivity_models =
//            properties_manager.get_tracer_diffusivity_vector();
//
//          diffusivity_models[0]->vector_value(fields, tracer_diffusivity_0);
//          diffusivity_models[1]->vector_value(fields, tracer_diffusivity_1);
//
//          // TODO Incomplete at the present time because the tracer VOF complete
//          // is not finished Blend the physical properties using the VOF field
//          for (unsigned int q = 0; q < this->n_q_points; ++q)
//            {
//              //          tracer_diffusivity[q] =
//              //            calculate_point_property(this->phase_values[q],
//              //                                     this->density_0[q],
//              //                                     this->density_1[q]);
//            }
//          break;
//        }
//      default:
//        throw std::runtime_error("Unsupported number of fluids (>2)");
//    }
//}


template class CahnHilliardScratchData<2>;
template class CahnHilliardScratchData<3>;
