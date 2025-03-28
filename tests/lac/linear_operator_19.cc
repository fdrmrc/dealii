// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2015 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// tests Tvmult() for two of the Trilinos direct solvers (KLU, MUMPS) that
// support it.

#include "deal.II/base/exception_macros.h"
#include <deal.II/base/exceptions.h>
#include <deal.II/base/function.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/generic_linear_algebra.h>
#include <deal.II/lac/linear_operator_tools.h>
#include <deal.II/lac/packaged_operation.h>
#include <deal.II/lac/trilinos_linear_operator.h>
#include <deal.II/lac/trilinos_precondition.h>
#include <deal.II/lac/trilinos_solver.h>
#include <deal.II/lac/trilinos_sparse_matrix.h>
#include <deal.II/lac/vector.h>

#include <deal.II/numerics/vector_tools.h>

#include "../tests.h"

using namespace dealii;



void
test(const unsigned int size)
{
  // Create a sparsity pattern
  TrilinosWrappers::SparsityPattern sparsity_pattern(size, size, size);
  for (unsigned int i = 0; i < size; ++i)
    {
      sparsity_pattern.add(i, i);
      if (i > 0)
        sparsity_pattern.add(i, i - 1);
      if (i < size - 1)
        sparsity_pattern.add(i, i + 1);
    }
  sparsity_pattern.compress();

  // Define the matrix
  TrilinosWrappers::SparseMatrix A;
  A.reinit(sparsity_pattern);

  for (unsigned int i = 0; i < size; ++i)
    {
      A.set(i, i, 2.0);
      if (i < size - 1)
        {
          A.set(i, i + 1, -1.0);
          A.set(i + 1, i, 2.0);
        }
    }

  A.compress(VectorOperation::insert);


  using VectorType = TrilinosWrappers::MPI::Vector;

  VectorType solution_transpose, residual;
  residual.reinit(complete_index_set(size));
  solution_transpose.reinit(complete_index_set(size));


  TrilinosWrappers::MPI::Vector system_rhs;
  system_rhs.reinit(complete_index_set(size));
  system_rhs = 1.0;
  system_rhs.compress(VectorOperation::insert);



  for (const std::string direct_solver_name : {"Amesos_Klu", "Amesos_Mumps"})
    {
      SolverControl solver_control(1000, 1e-15);
      {
        deallog << direct_solver_name << std::endl;
        residual = 0;

        TrilinosWrappers::SolverDirect::AdditionalData data;
        data.solver_type = direct_solver_name;
        TrilinosWrappers::SolverDirect solver(solver_control, data);
        solver.initialize(A);
        // First, solve A^Ty = b for y
        solver.Tvmult(solution_transpose, system_rhs);
        // Then, check the residual
        A.Tvmult(residual, solution_transpose);
        residual -= system_rhs;
        AssertThrow(residual.l2_norm() < 1e-14,
                    ExcMessage("Residual too large."));
        deallog << "OK" << std::endl;
      }

      // Do the same, but with a LinearOperator,

      {
        deallog << "LinearOperator with " + direct_solver_name << std::endl;
        residual           = 0.;
        solution_transpose = 0.;
        SolverControl solver_control(1000, 1e-15);
        TrilinosWrappers::SolverDirect::AdditionalData data;
        data.solver_type = direct_solver_name;
        TrilinosWrappers::SolverDirect direct_solver(solver_control, data);
        direct_solver.initialize(A);

        using PayloadType = dealii::TrilinosWrappers::internal::
          LinearOperatorImplementation::TrilinosPayload;
        auto lo_A = linear_operator<VectorType, VectorType, PayloadType>(A);
        auto lo_A_inv =
          linear_operator<VectorType, VectorType, PayloadType>(lo_A,
                                                               direct_solver);
        auto lo_A_inv_transpose = transpose_operator(lo_A_inv);
        solution_transpose      = lo_A_inv_transpose * system_rhs;
        lo_A.Tvmult(residual, solution_transpose);
        residual -= system_rhs;
        AssertThrow(residual.l2_norm() < 1e-14,
                    ExcMessage("Residual too large."));
        deallog << "OK (with LinearOperator)" << std::endl;
        // deallog << "Norm of residual: " << residual.l2_norm() << std::endl;
      }
    }
}


int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  initlog();
  deallog << std::setprecision(10);

  try
    {
      const unsigned int size = 5000;
      test(size);
    }
  catch (const std::exception &exc)
    {
      deallog << std::endl
              << std::endl
              << "----------------------------------------------------"
              << std::endl;
      deallog << "Exception on processing: " << std::endl
              << exc.what() << std::endl
              << "Aborting!" << std::endl
              << "----------------------------------------------------"
              << std::endl;

      return 1;
    }
  catch (...)
    {
      deallog << std::endl
              << std::endl
              << "----------------------------------------------------"
              << std::endl;
      deallog << "Unknown exception!" << std::endl
              << "Aborting!" << std::endl
              << "----------------------------------------------------"
              << std::endl;
      return 1;
    };
}