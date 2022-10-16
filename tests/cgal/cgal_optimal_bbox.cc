// ---------------------------------------------------------------------

// Copyright (C) 2022 by the deal.II authors

// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------

// Compute intersection of a square with a cube and return a Quadrature rule
// over its intersection.
#include <deal.II/fe/mapping_q1.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>

#include <deal.II/cgal/bounding_box_cgal.h>

#include "../tests.h"

int
main()
{
  initlog();
  Triangulation<3> tria;
  GridGenerator::hyper_cube(tria, -2, 2.5);
  GridTools::rotate(numbers::PI_4, 2, tria);



  using K = CGAL::Exact_predicates_inexact_constructions_kernel;
  // typename Custom_traits_BBox<K>::Vector x(3);
  // x.reinit(3);
  // x    = 0;
  // x[0] = 1.;


  OptimalBoundingBox<3> opt_bbox(tria.get_vertices());
  for (const auto &p : tria.get_vertices())
    deallog << p << std::endl;
  deallog << "Risultati:" << std::endl;
  GridOut       go;
  std::ofstream out("temp.vtk");
  go.write_vtk(tria, out);
  for (const auto &extre_pt : opt_bbox.get_extreme_pts())
    {
      deallog << extre_pt << std::endl;
    }



  // typename Custom_traits_BBox<K>::Vector test_vector(5);
  // // test_vector(2) = 1;
  // test_vector.set(2, 96);
  // test_vector[3] = 1;
  // test_vector(1) = 9000;
  // deallog << "La size è: " << test_vector.size() << std::endl;
  // for (const auto x : test_vector)
  //   deallog << x << std::endl;


  typename Custom_traits_BBox<K>::Matrix m(3, 3);
  typename Custom_traits_BBox<K>::Matrix n(3, 3);
  auto                                   mn = m * n;
  deallog << mn.number_of_columns() << std::endl;

  for (size_t i = 0; i < m.number_of_rows(); i++)
    {
      for (size_t j = 0; j < m.number_of_columns(); j++)
        {
          mn.set(i, j, i + i);
          deallog << mn(i, j) << "\t" << std::endl;
        }
      deallog << std::endl;
    }

  typename Custom_traits_BBox<K>::Matrix M(mn);

  for (size_t i = 0; i < m.number_of_rows(); i++)
    {
      for (size_t j = 0; j < m.number_of_columns(); j++)
        {
          deallog << "M(" << i << "," << j << "=)" << M(i, j) << "\t"
                  << std::endl;
        }
      deallog << std::endl;
    }

  // M.add(2, mn);
  auto test = M * M;
  for (size_t i = 0; i < test.number_of_rows(); i++)
    {
      for (size_t j = 0; j < test.number_of_columns(); j++)
        {
          deallog << " test(i, j) = " << test(i, j) << "\t" << std::endl;
        }
      deallog << std::endl;
    }

  // auto altrotest = M + M;
  // for (size_t i = 0; i < M.number_of_rows(); i++)
  //   {
  //     for (size_t j = 0; j < M.number_of_columns(); j++)
  //       {
  //         deallog << " altrotest(i, j) = " << altrotest(i, j) << "\t"
  //                 << std::endl;
  //       }
  //     deallog << std::endl;
  //   }

  // typename Custom_traits_BBox<K>::Vector miovec(altrotest.number_of_rows());
  // auto                                   vv = altrotest * miovec;


  // altrotest = mn;
  // deallog << test_vector.norm_sqr() << std::endl;
}
