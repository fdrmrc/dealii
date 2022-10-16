// ---------------------------------------------------------------------
//
// Copyright (C) 1999 - 2022 by the deal.II authors
//
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

#include <deal.II/cgal/bounding_box_cgal.h>



DEAL_II_NAMESPACE_OPEN
template <int spacedim, typename Number>
OptimalBoundingBox<spacedim, Number>::OptimalBoundingBox(
  const std::vector<Point<spacedim, Number>> &points)
{
#if DEAL_II_CGAL_VERSION_GTE(5, 1, 5)
  Assert(points.size() == Utilities::pow(2, spacedim),
         ExcMessage("Invalid number of points."));
  Assert(spacedim == 3, ExcNotImplemented("Not implemented in 1D and 2D."));
  using K          = CGAL::Exact_predicates_inexact_constructions_kernel;
  using CGALPoint3 = K::Point_3;

  std::vector<CGALPoint3>                             cgal_pts(points.size());
  std::array<CGALPoint3, Utilities::pow(2, spacedim)> cgal_out_pts;
  std::transform(points.begin(),
                 points.end(),
                 cgal_pts.begin(),
                 [&](const Point<spacedim, Number> &p) {
                   return CGALWrappers::dealii_point_to_cgal_point<CGALPoint3>(
                     p);
                 });
  Custom_traits_BBox<K> custom_traits;
  std::cout << "Pre crezione" << std::endl;
  CGAL::oriented_bounding_box(cgal_pts,
                              cgal_out_pts,
                              CGAL::parameters::geom_traits(custom_traits));
  std::cout << "Post crezione" << std::endl;
  std::transform(cgal_out_pts.begin(),
                 cgal_out_pts.end(),
                 extreme_pts.begin(),
                 [&](const CGALPoint3 &p) {
                   return CGALWrappers::cgal_point_to_dealii_point<spacedim>(p);
                 });


#else
  Assert(false, ExcNeedsCGAL());
#endif
}



template class OptimalBoundingBox<3, double>;
// template class OptimalBoundingBox<2, double>;
// template class OptimalBoundingBox<1, double>;

DEAL_II_NAMESPACE_CLOSE