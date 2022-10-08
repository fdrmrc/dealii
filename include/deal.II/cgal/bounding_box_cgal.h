// ---------------------------------------------------------------------
//
// Copyright (C) 2017 - 2021 by the deal.II authors
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

#ifndef dealii_base_bounding_box_cgal_h
#define dealii_base_bounding_box_cgal_h

#include <deal.II/base/bounding_box.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/base/point.h>
#include <deal.II/base/utilities.h>

#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/qr.h>

#include <CGAL/Aff_transformation_3.h>
#include <CGAL/Optimal_bounding_box/oriented_bounding_box.h>
#include <deal.II/cgal/point_conversion.h>

#include <algorithm>
#include <limits>

DEAL_II_NAMESPACE_OPEN



template <typename Kernel>
class Custom_traits_BBox : public Kernel
{
public:
  Custom_traits_BBox()
  {}
  /// The field number type
  using FT = typename Kernel::FT;

  /// The affine transformation type
  using Aff_transformation_3 = typename CGAL::Aff_transformation_3<Kernel>;

  // Matrix type
  class Matrix : public dealii::FullMatrix<FT>
  {
  public:
    Matrix(const std::size_t n_rows, const std::size_t n_cols)
      : dealii::FullMatrix<FT>(n_rows, n_cols)
    {}

    Matrix()
      : dealii::FullMatrix<FT>()
    {}

    inline std::size_t
    number_of_rows() const
    {
      return this->m();
    }

    inline std::size_t
    number_of_columns() const
    {
      return this->n();
    }

    Matrix
    operator*(const Matrix &B) const
    {
      dealii::FullMatrix<FT> C(B.number_of_rows(), this->number_of_columns());
      this->mmult(C, B);
      return *this; // TODO: convert properly
    }


    typename Custom_traits_BBox<Kernel>::Vector
    operator*(const typename Custom_traits_BBox<Kernel>::Vector &B) const
    {
      // dealii::FullMatrix<FT> C(B.number_of_rows(),
      // this->number_of_columns()); this->mmult(C, B);
      typename Custom_traits_BBox<Kernel>::Vector v;
      return v; // TODO: convert properly
    }


    Matrix
    operator+(const Matrix &B) const
    {
      // FullMatrix<FT> a;
      return *this; // TODO: convert properly
    }
  };

  // Vector type
  class Vector : public dealii::Vector<FT>
  {
  public:
    Vector() = default;

    Vector(std::size_t n)
    {
      this->reinit(n);
    }

    inline std::size_t
    size()
    {
      return this->size();
    }

    inline FT
    operator()(std::size_t i) const
    {
      return this->operator()(i);
    }

    void
    set(std::size_t i, const FT value)
    {
      this->operator[](i) = value;
    }

    FT *
    vector()
    {}
  };


public:
  static Matrix
  get_Q(const Matrix &m)
  {
    dealii::QR<Vector> qr;
    Vector             v0(3);
    v0.set(0, m(0, 0));
    v0.set(1, m(1, 0));
    v0.set(2, m(2, 0));
    Vector v1(3);
    v1.set(0, m(0, 1));
    v1.set(1, m(1, 1));
    v1.set(2, m(2, 1));
    Vector v2(3);
    v2.set(0, m(0, 2));
    v2.set(1, m(1, 2));
    v2.set(2, m(2, 2));
    [[maybe_unused]] bool c0 = qr.append_column(v0);
    [[maybe_unused]] bool c1 = qr.append_column(v1);
    [[maybe_unused]] bool c2 = qr.append_column(v2);
    std::vector<Vector>   Q(3); // here there's no defaults cstr for the
    // inner type Vector
    Matrix Q_matrix(3, 3);
    for (unsigned int j = 0; j < 3; ++j)
      {
        Vector x(3);
        x    = 0;
        x[j] = 1.;
        Q[j].reinit(3);
        qr.multiply_with_Q(Q[j], x);
      }

    for (unsigned int i = 0; i < 3; ++i)
      {
        for (unsigned int j = 0; j < 3; ++j)
          {
            Q_matrix(i, j) = Q[i][j];
          }
      }
    return Q_matrix;
  }
};



template <typename FT>
typename Custom_traits_BBox<
  CGAL::Exact_predicates_inexact_constructions_kernel>::Matrix
operator*(const FT &scalar,
          const typename Custom_traits_BBox<
            CGAL::Exact_predicates_inexact_constructions_kernel>::Matrix &rhs)
{
  typename Custom_traits_BBox<
    CGAL::Exact_predicates_inexact_constructions_kernel>::Matrix ciccio;
  return ciccio;
  // return rhs *= scalar; // calls rhs.operator*=(scalar);
}

template <int spacedim, typename Number = double>
class OptimalBoundingBox : public BoundingBox<spacedim, Number>
{
public:
  OptimalBoundingBox() = default;

  OptimalBoundingBox(const std::vector<Point<spacedim, Number>> &points);

private:
  std::array<Point<spacedim, Number>, Utilities::pow(2, spacedim)> extreme_pts;
};

DEAL_II_NAMESPACE_CLOSE

#endif
