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
  {
    std::cout << "Le entrà!" << std::endl;
  }
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
    {
      std::cout << "standard cstr matrix" << std::endl;
    }

    Matrix()
      : dealii::FullMatrix<FT>(3, 3) // TODO just testing
    {
      std::cout << "default matrix" << std::endl;
    }

    inline std::size_t
    number_of_rows() const
    {
      std::cout << "righe" << std::endl;

      return this->m();
    }

    inline std::size_t
    number_of_columns() const
    {
      std::cout << "colonne " << std::endl;

      return this->n();
    }

    Matrix
    operator*(const Matrix &B)
    {
      std::cout << "* Matrix in teoria è ok" << std::endl;

      dealii::FullMatrix<FT> C(B.number_of_rows(), this->number_of_columns());
      this->mmult(C, B);
      return *this;
    }


    typename Custom_traits_BBox<Kernel>::Vector
    operator*(const typename Custom_traits_BBox<Kernel>::Vector &v) const
    {
      deallog << "Matrix*vector" << std::endl;

      Custom_traits_BBox<Kernel>::Vector ret(this->number_of_rows());
      this->vmult(ret, v);
      return ret;
    }


    Matrix
    operator+(const Matrix &B) const
    {
      std::cout << "Devo ancora implementare + Matrix" << std::endl;
      Matrix ret(B.number_of_rows(), B.number_of_columns());
      for (unsigned int i = 0; i < B.number_of_rows(); ++i)
        {
          for (unsigned int j = 0; j < B.number_of_columns(); ++j)
            {
              ret.set(i, j, (*this)(i, j) + B(i, j));
            }
        }

      return ret;
    }
  };

  // Vector type
  class Vector : public dealii::Vector<FT>
  {
  public:
    Vector()
      : dealii::Vector<FT>(3)
    {
      std::cout << "Default vector." << std::endl;
    }

    Vector(std::size_t n)
    {
      std::cout << "Devo costruttore n vector" << std::endl;

      this->reinit(n);
    }

    // inline std::size_t
    // size()
    // {
    //   std::cout << "Devo Size vector" << std::endl;

    //   return this->size();
    // }

    // inline FT
    // operator()(std::size_t i) const
    // {
    //   std::cout << "Denttro FT ()" << std::endl;

    //   return this->operator()(i);
    // }

    void
    set(std::size_t i, const FT value)
    {
      std::cout << "dentro set" << std::endl;

      this->operator[](i) = value;
    }

    FT *
    vector()
    {
      std::cout << "Returing as array" << std::endl;
      FT *vec_as_array = new FT[this->size()];
      for (unsigned int i = 0; i < this->size(); ++i)
        {
          vec_as_array[i] = (*this)(i);
        }
      return vec_as_array;
    }
  };


public:
  static Matrix
  get_Q(const Matrix &m)
  {
    std::cout << "Entered in get_Q" << std::endl;
    // Assert(false, ExcMessage("Hello my friend."));
    dealii::QR<dealii::Vector<FT>> qr;

    std::cout << "Here I am" << std::endl;
    Vector v0(3);
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
    [[maybe_unused]] bool           c0 = qr.append_column(v0);
    [[maybe_unused]] bool           c1 = qr.append_column(v1);
    [[maybe_unused]] bool           c2 = qr.append_column(v2);
    std::vector<dealii::Vector<FT>> Q(
      qr.size()); // defaults cstr (empty vector)
    Matrix Q_matrix(3, qr.size());
    for (unsigned int j = 0; j < qr.size(); ++j)
      {
        dealii::Vector<FT> x(qr.size());
        x = 0; //
        // x.reinit(3);
        std::cout << "Prima di accedere a j." << std::endl;
        x[j] = 1.;
        std::cout << "Dopo acceduto a j." << std::endl;
        Q[j].reinit(3);
        std::cout << "Here I am also." << std::endl;
        qr.multiply_with_Q(Q[j], x);
      }

    for (unsigned int i = 0; i < 3; ++i)
      {
        for (unsigned int j = 0; j < qr.size(); ++j)
          {
            Q_matrix.set(i, j, Q[i][j]);
          }
      }
    return Q_matrix;
  }
};



template <typename Matrix = typename Custom_traits_BBox<
            CGAL::Exact_predicates_inexact_constructions_kernel>::Matrix,
          typename FT = double>
Matrix
operator*(const double &scalar, const Matrix &B)
{
  deallog << "Scalar from left" << std::endl;
  Matrix ret(B.number_of_rows(), B.number_of_columns());

  for (unsigned int i = 0; i < B.number_of_rows(); ++i)
    {
      for (unsigned int j = 0; j < B.number_of_columns(); ++j)
        {
          const double Bij        = B(i, j);
          const double Bij_scalar = Bij * scalar; // TODO scalar is missing
          ret.set(i, j, Bij_scalar);
        }
    }
  return ret;
}

template <typename Matrix = typename Custom_traits_BBox<
            CGAL::Exact_predicates_inexact_constructions_kernel>::Matrix>
Matrix
operator+(const Matrix &A, const Matrix &B)
{
  std::cout << " + Matrix" << std::endl;
  Matrix ret(B.number_of_rows(), B.number_of_columns());
  // for (unsigned int i = 0; i < B.number_of_rows(); ++i)
  //   {
  //     for (unsigned int j = 0; j < B.number_of_columns(); ++j)
  //       {
  //         ret.set(i, j, A(i, j) + B(i, j));
  //       }
  //   }
  A.add(ret, B);
  return ret;
}

typename Custom_traits_BBox<
  CGAL::Exact_predicates_inexact_constructions_kernel>::Matrix
operator*(const typename Custom_traits_BBox<
            CGAL::Exact_predicates_inexact_constructions_kernel>::Matrix &A,
          const typename Custom_traits_BBox<
            CGAL::Exact_predicates_inexact_constructions_kernel>::Matrix &B)
{
  std::cout << "Matrix" << std::endl;
  typename Custom_traits_BBox<
    CGAL::Exact_predicates_inexact_constructions_kernel>::Matrix
    ret(B.number_of_rows(), B.number_of_columns());
  A.mmult(ret, B);
  // TODO: use matrix-matrix multiplication
  return ret;
}

// Controlla numero di righe e colonne nelle varia operazioni... magari sta
// sforando nell'algoritmo.

template <int spacedim, typename Number = double>
class OptimalBoundingBox
{
public:
  OptimalBoundingBox() = default;

  OptimalBoundingBox(const std::vector<Point<spacedim, Number>> &points);

  std::array<Point<spacedim, Number>, Utilities::pow(2, spacedim)>
  get_extreme_pts()
  {
    return extreme_pts;
  }

private:
  std::array<Point<spacedim, Number>, Utilities::pow(2, spacedim)> extreme_pts;
};

DEAL_II_NAMESPACE_CLOSE

#endif
