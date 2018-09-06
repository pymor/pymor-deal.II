// This file is part of the pyMOR project (http://www.pymor.org).
// Copyright 2013-2018 pyMOR developers and contributors. All rights reserved.
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef RHS_HH
#define RHS_HH

namespace dealii {
// @sect3{Right hand side values}

// Before going over to the implementation of the main class, we declare and
// define the class which describes the right hand side. This time, the
// right hand side is vector-valued, as is the solution, so we will describe
// the changes required for this in some more detail.
//
// The first thing is that vector-valued functions have to have a
// constructor, since they need to pass down to the base class of how many
// components the function consists. The default value in the constructor of
// the base class is one (i.e.: a scalar function), which is why we did not
// need not define a constructor for the scalar function used in previous
// programs.
template <int dim>
class RightHandSide : public Function<dim> {
public:
  RightHandSide();

  // The next change is that we want a replacement for the
  // <code>value</code> function of the previous examples. There, a second
  // parameter <code>component</code> was given, which denoted which
  // component was requested. Here, we implement a function that returns the
  // whole vector of values at the given place at once, in the second
  // argument of the function. The obvious name for such a replacement
  // function is <code>vector_value</code>.
  //
  // Secondly, in analogy to the <code>value_list</code> function, there is
  // a function <code>vector_value_list</code>, which returns the values of
  // the vector-valued function at several points at once:
  virtual void vector_value(const Point<dim>& p, Vector<double>& values) const;

  virtual void vector_value_list(const std::vector<Point<dim>>& points, std::vector<Vector<double>>& value_list) const;
};

// This is the constructor of the right hand side class. As said above, it
// only passes down to the base class the number of components, which is
// <code>dim</code> in the present case (one force component in each of the
// <code>dim</code> space directions).
//
// Some people would have moved the definition of such a short function
// right into the class declaration. We do not do that, as a matter of
// style: the deal.II style guides require that class declarations contain
// only declarations, and that definitions are always to be found
// outside. This is, obviously, as much as matter of taste as indentation,
// but we try to be consistent in this direction.
template <int dim>
RightHandSide<dim>::RightHandSide()
  : Function<dim>(dim) {}

// Next the function that returns the whole vector of values at the point
// <code>p</code> at once.
//
// To prevent cases where the return vector has not previously been set to
// the right size we test for this case and otherwise throw an exception at
// the beginning of the function. Note that enforcing that output arguments
// already have the correct size is a convention in deal.II, and enforced
// almost everywhere. The reason is that we would otherwise have to check at
// the beginning of the function and possibly change the size of the output
// vector. This is expensive, and would almost always be unnecessary (the
// first call to the function would set the vector to the right size, and
// subsequent calls would only have to do redundant checks). In addition,
// checking and possibly resizing the vector is an operation that can not be
// removed if we can't rely on the assumption that the vector already has
// the correct size; this is in contract to the <code>Assert</code> call
// that is completely removed if the program is compiled in optimized mode.
//
// Likewise, if by some accident someone tried to compile and run the
// program in only one space dimension (in which the elastic equations do
// not make much sense since they reduce to the ordinary Laplace equation),
// we terminate the program in the second assertion. The program will work
// just fine in 3d, however.
template <int dim>
inline void RightHandSide<dim>::vector_value(const Point<dim>& p, Vector<double>& values) const {
  Assert(values.size() == dim, ExcDimensionMismatch(values.size(), dim));
  Assert(dim >= 2, ExcNotImplemented());

  // The rest of the function implements computing force values. We will use
  // a constant (unit) force in x-direction located in two little circles
  // (or spheres, in 3d) around points (0.5,0) and (-0.5,0), and y-force in
  // an area around the origin; in 3d, the z-component of these centers is
  // zero as well.
  //
  // For this, let us first define two objects that denote the centers of
  // these areas. Note that upon construction of the <code>Point</code>
  // objects, all components are set to zero.
  Point<dim> point_1, point_2;
  point_1(0) = 0.5;
  point_2(0) = -0.5;

  // If now the point <code>p</code> is in a circle (sphere) of radius 0.2
  // around one of these points, then set the force in x-direction to one,
  // otherwise to zero:
  if ((static_cast<Point<dim>>(p - point_1).square() < 0.2 * 0.2)
      || (static_cast<Point<dim>>(p - point_2).norm_square() < 0.2 * 0.2))
    values(0) = 1;
  else
    values(0) = 0;

  // Likewise, if <code>p</code> is in the vicinity of the origin, then set
  // the y-force to 1, otherwise to zero:
  if (p.square() < 0.2 * 0.2)
    values(1) = 1;
  else
    values(1) = 0;
}

// Now, this is the function of the right hand side class that returns the
// values at several points at once. The function starts out with checking
// that the number of input and output arguments is equal (the sizes of the
// individual output vectors will be checked in the function that we call
// further down below). Next, we define an abbreviation for the number of
// points which we shall work on, to make some things simpler below.
template <int dim>
void RightHandSide<dim>::vector_value_list(const std::vector<Point<dim>>& points,
                                           std::vector<Vector<double>>& value_list) const {
  Assert(value_list.size() == points.size(), ExcDimensionMismatch(value_list.size(), points.size()));

  const unsigned int n_points = points.size();

  // Finally we treat each of the points. In one of the previous examples,
  // we have explained why the
  // <code>value_list</code>/<code>vector_value_list</code> function had
  // been introduced: to prevent us from calling virtual functions too
  // frequently. On the other hand, we now need to implement the same
  // function twice, which can lead to confusion if one function is changed
  // but the other is not.
  //
  // We can prevent this situation by calling
  // <code>RightHandSide::vector_value</code> on each point in the input
  // list. Note that by giving the full name of the function, including the
  // class name, we instruct the compiler to explicitly call this function,
  // and not to use the virtual function call mechanism that would be used
  // if we had just called <code>vector_value</code>. This is important,
  // since the compiler generally can't make any assumptions which function
  // is called when using virtual functions, and it therefore can't inline
  // the called function into the site of the call. On the contrary, here we
  // give the fully qualified name, which bypasses the virtual function
  // call, and consequently the compiler knows exactly which function is
  // called and will inline above function into the present location. (Note
  // that we have declared the <code>vector_value</code> function above
  // <code>inline</code>, though modern compilers are also able to inline
  // functions even if they have not been declared as inline).
  //
  // It is worth noting why we go to such length explaining what we
  // do. Using this construct, we manage to avoid any inconsistency: if we
  // want to change the right hand side function, it would be difficult to
  // always remember that we always have to change two functions in the same
  // way. Using this forwarding mechanism, we only have to change a single
  // place (the <code>vector_value</code> function), and the second place
  // (the <code>vector_value_list</code> function) will always be consistent
  // with it. At the same time, using virtual function call bypassing, the
  // code is no less efficient than if we had written it twice in the first
  // place:
  for (unsigned int p = 0; p < n_points; ++p)
    RightHandSide<dim>::vector_value(points[p], value_list[p]);
}

} // namespace dealii
#endif // RHS_HH
