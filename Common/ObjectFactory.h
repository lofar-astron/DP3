//# ObjectFactory.h: Generic object factory
//#
//# Copyright (C) 2006
//# ASTRON (Netherlands Institute for Radio Astronomy)
//# P.O.Box 2, 7990 AA Dwingeloo, The Netherlands
//#
//# This file is part of the LOFAR software suite.
//# The LOFAR software suite is free software: you can redistribute it and/or
//# modify it under the terms of the GNU General Public License as published
//# by the Free Software Foundation, either version 3 of the License, or
//# (at your option) any later version.
//#
//# The LOFAR software suite is distributed in the hope that it will be useful,
//# but WITHOUT ANY WARRANTY; without even the implied warranty of
//# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//# GNU General Public License for more details.
//#
//# You should have received a copy of the GNU General Public License along
//# with the LOFAR software suite. If not, see <http://www.gnu.org/licenses/>.
//#
//# $Id: ObjectFactory.h 32489 2015-09-26 18:38:13Z dijkema $

// \file
// Generic object factory.
//
// The implementation of this object factory was inspired by the article
// <em>Creating a Generic Object Factory</em>, written by Robert Geiman on
// GameDev.net
// (http://www.gamedev.net/reference/articles/article2097.asp). His
// implementation was taken as a starting point and adapted to work with the
// Boost.Preprocessor library. Some code refactoring was performed, and the
// code was made compliant with the %LOFAR coding standards.
//
// The structure of this header file may look intimidating at first; it
// performs almost magic. Keep in mind, though, that all \c BOOST_PP macros
// are simply used to generate repetitive code. In this case template
// specializations of the primary template
// \code 
// template <typename Signature, typename TypeId> class ObjectFactory;
// \endcode
// where \c Signature is the signature of the class' constructor, and \c
// TypeId is the type of its unique id. These template specializations define
// object factories for classes with up to \c OBJECT_FACTORY_MAX_CTOR_ARG
// constructor arguments.
//
// The preprocessor constructs below were taken from "Appendix A: An
// Introduction to Preprocessor Metaprogramming" in the book <em>C++ Template
// Metaprogramming</em>, by David Abrahams and Aleksey Gurtovoy -- available
// online at http://boost-consulting.com/tmpbook/preprocessor.html

#ifndef BOOST_PP_IS_ITERATING

#  ifndef LOFAR_COMMON_OBJECT_FACTORY_H
#    define LOFAR_COMMON_OBJECT_FACTORY_H

#    include <boost/preprocessor/repetition.hpp>
#    include <boost/preprocessor/iteration/iterate.hpp>

#include <map>
#include <vector>

#    ifndef OBJECT_FACTORY_MAX_CTOR_ARG
#      define OBJECT_FACTORY_MAX_CTOR_ARG 8
#    endif

namespace DP3
{
  // \addtogroup Common
  // @{
  // Primary template. \c Signature is the signature of the class'
  // constructor; \c TypeId is the type of its unique id.
  template <typename Signature, typename TypeId> class ObjectFactory;
  // @}
}

#    define BOOST_PP_ITERATION_LIMITS (0, OBJECT_FACTORY_MAX_CTOR_ARG)
#    define BOOST_PP_FILENAME_1 "ObjectFactory.h"
#    include BOOST_PP_ITERATE()
#  endif

#else

// \cond
#  define n BOOST_PP_ITERATION()
// \endcond

namespace DP3
{
  // \addtogroup Common
  // @{

  // Template specialization.
  // Expansion of the \c BOOST_PP_ENUM macros leads to the generation of the
  // following template specializations:
  // \code
  // template<typename Base, typename TypeId>
  // class ObjectFactory<Base (), TypeId>
  // \endcode
  // \code
  // template<typename Base, typename A0, typename TypeId>
  // class ObjectFactory<Base (A0), TypeId>
  // \endcode
  // \code
  // template<typename Base, typename A0, typename A1, typename TypeId>
  // class ObjectFactory<Base (A0, A1), TypeId>
  // \endcode
  // etc.
  //
  // The number of specializations generated depends on the value of the
  // preprocessor macro \c OBJECT_FACTORY_MAX_CTOR_ARG, which defines the
  // maximum number of constructor arguments for the class \c Base. It
  // defaults to 8.
  //
  template<typename Base BOOST_PP_ENUM_TRAILING_PARAMS(n, typename A), typename TypeId>
  class ObjectFactory<Base* (BOOST_PP_ENUM_PARAMS(n, A)), TypeId>
  {
  private:
    // Typedef for the function that creates an instance of a class that
    // inherits from \c Base.
    typedef Base* (*CreatorFunc)(BOOST_PP_ENUM_PARAMS(n, A));

    // Map associating \c TypeId, used to uniquely identify a class, and \c
    // CreatorFunc, a pointer to a function that creates an instance of this
    // class.
    typedef /*typename*/ std::map<TypeId, CreatorFunc> CreatorMap;

  public:
    typedef typename CreatorMap::const_iterator const_iterator;
    typedef typename CreatorMap::value_type value_type;

    // Register the class \c Derived using \a id as its unique identifier. The
    // static method doCreate<Derived>() will be registered as the method for
    // creating objects of type \c Derived. Registration will fail when \c
    // itsCreatorMap already contains a key with value \a id.
    // \return \c true is registration was successful; otherwise \c false.
    // \note \c Derived must inherit from the class \c Base. 
    template<typename Derived>
      bool registerClass(TypeId id)
    {
      return itsCreatorMap.insert(value_type(id, static_cast<CreatorFunc>(&doCreate<Derived>))).second;
    }
  
    // Unregister the class identified by \a id. The class identified by \a id
    // will be removed from \c itsCreatorMap.
    // \return \c true if deregistration was successful; otherwise \c false.
    bool unregisterClass(TypeId id)
    {
      return (itsCreatorMap.erase(id));
    }

    // Return a vector of TypeId containing the ID's of all registered
    // classes.
    std::vector<TypeId> registeredClassIds()
    {
      std::vector<TypeId> ids;
      const_iterator end = itsCreatorMap.end();
      for (const_iterator it = itsCreatorMap.begin(); it != end; ++it) {
	ids.push_back(it->first);
      }
      return ids;
    }

    // Create a new class instance of the class identified by \a id. 
    // If \a id is not found in \c itsCreatorMap a null pointer is returned.
    Base* create(TypeId id BOOST_PP_ENUM_TRAILING_BINARY_PARAMS(n, A, a))
    {
      const_iterator iter = itsCreatorMap.find(id);
      if (iter == itsCreatorMap.end()) return 0;
      else return iter->second(BOOST_PP_ENUM_PARAMS(n, a));
    }

    // Return a \c const iterator for the first element in \c itsCreatorMap.
    const_iterator begin() const
    {
      return itsCreatorMap.begin();
    }
  
    // Return a \c const iterator for the position after the last element in
    // \c itsCreatorMap.
    const_iterator end() const
    {
      return itsCreatorMap.end();
    }

  private:

    // Map binding the unique class id and its creator function.
    CreatorMap itsCreatorMap;

    // Return a pointer to a new instance of class \c Derived. This is the
    // method that's being registered as the creator function when a class of
    // type \c Derived registers itself with the object factory.
    // \note \c Derived must inherit from \c Base.
    template<typename Derived>
      static Base* doCreate(BOOST_PP_ENUM_BINARY_PARAMS(n, A, a))
    {
      return new Derived (BOOST_PP_ENUM_PARAMS(n, a));
    }

  };

  // @}

} // namespace LOFAR

#  undef n

#endif
