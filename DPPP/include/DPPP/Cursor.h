//# Cursor.h: Multi-dimensional iterators.
//#
//# Copyright (C) 2012
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
//# $Id$

#ifndef DPPP_CURSOR_H
#define DPPP_CURSOR_H

// \file
// Multi-dimensional iterators.

#include <Common/lofar_algorithm.h>
#include <Common/LofarLogger.h>

namespace LOFAR
{
namespace DPPP
{

// \addtogroup NDPPP
// @{

template <typename T>
class cursor
{
public:
    cursor()
        :   itsPointer(0),
            itsRank(0)
    {
    }

    cursor(T *pointer)
        :   itsPointer(pointer),
            itsRank(1)
    {
        fill(itsStrides, itsStrides + MAX_RANK, 0);
        itsStrides[0] = 1;
    }

    template <typename T_STRIDE>
    cursor(T *pointer, size_t rank, const T_STRIDE *strides)
        :   itsPointer(pointer),
            itsRank(rank)
    {
        DBGASSERT (rank <= MAX_RANK);
        copy(strides, strides + itsRank, itsStrides);
        fill(itsStrides + itsRank, itsStrides + MAX_RANK, 0);
    }

    size_t rank() const
    {
        return itsRank;
    }

    cursor &operator++()
    {
        itsPointer += itsStrides[0];
        return *this;
    }

    cursor operator++(int)
    {
        cursor tmp = *this;
        itsPointer += itsStrides[0];
        return tmp;
    }

    cursor &operator+=(size_t n)
    {
        itsPointer += n * itsStrides[0];
        return *this;
    }

    cursor &operator-=(size_t n)
    {
        itsPointer -= n * itsStrides[0];
        return *this;
    }

    T &operator*()
    {
        return *itsPointer;
    }

    const T &operator*() const
    {
        return *itsPointer;
    }

    T *operator->()
    {
        return itsPointer;
    }

    const T *operator->() const
    {
        return itsPointer;
    }

    T &operator[](size_t n)
    {
        return *(itsPointer + n * itsStrides[0]);
    }

    const T &operator[](size_t n) const
    {
        return *(itsPointer + n * itsStrides[0]);
    }

    void forward(size_t i)
    {
        itsPointer += itsStrides[i];
    }

    void forward(size_t i, size_t n)
    {
        itsPointer += n * itsStrides[i];
    }

    void backward(size_t i)
    {
        itsPointer -= itsStrides[i];
    }

    void backward(size_t i, size_t n)
    {
        itsPointer -= n * itsStrides[i];
    }

    T *address()
    {
        return itsPointer;
    }

    const T *address() const
    {
        return itsPointer;
    }

    size_t stride(size_t i) const
    {
        return itsStrides[i];
    }

private:
    static const size_t MAX_RANK = 5;

    T*      itsPointer;
    size_t  itsRank;
    size_t  itsStrides[MAX_RANK];
};

template <typename T>
class const_cursor
{
public:
    const_cursor()
        :   itsPointer(0),
            itsRank(0)
    {
    }

    const_cursor(const T *pointer)
        :   itsPointer(pointer),
            itsRank(1)
    {
        fill(itsStrides, itsStrides + MAX_RANK, 0);
        itsStrides[0] = 1;
    }

    template <typename T_STRIDE>
    const_cursor(const T *pointer, size_t rank, const T_STRIDE *strides)
        :   itsPointer(pointer),
            itsRank(rank)
    {
        copy(strides, strides + itsRank, itsStrides);
        fill(itsStrides + itsRank, itsStrides + MAX_RANK, 0);
    }

    const_cursor(const cursor<T> &other)
        :   itsPointer(other.address()),
            itsRank(other.rank())
    {
        for(size_t i = 0; i < itsRank; ++i)
        {
            itsStrides[i] = other.stride(i);
        }
        fill(itsStrides + itsRank, itsStrides + MAX_RANK, 0);
    }

    size_t rank() const
    {
        return itsRank;
    }

    const_cursor &operator++()
    {
        itsPointer += itsStrides[0];
        return *this;
    }

    const_cursor operator++(int)
    {
        const_cursor tmp = *this;
        itsPointer += itsStrides[0];
        return tmp;
    }

    const_cursor &operator+=(size_t n)
    {
        itsPointer += n * itsStrides[0];
        return *this;
    }

    const_cursor &operator-=(size_t n)
    {
        itsPointer -= n * itsStrides[0];
        return *this;
    }

    const T &operator*() const
    {
        return *itsPointer;
    }

    const T *operator->() const
    {
        return itsPointer;
    }

    const T &operator[](size_t n) const
    {
        return *(itsPointer + n * itsStrides[0]);
    }

    void forward(size_t i)
    {
        itsPointer += itsStrides[i];
    }

    void forward(size_t i, size_t n)
    {
        itsPointer += n * itsStrides[i];
    }

    void backward(size_t i)
    {
        itsPointer -= itsStrides[i];
    }

    void backward(size_t i, size_t n)
    {
        itsPointer -= n * itsStrides[i];
    }

    const T *address() const
    {
        return itsPointer;
    }

    size_t stride(size_t i) const
    {
        return itsStrides[i];
    }

private:
    static const size_t MAX_RANK = 5;

    const T*    itsPointer;
    size_t      itsRank;
    size_t      itsStrides[MAX_RANK];
};

// @}

} //# namespace DPPP
} //# namespace LOFAR

#endif
