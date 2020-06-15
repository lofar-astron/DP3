// Barrier.h: Barrier operation for synchronizing threads.
//
// Copyright (C) 2020
// ASTRON (Netherlands Institute for Radio Astronomy)
// P.O.Box 2, 7990 AA Dwingeloo, The Netherlands
//
// This file is part of the LOFAR software suite.
// The LOFAR software suite is free software: you can redistribute it and/or
// modify it under the terms of the GNU General Public License as published
// by the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// The LOFAR software suite is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License along
// with the LOFAR software suite. If not, see <http://www.gnu.org/licenses/>.

/// \file
/// \brief Barrier operation for synchronizing threads.

#ifndef BARRIER_H
#define BARRIER_H

#include <condition_variable>
#include <functional>
#include <mutex>

namespace DP3 {
  
  /**
  * This class is unfortunately necessary because boost::barrier had a 
  * bug in completion functions ( https://svn.boost.org/trac10/ticket/13561 ),
  * and std::barrier is still experimental.
  */
  class Barrier
  {
  public:
    /**
    * Construct barrier for n threads with the given completion function.
    * @param n Number of threads to wait for
    * @param completionFunction void function that is called when all threads have
    * arrived, just before the threads are released.
    */
    Barrier(size_t n, std::function<void()> completionFunction) : _n(n), _count(_n), _cycle(0), _completionFunction(completionFunction)
    {
    }
    
    Barrier& operator=(Barrier&& rhs)
    {
      /// count could be checked for != _n to assert no threads
      /// are waiting, but I leave this the responsibility of the caller.
      _n = rhs._n;
      _count = _n;
      _cycle = 0;
      _completionFunction = rhs._completionFunction;
      return *this;
    }
    
    /**
    * Wait until all threads are waiting for the barrier.
    */
    void wait()
    {
      std::unique_lock<std::mutex> lock(_mutex);
      --_count;
      if(_count == 0)
      {
        ++_cycle;
        _count = _n;
        _completionFunction();
        _condition.notify_all();
      }
      else {
        size_t cycle = _cycle;
        while(cycle == _cycle) 
          _condition.wait(lock);
      }
    }
    
  private:
    std::mutex _mutex;
    std::condition_variable _condition;
    size_t _n, _count, _cycle;
    std::function<void()> _completionFunction;
  };
}

#endif
