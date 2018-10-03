#ifndef PARALLEL_FOR_H
#define PARALLEL_FOR_H

#include <cstring>
#include <thread>
#include <mutex>
#include <vector>

#include <sched.h>

namespace DP3 {

  template<typename Iter>
  class ParallelFor
  {
  public:
    ParallelFor(size_t nThreads) : _nThreads(nThreads)
    { }
    
  /**
   * Iteratively call a function in parallel.
   * 
   * The function is expected to accept two size_t parameters, the loop
   * index and the thread id, e.g.:
   *   void loopFunction(size_t iteration, size_t threadID);
   * It is called (end-start) times.
   * 
   * This function is very similar to ThreadPool::For(), but does not
   * support recursion. For non-recursive loop, this function will be
   * faster.
   */
    template<typename Function>
    void Run(Iter start, Iter end, Function function)
    {
      _current = start;
      _end = end;
      std::vector<std::thread> threads;
      if(_nThreads > 1)
      {
        threads.reserve(_nThreads-1);
        for(unsigned t=1; t!=_nThreads; ++t)
          threads.emplace_back(&ParallelFor::run<Function>, this, t, function);
      }
      run<Function>(0, function);
      for(std::thread& thr : threads)
        thr.join();
    }
    
    size_t NThreads() const { return _nThreads; }
  private:
    template<typename Function>
    void run(size_t thread, Function function)
    {
      Iter iter;
      while(next(iter)) {
        function(iter, thread);
      }
    }
    
    bool next(Iter& iter)
    {
      std::lock_guard<std::mutex> lock(_mutex);
      if(_current == _end)
        return false;
      else {
        iter = _current;
        ++_current;
        return true;
      }
    }
    
    Iter _current, _end;
    std::mutex _mutex;
    size_t _nThreads;
  };
}

#endif
