// Copyright (C) 2002-2003
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

#ifndef THREAD_POOL_H
#define THREAD_POOL_H

#include <sched.h>

#include <condition_variable>
#include <functional>
#include <map>
#include <mutex>
#include <thread>
#include <vector>

namespace DP3
{

/// \brief Defines a class for managing threads.
class ThreadPool
{
public:
	/**
	 * Returns the number of CPUs in this machine.
	 */
	static size_t NCPUs()
	{
#ifdef __APPLE__
		return sysconf(_SC_NPROCESSORS_ONLN);
#else
		cpu_set_t cs;
		CPU_ZERO(&cs);
		sched_getaffinity(0, sizeof cs , &cs);
		return CPU_COUNT(&cs);
#endif
	}

	/**
	 * Create a thread pool with NThreads()==NCPUs().
	 */
	ThreadPool() :
		_isStopped(false),
		_priority(0)
	{
		size_t nThreads = NCPUs();
		if(nThreads == 0)
			nThreads = 1;
		// We reserve one thread less, because we always want a new For loop
		// to be able to add a new thread (with index 0).
		_threads.reserve(nThreads-1);
		for(size_t i=1; i!=nThreads; ++i)
			_threads.emplace_back(&ThreadPool::threadFunc, this, i);
	}

	/**
	 * Create a thread pool with the specified number of threads.
	 */
	ThreadPool(size_t nThreads) :
		_isStopped(false),
		_priority(0)
	{
		if(nThreads == 0)
			throw std::runtime_error("A threadPool was created with nThreads=0");
		// We reserve one thread less, because we always want a new For loop
		// to be able to add a new thread (with index 0).
		_threads.reserve(nThreads-1);
		for(size_t i=1; i!=nThreads; ++i)
			_threads.emplace_back(&ThreadPool::threadFunc, this, i);
	}

	ThreadPool(const ThreadPool&) = delete;
	ThreadPool& operator=(const ThreadPool&) = delete;

	~ThreadPool()
	{
		std::unique_lock<std::mutex> lock(_mutex);
		_isStopped = true;
		_onProgress.notify_all();
		lock.unlock();
		for(std::thread& t : _threads)
			t.join();
	}

	size_t NThreads() const
	{
		return _threads.size()+1;
	}

  void SetNThreads(size_t nThreads)
  {
    if(nThreads != NThreads())
    {
      std::unique_lock<std::mutex> lock(_mutex);
      _isStopped = true;
      _onProgress.notify_all();
      lock.unlock();
      for(std::thread& t : _threads)
        t.join();

      _threads.clear();

      _threads.reserve(nThreads-1);
      for(size_t i=1; i!=nThreads; ++i)
        _threads.emplace_back(&ThreadPool::threadFunc, this, i);
    }
  }

	/**
	 * Iteratively call a function in parallel.
	 *
	 * The function is expected to accept two size_t parameters, the loop
	 * index and the thread id, e.g.:
	 *   void loopFunction(size_t iteration, size_t threadID);
	 * It is called (end-start) times.
	 */
	template<typename Func>
	void For(size_t start, size_t end, Func func)
	{
		std::unique_lock<std::mutex> lock(_mutex);
		size_t thisPriority = _priority;
		++_priority;
		lock.unlock();

		size_t progress = end-start;

		std::thread localThread(&ThreadPool::threadSpecificPriorityFunc, this, 0, thisPriority, &progress);

		// Queue tasks for all iterations
		while(start!=end)
		{
			write(thisPriority, std::bind(func, start, std::placeholders::_1), &progress);
			++start;
		}

		localThread.join();
	}

private:
	void threadFunc(size_t threadId)
	{
		std::pair<std::function<void(size_t)>, size_t*> func;
		while(read_highest_priority(func))
		{
			func.first(threadId);

			std::unique_lock<std::mutex> lock(_mutex);
			--(*func.second); // decrease progress counter (requires lock)
			_onProgress.notify_all();
		}
	}

	void threadSpecificPriorityFunc(size_t threadId, size_t priority, size_t* progressPtr)
	{
		std::pair<std::function<void(size_t)>, size_t*> func;
		while(read_specific_priority(priority, func, progressPtr))
		{
			func.first(threadId);

			std::unique_lock<std::mutex> lock(_mutex);
			--(*progressPtr);
			_onProgress.notify_all();
		}
	}

	bool read_highest_priority(std::pair<std::function<void(size_t)>, size_t*>& func)
	{
		std::unique_lock<std::mutex> lock(_mutex);
		while(!_isStopped && _tasks.empty())
			_onProgress.wait(lock);
		if(!_tasks.empty())
		{
			func = std::move(_tasks.begin()->second);
			_tasks.erase(_tasks.begin());
			_onProgress.notify_all();
			return true;
		}
		else {
			return false;
		}
	}

	bool read_specific_priority(size_t priority, std::pair<std::function<void(size_t)>, size_t*>& func, size_t* progress)
	{
		std::unique_lock<std::mutex> lock(_mutex);
		auto iter = _tasks.find(priority);
		while(!_isStopped && (*progress)>0 && iter == _tasks.end())
		{
			_onProgress.wait(lock);
			iter = _tasks.find(priority);
		}
		if(iter != _tasks.end())
		{
			func = std::move(iter->second);
			_tasks.erase(iter);
			_onProgress.notify_all();
			return true;
		}
		else {
			return false;
		}
	}

	void write(size_t priority, std::function<void(size_t)>&& func, size_t* progressPtr)
	{
		// Wait until there is space in the map (so that the map
		// doesn't get too large)
		std::unique_lock<std::mutex> lock(_mutex);
		while(_tasks.count(priority) >= NThreads())
		{
			_onProgress.wait(lock);
		}
		_tasks.emplace(priority, std::make_pair(std::move(func), progressPtr));
		_onProgress.notify_all();
	}

	// Priority, (function, progress*)
	bool _isStopped;
	size_t _priority;
	std::multimap<
		size_t,
		std::pair<std::function<void(size_t)>, size_t*>,
		std::greater<size_t>
	> _tasks;
	std::vector<std::thread> _threads;
	std::mutex _mutex;
	std::condition_variable _onProgress;
};

};

#endif
