#ifndef THREAD_POOL_H
#define THREAD_POOL_H

#include <sched.h>

#include <condition_variable>
#include <functional>
#include <mutex>
#include <thread>
#include <vector>

#include "Lane.h"

namespace DP3
{

class ThreadPool
{
public:
	ThreadPool()
	{
		size_t nthreads = cpus();
		_threads.reserve(nthreads);
		for(size_t i=0; i!=nthreads; ++i)
			_threads.emplace_back(&ThreadPool::threadFunc, this, i);
		_tasks.resize(nthreads);
	}
	
	~ThreadPool()
	{
		_tasks.write_end();
		for(std::thread& t : _threads)
			t.join();
	}
	
	size_t NThreads() const
	{
		return _threads.size();
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
		while(start!=end)
		{
			_tasks.emplace(std::bind(func, start, std::placeholders::_1));
			++start;
		}
	}
	
private:
	void threadFunc(size_t threadId)
	{
		std::function<void(size_t)> func;
		while(_tasks.read(func))
		{
			func(threadId);
		}
	}
	
	static unsigned cpus()
	{
#ifdef __APPLE__
		return sysconf(_SC_NPROCESSORS_ONLN);
#else
		cpu_set_t cs;
		CPU_ZERO(&cs);
		sched_getaffinity(0, sizeof cs , &cs);

		int count = 0;
		for (int i = 0; i < CPU_SETSIZE; i++)
		{
			if (CPU_ISSET(i, &cs))
				++count;
		}
		return count;
#endif
	}
	
	ao::lane<std::function<void(size_t)>> _tasks;
	std::vector<std::thread> _threads;
};

};

#endif
