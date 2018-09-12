#ifndef THREAD_POOL_H
#define THREAD_POOL_H

#include <sched.h>

#include <condition_variable>
#include <functional>
#include <map>
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
		_freeThreads = nthreads;
		_tasks.resize(nthreads);
		_threads.reserve(nthreads);
		for(size_t i=0; i!=nthreads; ++i)
			_threads.emplace_back(&ThreadPool::threadFunc, this, i);
	}
	
	ThreadPool(const ThreadPool&) = delete;
	ThreadPool& operator=(const ThreadPool&) = delete;
	
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
	  size_t progress = 0, n = end-start;
		
		std::unique_lock<std::mutex> lock(_mutex);
		if(_freeThreads == 0)
		{
			lock.unlock();
			// All our threads are busy
			while(start!=end)
			{
				func(start, 0);
				++start;
			}
		}
		else {
			lock.unlock();
			// Queue tasks for all iterations
			while(start!=end)
			{
				_tasks.emplace(std::bind(func, start, std::placeholders::_1), &progress);
				++start;
			}
			
			// Wait untill we have performed all iterations
			std::unique_lock<std::mutex> lock(_mutex);
			while(progress != n)
			{
				_onProgress.wait(lock);
			}
		}
	}
	
private:
	void threadFunc(size_t threadId)
	{
		std::pair<std::function<void(size_t)>, size_t*> func;
		while(_tasks.read(func))
		{
			std::unique_lock<std::mutex> lock(_mutex);
			_freeThreads--;
			lock.unlock();
			
			func.first(threadId);
			
			lock.lock();
			++_freeThreads;
			++(*func.second);
			_onProgress.notify_all();
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
	
	ao::lane<std::pair<std::function<void(size_t)>, size_t*>> _tasks;
	std::vector<std::thread> _threads;
	std::mutex _mutex;
	std::condition_variable _onProgress;
	size_t _freeThreads;
};

};

#endif
