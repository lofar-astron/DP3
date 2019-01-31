#ifndef BARRIER_H
#define BARRIER_H

#include <condition_variable>
#include <functional>
#include <mutex>

namespace DP3 {
	
	/**
	 * Simple barrier class.
	 * 
	 * This class is unfortunately necessary because boost::barrier had a 
	 * bug in completion functions, and std::barrier is still experimental.
	 * 
	 * It implements a barrier similar to boost::barrier.
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
		Barrier(size_t n, std::function<void()> completionFunction) : _n(n), _count(_n), _completionFunction(completionFunction)
		{
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
				_count = _n;
				_completionFunction();
				_condition.notify_all();
			}
			else {
				while(_count != _n)
					_condition.wait(lock);
			}
		}
		
	private:
		std::mutex _mutex;
		std::condition_variable _condition;
		size_t _n, _count;
		std::function<void()> _completionFunction;
	};
}

#endif
