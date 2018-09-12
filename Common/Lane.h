#ifndef AO_LANE_11_H
#define AO_LANE_11_H

#include <cstring>
#include <deque>
#include <mutex>
#include <condition_variable>

/**
 * @file
 * Internal header file for the lane.
 * @headername{lane.h}
 */

//#define LANE_DEBUG_MODE

#ifdef LANE_DEBUG_MODE
#include <string>
#include <iostream>
#include <sstream>
#include <cmath>
#endif

namespace ao
{

#ifdef LANE_DEBUG_MODE
#define set_lane_debug_name(lane, str) (lane).setDebugName(str)
#define LANE_REGISTER_DEBUG_INFO registerDebugInfo()
#define LANE_REGISTER_DEBUG_WRITE_WAIT registerDebugWriteWait()
#define LANE_REGISTER_DEBUG_READ_WAIT registerDebugReadWait()
#define LANE_REPORT_DEBUG_INFO reportDebugInfo()
	
#else

#define set_lane_debug_name(lane, str)
#define LANE_REGISTER_DEBUG_INFO
#define LANE_REGISTER_DEBUG_WRITE_WAIT
#define LANE_REGISTER_DEBUG_READ_WAIT
#define LANE_REPORT_DEBUG_INFO

#endif

/**
 * @brief The lane is an efficient cyclic buffer that is synchronized.
 * @details
 * A lane can typically be used in a multi-threaded producer-consumer
 * situation. The lane also holds a state which allows for
 * an ellegant way of communicating from producer(s) to
 * consumer(s) that all data has been produced.
 * 
 * A simple example:
 * @code
 * void producer(lane<Task>* taskLane)
 * {
 *   while(moreTasks)
 *     taskLane->write(nextTask());
 *   taskLane->write_end();
 * }
 * 
 * void consumer(lane<Task>* taskLane)
 * {
 *   Task task;
 *   while(taskLane->read(task))
 *     processTask(task);
 * }
 * 
 * void run()
 * {
 *   lane<Task> taskLane;
 *   std::thread consumerThread(&consumer(), &taskLane);
 *   producer(&taskLane);
 *   consumerThread.join();
 * }
 * @endcode
 * 
 * The various read and write methods, as well as the empty(),
 * capacity() and size() methods are always thread safe. The other
 * methods are not: assignment, swap(), clear() and resize() can not
 * be called from a different thread while another thread is also
 * accessing the lane. The same holds obviously for the constructors
 * and destructor. This is chosen because these methods should almost never
 * be called in parallel with other methods, and hence it is not worth
 * to increase every call with extra locks to make this possible.
 * 
 * With one reader and one writer, the order is guaranteed to be consistent.
 * With multiple readers or writers in combination with multi-element
 * write or read functions, a sequence of symbols might be interrupted. For
 * example, if a multi-element write() won't fit completely in the buffer,
 * the thread will wait for free space. Another thread might get now write
 * access first, causing the single call to the multi-element write to be
 * "split up".
 * 
 * @author Andre Offringa
 * @tparam Tp Type of elements to be stored in the lane.
 */
template<typename Tp>
class lane
{
	public:
		/** @brief Integer type used to store size types. */
		typedef std::size_t size_type;
		
		/** @brief Type of elements stored in the lane. */
		typedef Tp value_type;
		
		/** @brief Construct a lane with zero elements.
		 * @details A lane with zero elements can not be written to or read to
		 * (both operations will wait forever).
		 * 
		 * This constructor makes it easy to construct e.g. a container
		 * of lanes. After the container is created, the lanes can be
		 * resized with @ref resize().
		 */
		lane() noexcept :
			_buffer(0),
			_capacity(0),
			_write_position(0),
			_free_write_space(0),
			_status(status_normal)
		{
		}
		
		/** @brief Construct a lane with the given capacity.
		 * @details After construction, the lane is ready for writing to and reading from.
		 * @param capacity Number of elements that the lane can hold at once.
		 */
		explicit lane(size_t capacity) :
			_buffer(new Tp[capacity]),
			_capacity(capacity),
			_write_position(0),
			_free_write_space(_capacity),
			_status(status_normal)
		{
		}
		
		lane(const lane<Tp>& source) = delete;
		
		/** @brief Move construct a lane.
		 * @details This operation is not thread safe: the behaviour is undefined when
		 * other threads access the source lane.
		 * @param source Original lane to be moved from.
		 */
		lane(lane<Tp>&& source) noexcept :
			_buffer(0),
			_capacity(0),
			_write_position(0),
			_free_write_space(0),
			_status(status_normal)
		{
			swap(source);
		}
		
		/** @brief Destructor.
		 * @details The destructor is not synchronized.
		 */
		~lane()
		{
			LANE_REPORT_DEBUG_INFO;
			delete[] _buffer;
		}
		
		lane<Tp>& operator=(const lane<Tp>& source) = delete;
		
		/** @brief Move assignment.
		 * @details This operation is not thread safe: the behaviour is undefined when
		 * other threads access the source lane.
		 * @param source Original lane to be moved from.
		 * @returns This lane.
		 */
		lane<Tp>& operator=(lane<Tp>&& source) noexcept
		{
			swap(source);
			return *this;
		}
		
		/** @brief Swap the contents of this lane with another.
		 * @details This operation is not thread safe: the behaviour is undefined when
		 * other threads access either lane.
		 */
		void swap(lane<Tp>& other) noexcept
		{
			std::swap(_buffer, other._buffer);
			std::swap(_capacity, other._capacity);
			std::swap(_write_position, other._write_position);
			std::swap(_free_write_space, other._free_write_space);
			std::swap(_status, other._status);
		}
		
		/** @brief Clear the contents and reset the state of the lane.
		 * @details After calling clear(), the lane is in the same state as after
		 * construction. This also means that after clearing the lane, it
		 * is as if write_end() has not been called yet.
		 * 
		 * This method is not thread safe.
		 */
		void clear() noexcept
		{
			_write_position = 0;
			_free_write_space = _capacity;
			_status = status_normal;
		}
		
		/** @brief Write a single element.
		 * @details This method is thread safe, and can be called together with
		 * other write and read methods from different threads.
		 * 
		 * If this call comes after a call to write_end(), the call
		 * will be ignored.
		 * @param element Object to be copied into the cyclic buffer.
		 */
		void write(const value_type& element)
		{
			std::unique_lock<std::mutex> lock(_mutex);
			LANE_REGISTER_DEBUG_INFO;
			
			if(_status == status_normal)
			{
				while(_free_write_space == 0)
				{
					LANE_REGISTER_DEBUG_WRITE_WAIT;
					_writing_possible_condition.wait(lock);
				}
				
				_buffer[_write_position] = element;
				_write_position = (_write_position+1) % _capacity;
				--_free_write_space;
				// Now that there is less free write space, there is more free read
				// space and thus readers can possibly continue.
				_reading_possible_condition.notify_all();
			}
		}

		/** @brief Write a single element by constructing it.
		 * @details This method is thread safe, and can be called together with
		 * other write and read methods from different threads.
		 * 
		 * If this call comes after a call to write_end(), the call
		 * will be ignored. The implementation does not construct the value
		 * in place, but rather constructs the value and then move assigns it.
		 * This is because the value that it is moved into has already been
		 * constructed (in the current implementation).
		 * @param element Object to be moved into the cyclic buffer.
		 */
		template<typename... Args>
		void emplace(Args&&... args)
		{
			std::unique_lock<std::mutex> lock(_mutex);
			LANE_REGISTER_DEBUG_INFO;

			if(_status == status_normal)
			{
				while(_free_write_space == 0)
				{
					LANE_REGISTER_DEBUG_WRITE_WAIT;
					_writing_possible_condition.wait(lock);
				}

				_buffer[_write_position] = value_type(args...);
				_write_position = (_write_position+1) % _capacity;
				--_free_write_space;
				// Now that there is less free write space, there is more free read
				// space and thus readers can possibly continue.
				_reading_possible_condition.notify_all();
			}
		}
		
		/** @brief Write a single element by moving it in.
		 * @details This method is thread safe, and can be called together with
		 * other write and read methods from different threads.
		 * 
		 * If this call comes after a call to write_end(), the call
		 * will be ignored.
		 * @param element Object to be moved into the cyclic buffer.
		 */
		void write(value_type&& element)
		{
			std::unique_lock<std::mutex> lock(_mutex);
			LANE_REGISTER_DEBUG_INFO;
			
			if(_status == status_normal)
			{
				while(_free_write_space == 0)
				{
					LANE_REGISTER_DEBUG_WRITE_WAIT;
					_writing_possible_condition.wait(lock);
				}
				
				_buffer[_write_position] = std::move(element);
				_write_position = (_write_position+1) % _capacity;
				--_free_write_space;
				// Now that there is less free write space, there is more free read
				// space and thus readers can possibly continue.
				_reading_possible_condition.notify_all();
			}
		}
		
		void write(const value_type* elements, size_t n)
		{
			write_generic(elements, n);
		}
		
		void move_write(value_type* elements, size_t n)
		{
			write_generic(elements, n);
		}
		
		bool read(value_type& destination)
		{
			std::unique_lock<std::mutex> lock(_mutex);
			LANE_REGISTER_DEBUG_INFO;
			while(free_read_space() == 0 && _status == status_normal)
			{
				LANE_REGISTER_DEBUG_READ_WAIT;
				_reading_possible_condition.wait(lock);
			}
			if(free_read_space() == 0)
				return false;
			else
			{
				destination = std::move(_buffer[read_position()]);
				++_free_write_space;
				// Now that there is more free write space, writers can possibly continue.
				_writing_possible_condition.notify_all();
				return true;
			}
		}
		
		size_t read(value_type* destinations, size_t n)
		{
			size_t n_left = n;
			
			std::unique_lock<std::mutex> lock(_mutex);
			LANE_REGISTER_DEBUG_INFO;
			
			size_t free_space = free_read_space();
			size_t read_size = free_space > n ? n : free_space;
			immediate_read(destinations, read_size);
			n_left -= read_size;
			
			while(n_left != 0 && _status == status_normal)
			{
				destinations += read_size;
				
				do {
					LANE_REGISTER_DEBUG_READ_WAIT;
					_reading_possible_condition.wait(lock);
				} while(free_read_space() == 0 && _status == status_normal);
				
				free_space = free_read_space();
				read_size = free_space > n_left ? n_left : free_space;
				immediate_read(destinations, read_size);
				n_left -= read_size;
			}
			return n - n_left;
		}
		
		void write_end()
		{
			std::lock_guard<std::mutex> lock(_mutex);
			LANE_REGISTER_DEBUG_INFO;
			_status = status_end;
			_writing_possible_condition.notify_all();
			_reading_possible_condition.notify_all();
		}
		
		size_t capacity() const noexcept
		{
			return _capacity;
		}
		
		size_t size() const
		{
			std::lock_guard<std::mutex> lock(_mutex);
			return _capacity - _free_write_space;
		}
		
		bool empty() const
		{
			std::lock_guard<std::mutex> lock(_mutex);
			return _capacity == _free_write_space;
		}
		
		/**
		 * Change the capacity of the lane. This will erase all data in the lane.
		 */
		void resize(size_t new_capacity)
		{
			Tp *new_buffer = new Tp[new_capacity];
			delete[] _buffer;
			_buffer = new_buffer;
			_capacity = new_capacity;
			_write_position = 0;
			_free_write_space = new_capacity;
			_status = status_normal;
		}
		
		/**
		 * Wait until this lane is empty.
		 */
		void wait_for_empty()
		{
			std::unique_lock<std::mutex> lock(_mutex);
			while(_capacity != _free_write_space)
			{
				_reading_possible_condition.wait(lock);
			}
		}
		
#ifdef LANE_DEBUG_MODE
		/**
		 * Change the name of this lane to make it appear in the output along
		 * with statistics. Do not use this function directly; use the
		 * set_lane_debug_name() macro instead.
		 * @param nameStr New debug description of this lane.
		 */
		void setDebugName(const std::string& nameStr)
		{
			_debugName = nameStr;
		}
#endif
	private:
		Tp* _buffer;
		
		size_t _capacity;
		
		size_t _write_position;
		
		size_t _free_write_space;
		
		enum { status_normal, status_end } _status;
		
		mutable std::mutex _mutex;
		
		std::condition_variable _writing_possible_condition, _reading_possible_condition;
		
		size_t read_position() const noexcept
		{
			return (_write_position + _free_write_space) % _capacity;
		}
		
		size_t free_read_space() const noexcept
		{
			return _capacity - _free_write_space;
		}
		
		// This is a template to allow const and non-const (to be able to move)
		template<typename T>
		void write_generic(T* elements, size_t n)
		{
			std::unique_lock<std::mutex> lock(_mutex);
			LANE_REGISTER_DEBUG_INFO;
			
			if(_status == status_normal)
			{
				size_t write_size = _free_write_space > n ? n : _free_write_space;
				immediate_write(elements, write_size);
				n -= write_size;
				
				while(n != 0) {
					elements += write_size;
				
					do {
						LANE_REGISTER_DEBUG_WRITE_WAIT;
						_writing_possible_condition.wait(lock);
					} while(_free_write_space == 0 && _status == status_normal);
					
					write_size = _free_write_space > n ? n : _free_write_space;
					immediate_write(elements, write_size);
					n -= write_size;
				} while(n != 0);
			}
		}
		
		// This is a template to allow const and non-const (to be able to move)
		template<typename T>
		void immediate_write(T *elements, size_t n) noexcept
		{
			// Split the writing in two ranges if needed. The first range fits in
			// [_write_position, _capacity), the second range in [0, end). By doing
			// so, we only have to calculate the modulo in the write position once.
			if(n > 0)
			{
				size_t nPart;
				if(_write_position + n > _capacity)
				{
					nPart = _capacity - _write_position;
				} else {
					nPart = n;
				}
				for(size_t i = 0; i < nPart ; ++i, ++_write_position)
				{
					_buffer[_write_position] = std::move(elements[i]);
				}
				
				_write_position = _write_position % _capacity;
				
				for(size_t i = nPart; i < n ; ++i, ++_write_position)
				{
					_buffer[_write_position] = std::move(elements[i]);
				}
				
				_free_write_space -= n;
				
				// Now that there is less free write space, there is more free read
				// space and thus readers can possibly continue.
				_reading_possible_condition.notify_all();
			}
		}
		
		void immediate_read(value_type *elements, size_t n) noexcept
		{
			// As with write, split in two ranges if needed. The first range fits in
			// [read_position(), _capacity), the second range in [0, end).
			if(n > 0)
			{
				size_t nPart;
				size_t position = read_position();
				if(position + n > _capacity)
				{
					nPart = _capacity - position;
				} else {
					nPart = n;
				}
				for(size_t i = 0; i < nPart ; ++i, ++position)
				{
					elements[i] = std::move(_buffer[position]);
				}
				
				position = position % _capacity;
				
				for(size_t i = nPart; i < n ; ++i, ++position)
				{
					elements[i] = std::move(_buffer[position]);
				}
				
				_free_write_space += n;
				
				// Now that there is more free write space, writers can possibly continue.
				_writing_possible_condition.notify_all();
			}
		}
#ifdef LANE_DEBUG_MODE
		void registerDebugInfo() noexcept
		{
			_debugSummedSize += _capacity - _free_write_space;
			_debugMeasureCount++;
		}
		void registerDebugReadWait() noexcept
		{
			++_debugReadWaitCount;
		}
		void registerDebugWriteWait() noexcept
		{
			++_debugWriteWaitCount;
		}
		void reportDebugInfo()
		{
			if(!_debugName.empty())
			{
				std::stringstream str;
				str
					<< "*** Debug report for the following lane: ***\n"
					<< "\"" << _debugName << "\"\n"
					<< "Capacity: " << _capacity << '\n'
					<< "Total read/write ops: " << _debugMeasureCount << '\n'
					<< "Average size of buffer, measured per read/write op.: " << round(double(_debugSummedSize)*100.0/_debugMeasureCount)/100.0 << '\n'
					<< "Number of wait events during reading: " << _debugReadWaitCount << '\n'
					<< "Number of wait events during writing: " << _debugWriteWaitCount << '\n';
				std::cout << str.str();
			}
		}
		std::string _debugName;
		size_t
			_debugSummedSize = 0, _debugMeasureCount = 0,
			_debugReadWaitCount = 0, _debugWriteWaitCount = 0;
#endif
};

template<typename Tp>
void swap(ao::lane<Tp>& first, ao::lane<Tp>& second) noexcept
{
	first.swap(second);
}

} // end of namespace

#endif // AO_LANE11_H
