#ifndef BUFFER_LANE_H
#define BUFFER_LANE_H

#include <vector>

#include "lane.h"

template<typename Tp>
class lane_write_buffer 
{
public:
	typedef typename ao::lane<Tp>::size_type size_type;
	typedef typename ao::lane<Tp>::value_type value_type;
	
	lane_write_buffer() : _buffer_size(0), _lane(0)
	{ }
	
	lane_write_buffer(ao::lane<Tp>* lane, size_type buffer_size) : _buffer_size(buffer_size), _lane(lane)
	{
		_buffer.reserve(buffer_size);
	}
	
	~lane_write_buffer()
	{
		flush();
	}
	
	void reset(ao::lane<Tp>* lane, size_type buffer_size)
	{
		_buffer.clear();
		_buffer.reserve(buffer_size);
		_buffer_size = buffer_size;
		_lane = lane;
	}
	
	void clear()
	{
		_lane->clear();
		_buffer.clear();
	}
		
	void write(const value_type& element)
	{
		_buffer.push_back(element);
		if(_buffer.size() == _buffer_size)
			flush();
	}
	
	void write(value_type&& element)
	{
		_buffer.push_back(std::move(element));
		if(_buffer.size() == _buffer_size)
			flush();
	}
	
	template<typename... Args>
	void emplace(Args&&... args)
	{
		_buffer.emplace_back(args...);
		if(_buffer.size() == _buffer_size)
			flush();
	}
	
	void write_end()
	{
		flush();
		_lane->write_end();
	}
		
	void flush()
	{
		_lane->move_write(&_buffer[0], _buffer.size());
		_buffer.clear();
	}
private:
	size_type _buffer_size;
	std::vector<value_type> _buffer;
	ao::lane<Tp>* _lane;
};

template<typename Tp>
class lane_read_buffer 
{
public:
	lane_read_buffer(ao::lane<Tp>* lane, size_t buffer_size) :
		_buffer(new Tp[buffer_size]),
		_buffer_size(buffer_size),
		_buffer_pos(0),
		_buffer_fill_count(0),
		_lane(lane)
	{
	}
	
	~lane_read_buffer()
	{
		delete[] _buffer;
	}
	
	bool read(Tp& element)
	{
		if(_buffer_pos == _buffer_fill_count)
		{
			_buffer_fill_count = _lane->read(_buffer, _buffer_size);
			_buffer_pos = 0;
			if(_buffer_fill_count == 0)
				return false;
		}
		element = std::move(_buffer[_buffer_pos]);
		++_buffer_pos;
		return true;
	}
	
private:
	lane_read_buffer(const lane_read_buffer&) = delete;
	lane_read_buffer& operator=(const lane_read_buffer&) = delete;

	Tp* _buffer;
	size_t _buffer_size, _buffer_pos, _buffer_fill_count;
	ao::lane<Tp>* _lane;
};

#endif
