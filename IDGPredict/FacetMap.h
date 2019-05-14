#ifndef FACET_MAP_H
#define FACET_MAP_H

#include <string>
#include <vector>
#include <cmath>

struct Vertex {
	Vertex(int _x, int _y) : x(_x), y(_y) { }
	int x, y;
};

class Facet
{
public:
  Facet() : _dirRA(0.0), _dirDec(0.0) { }
  
	typedef std::vector<Vertex>::iterator iterator;
	typedef std::vector<Vertex>::const_iterator const_iterator;
	
	iterator begin() { return _vertices.begin(); }
	iterator end() { return _vertices.end(); }
	const_iterator begin() const { return _vertices.begin(); }
	const_iterator end() const { return _vertices.end(); }
	
	void AddVertex(int x, int y)
	{
		_vertices.emplace_back(x, y);
	}
	
	void BoundingBox(int& x1, int& x2, int& y1, int& y2) const
	{
		if(_vertices.empty())
		{
			x1 = 0; x2 = 0;
			y1 = 0; y2 = 0;
		}
		else {
			x1 = _vertices.front().x; x2 = _vertices.front().x;
			y1 = _vertices.front().y; y2 = _vertices.front().y;
			for(auto i = _vertices.begin()+1; i!=_vertices.end(); ++i)
			{
				x1 = std::min(x1, i->x); x2 = std::max(x2, i->x);
				y1 = std::min(y1, i->y); y2 = std::max(y2, i->y);
			}
		}
	}
	
	bool HorizontalIntersections(int yIntersect, int& x1, int& x2) const
	{
		size_t nInts = 0;
		x1 = 0; x2 = 0;
		for(size_t i=0; i!=_vertices.size(); ++i)
		{
			Vertex
				v1 = _vertices[i],
				v2 = _vertices[(i+1)%_vertices.size()];
			if(v1.y > v2.y) std::swap(v1, v2);
			if(v1.y <= yIntersect && v2.y > yIntersect)
			{
				size_t x;
				if(v1.y == v2.y)
					x = std::min(v1.x, v2.x);
				else {
					double beta = double(v2.x - v1.x) / double(v2.y - v1.y);
					double xfl = v1.x + beta * (yIntersect - v1.y);
					x = round(xfl);
				}
				if(nInts == 0)
				{
					x1 = x;
					++nInts;
				}
				else {
					x2 = x;
					if(x1 > x2) std::swap(x1, x2);
					return true;
				}
			}
		}
		return false;
	}
	
	double RA() const { return _dirRA; }
	double Dec() const { return _dirDec; }
	
	void SetRA(double dirRA) { _dirRA = dirRA; }
	void SetDec(double dirDec) { _dirDec = dirDec; }
	
private:
	std::vector<Vertex> _vertices;
  double _dirRA, _dirDec;
};

class FacetMap
{
public:
	typedef std::vector<Facet>::const_iterator const_iterator;
	
	const_iterator begin() const { return _facets.begin(); }
	
	const_iterator end() const { return _facets.end(); }
	
	Facet& operator[](size_t index) { return _facets[index]; }
	
	Facet& AddFacet()
	{
		_facets.emplace_back();
		return _facets.back();
	}
	
	size_t NFacets() const { return _facets.size(); }
	
private:
	std::vector<Facet> _facets;
};

#endif
