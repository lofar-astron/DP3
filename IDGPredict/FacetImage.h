#ifndef FACET_IMAGE_H
#define FACET_IMAGE_H

#include "FacetMap.h"

#include "../Common/UVector.h"

class FacetImage
{
public:
	FacetImage() :
		_data(),
		_offsetX(0),
		_offsetY(0)
	{ }
	
	size_t Width() const { return _width; }
	
	size_t Height() const { return _height; }
	
	int OffsetX() const { return _offsetX; }
	
	int OffsetY() const { return _offsetY; }
	
	void Set(size_t width, size_t height, double value)
	{
		_width = width;
		_height = height;
		_data.assign(width * height, value);
	}
	
	void CopyFacetPart(const Facet& facet, const double* input, size_t inputWidth, size_t inputHeight, double padding, bool makeSquare)
	{
		Facet cFacet = clippedFacet(facet, inputWidth, inputHeight);
		int x1, x2, y1, y2;
		cFacet.BoundingBox(x1, x2, y1, y2);
		size_t width = x2-x1, height = y2-y1;
		size_t paddedWidth = (size_t) ceil(width * padding);
		size_t paddedHeight = (size_t) ceil(height * padding);
		if(makeSquare)
		{
			paddedWidth = std::max(paddedWidth, paddedHeight);
			paddedHeight = paddedWidth;
		}
		// Make the width and height divisable by four.
		paddedWidth += (4-(paddedWidth%4))%4;
		paddedHeight += (4-(paddedHeight%4))%4;
		int padX = (paddedWidth-width)/2;
		int padY = (paddedHeight-height)/2;
    Set(paddedWidth, paddedHeight, 0.0);
		_offsetX = x1-padX;
		_offsetY = y1-padY;
		for(int y=y1; y!=y2; ++y)
		{
			int xi1, xi2;
			if(cFacet.HorizontalIntersections(y, xi1, xi2))
			{
				if(xi1<0) xi1=0;
				if(xi2<0) xi2=0;
				if(xi2-xi1>int(width)) xi2=xi1+width;
				for(int x=xi1; x!=xi2; ++x)
				{
					_data[(y-y1+padY)*paddedWidth + x-x1+padX] = input[y*inputWidth + x];
				}
			}
		}
	}
	
	void FillFacet(const Facet& facet, int colour)
	{
		Facet cFacet = clippedFacet(facet, _width, _height);
		int x1, x2, y1, y2;
		cFacet.BoundingBox(x1, x2, y1, y2);
		size_t width = x2-x1;
		for(int y=y1; y!=y2; ++y)
		{
			int xi1, xi2;
			if(cFacet.HorizontalIntersections(y, xi1, xi2))
			{
				if(xi1<0) xi1=0;
				if(xi2<0) xi2=0;
				if(xi2-xi1>int(width)) xi2=xi1+width;
				for(int x=xi1; x!=xi2; ++x)
				{
					_data[y*_width + x] = colour;
				}
			}
		}
	}
	
	double* Data() { return _data.data(); }
	
	ao::uvector<double> AcquireData() { return std::move(_data); }
	
private:
	Facet clippedFacet(const Facet& input, int width, int height)
	{
		Facet clippedFacet(input);
		
		for(auto& v : clippedFacet)
		{
			if(v.x < 0) v.x = 0;
			if(v.y < 0) v.y = 0;
			if(v.x > width) v.x = width;
			if(v.y > height) v.y = height;
		}
		return clippedFacet;
	}
	
	ao::uvector<double> _data;
	size_t _width, _height;
	int _offsetX, _offsetY;
};

#endif

