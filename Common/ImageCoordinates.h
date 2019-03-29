#ifndef IMAGE_COORDINATES_H
#define IMAGE_COORDINATES_H

#include <algorithm>
#include <cmath>
#include <vector>

/**
 * This class collects all the LM coordinate transform as defined in
 * Perley (1999)'s "imaging with non-coplaner arrays".
 */
class ImageCoordinates
{
	public:
		template<typename T>
		static void RaDecToLM(T ra, T dec, T phaseCentreRa, T phaseCentreDec, T &destL, T &destM)
		{
			const T
				deltaAlpha = ra - phaseCentreRa,
				sinDeltaAlpha = sin(deltaAlpha),
				cosDeltaAlpha = cos(deltaAlpha),
				sinDec = sin(dec),
				cosDec = cos(dec),
				sinDec0 = sin(phaseCentreDec),
				cosDec0 = cos(phaseCentreDec);
			
			destL = cosDec * sinDeltaAlpha;
			destM = sinDec*cosDec0 - cosDec*sinDec0*cosDeltaAlpha;
		}
		
		template<typename T>
		static T RaDecToN(T ra, T dec, T phaseCentreRa, T phaseCentreDec)
		{
			const T
				cosDeltaAlpha = cos(ra - phaseCentreRa),
				sinDec = sin(dec),
				cosDec = cos(dec),
				sinDec0 = sin(phaseCentreDec),
				cosDec0 = cos(phaseCentreDec);
			
			return sinDec*sinDec0 + cosDec*cosDec0*cosDeltaAlpha;
		}
		
		template<typename T>
		static void LMToRaDec(T l, T m, T phaseCentreRa, T phaseCentreDec, T &destRa, T &destDec)
		{
			const T
				cosDec0 = cos(phaseCentreDec),
				sinDec0 = sin(phaseCentreDec),
				lmTerm = sqrt((T) 1.0 - l*l - m*m),
				deltaAlpha = atan2(l, lmTerm*cosDec0 - m*sinDec0);
				
			destRa = deltaAlpha + phaseCentreRa;
			destDec = asin(m*cosDec0 + lmTerm*sinDec0);
		}
		
		template<typename T>
		static void XYToLM(size_t x, size_t y, T pixelSizeX, T pixelSizeY, size_t width, size_t height, T &l, T &m)
		{
			T midX = (T) width / 2.0, midY = (T) height / 2.0;
			l = (midX - (T) x) * pixelSizeX;
			m = ((T) y - midY) * pixelSizeY;
		}
		
		template<typename T>
		static void LMToXY(T l, T m, T pixelSizeX, T pixelSizeY, size_t width, size_t height, int &x, int &y)
		{
			T midX = (T) width / 2.0, midY = (T) height / 2.0;
			x = round(-l / pixelSizeX) + midX;
			y = round(m / pixelSizeY) + midY;
		}
		
		template<typename T>
		static void LMToXYfloat(T l, T m, T pixelSizeX, T pixelSizeY, size_t width, size_t height, T &x, T &y)
		{
			T midX = (T) width / 2.0, midY = (T) height / 2.0;
			x = -l / pixelSizeX + midX;
			y = m / pixelSizeY + midY;
		}
		
		template<typename T>
		static T AngularDistance(T ra1, T dec1, T ra2, T dec2)
		{
			T sinDec1, sinDec2, cosDec1, cosDec2;
			SinCos(dec1, &sinDec1, &cosDec1);
			SinCos(dec2, &sinDec2, &cosDec2);
			T cosVal = sinDec1*sinDec2 + cosDec1*cosDec2*std::cos(ra1 - ra2);
			// Rounding errors sometimes cause cosVal to be slightly larger than 1, which would cause
			// an NaN return value.
			return cosVal <= 1.0 ? std::acos(cosVal) : 0.0;
		}
		
		template<typename T>
		static T MeanRA(const std::vector<T>& raValues)
		{
			std::vector<T> sorted(raValues);
			for(size_t i=0; i!=sorted.size(); ++i) {
				while(sorted[i] >= 2*M_PI) sorted[i] -= 2.0*M_PI;
				while(sorted[i] < 0.0) sorted[i] += 2.0*M_PI;
			}
			std::sort(sorted.begin(), sorted.end());
			T gapSize = 0.0, gapCentre = 0.0;
			for(size_t i=0; i!=sorted.size(); ++i)
			{
				double dist;
				if(i == sorted.size()-1)
					dist = 2.0*M_PI + sorted.front() - sorted.back();
				else
					dist = sorted[i+1] - sorted[i];
				if(dist > gapSize)
				{
					gapSize = dist;
					gapCentre = sorted[i] + gapSize*0.5;
				}
			}
			if(gapCentre >= 2.0*M_PI) gapCentre-=2.0*M_PI;
			T sum = 0.0;
			for(size_t i=0; i!=sorted.size(); ++i)
			{
				if(sorted[i] < gapCentre)
					sum += sorted[i];
				else
					sum += sorted[i] - 2.0*M_PI;
			}
			sum /= sorted.size();
			if(sum < 0.0)
				return sum + 2.0*M_PI;
			else
				return sum;
		}
	private:
		static void SinCos(double angle, double* sinAngle, double* cosAngle)
		{
			*sinAngle = std::sin(angle);
			*cosAngle = std::cos(angle);
		}
		
		static void SinCos(long double angle, long double* sinAngle, long double* cosAngle)
		{
			*sinAngle = std::sin(angle);
			*cosAngle = std::cos(angle);
		}
		
		static void SinCos(float angle, float* sinAngle, float* cosAngle)
		{
			*sinAngle = std::sin(angle);
			*cosAngle = std::cos(angle);
		}
		
		ImageCoordinates();
};

#endif
