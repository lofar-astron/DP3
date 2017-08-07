#ifndef MATRIX_2X2_H
#define MATRIX_2X2_H

#include <complex>
#include <limits>
#include <sstream>

class Matrix2x2
{
public:
	template<typename LHS_T, typename RHS_T>
	static void Assign(std::complex<LHS_T>* dest, const std::complex<RHS_T>* source)
	{
		for(size_t p=0; p!=4; ++p)
			dest[p] = source[p];
	}
	
	template<typename LHS_T, typename RHS_T>
	static void Assign(LHS_T* dest, const RHS_T* source)
	{
		for(size_t p=0; p!=4; ++p)
			dest[p] = source[p];
	}
	
	template<typename T, typename RHS_T>
	static void Add(std::complex<T>* dest, const RHS_T * rhs)
	{
		for(size_t p=0; p!=4; ++p)
			dest[p] += rhs[p];
	}
	
	template<typename T>
	static void Subtract(std::complex<T>* dest, const std::complex<T>* rhs)
	{
		for(size_t p=0; p!=4; ++p)
			dest[p] -= rhs[p];
	}
	
	template<typename T>
	static bool IsFinite(const std::complex<T>* matrix)
	{
		return
			std::isfinite(matrix[0].real()) && std::isfinite(matrix[0].imag()) &&
			std::isfinite(matrix[1].real()) && std::isfinite(matrix[1].imag()) &&
			std::isfinite(matrix[2].real()) && std::isfinite(matrix[2].imag()) &&
			std::isfinite(matrix[3].real()) && std::isfinite(matrix[3].imag());
	}
	
	template<typename T>
	static void ScalarMultiply(std::complex<T>* dest, T factor)
	{
		for(size_t p=0; p!=4; ++p)
			dest[p] *= factor;
	}
	
	template<typename T>
	static void ScalarMultiply(T* dest, T factor)
	{
		for(size_t p=0; p!=4; ++p)
			dest[p] *= factor;
	}
	
	template<typename T, typename RHS>
	static void MultiplyAdd(std::complex<T>* dest, const RHS* rhs, T factor)
	{
		for(size_t p=0; p!=4; ++p)
			dest[p] += rhs[p] * factor;
	}
	
	template<typename ComplType, typename LHS_T, typename RHS_T>
	static void ATimesB(std::complex<ComplType>* dest, const LHS_T* lhs, const RHS_T* rhs)
	{
		dest[0] = lhs[0] * rhs[0] + lhs[1] * rhs[2];
		dest[1] = lhs[0] * rhs[1] + lhs[1] * rhs[3];
		dest[2] = lhs[2] * rhs[0] + lhs[3] * rhs[2];
		dest[3] = lhs[2] * rhs[1] + lhs[3] * rhs[3];
	}
	
	static void PlusATimesB(std::complex<double> *dest, const std::complex<double> *lhs, const std::complex<double> *rhs)
	{
		dest[0] += lhs[0] * rhs[0] + lhs[1] * rhs[2];
		dest[1] += lhs[0] * rhs[1] + lhs[1] * rhs[3];
		dest[2] += lhs[2] * rhs[0] + lhs[3] * rhs[2];
		dest[3] += lhs[2] * rhs[1] + lhs[3] * rhs[3];
	}
	
	template<typename ComplType, typename LHS_T, typename RHS_T>
	static void ATimesHermB(std::complex<ComplType> *dest, const LHS_T* lhs, const RHS_T* rhs)
	{
		dest[0] = lhs[0] * std::conj(rhs[0]) + lhs[1] * std::conj(rhs[1]);
		dest[1] = lhs[0] * std::conj(rhs[2]) + lhs[1] * std::conj(rhs[3]);
		dest[2] = lhs[2] * std::conj(rhs[0]) + lhs[3] * std::conj(rhs[1]);
		dest[3] = lhs[2] * std::conj(rhs[2]) + lhs[3] * std::conj(rhs[3]);
	}

	template<typename ComplType, typename LHS_T, typename RHS_T>
	static void PlusATimesHermB(std::complex<ComplType> *dest, const LHS_T* lhs, const RHS_T* rhs)
	{
		dest[0] += lhs[0] * std::conj(rhs[0]) + lhs[1] * std::conj(rhs[1]);
		dest[1] += lhs[0] * std::conj(rhs[2]) + lhs[1] * std::conj(rhs[3]);
		dest[2] += lhs[2] * std::conj(rhs[0]) + lhs[3] * std::conj(rhs[1]);
		dest[3] += lhs[2] * std::conj(rhs[2]) + lhs[3] * std::conj(rhs[3]);
	}

	template<typename ComplType, typename LHS_T, typename RHS_T>
	static void HermATimesB(std::complex<ComplType> *dest, const LHS_T* lhs, const RHS_T* rhs)
	{
		dest[0] = std::conj(lhs[0]) * rhs[0] + std::conj(lhs[2]) * rhs[2];
		dest[1] = std::conj(lhs[0]) * rhs[1] + std::conj(lhs[2]) * rhs[3];
		dest[2] = std::conj(lhs[1]) * rhs[0] + std::conj(lhs[3]) * rhs[2];
		dest[3] = std::conj(lhs[1]) * rhs[1] + std::conj(lhs[3]) * rhs[3];
	}

	static void HermATimesHermB(std::complex<double> *dest, const std::complex<double> *lhs, const std::complex<double> *rhs)
	{
		dest[0] = std::conj(lhs[0]) * std::conj(rhs[0]) + std::conj(lhs[2]) * std::conj(rhs[1]);
		dest[1] = std::conj(lhs[0]) * std::conj(rhs[2]) + std::conj(lhs[2]) * std::conj(rhs[3]);
		dest[2] = std::conj(lhs[1]) * std::conj(rhs[0]) + std::conj(lhs[3]) * std::conj(rhs[1]);
		dest[3] = std::conj(lhs[1]) * std::conj(rhs[2]) + std::conj(lhs[3]) * std::conj(rhs[3]);
	}
	
	template<typename ComplType, typename LHS_T, typename RHS_T>
	static void PlusHermATimesB(std::complex<ComplType> *dest, const LHS_T* lhs, const RHS_T* rhs)
	{
		dest[0] += std::conj(lhs[0]) * rhs[0] + std::conj(lhs[2]) * rhs[2];
		dest[1] += std::conj(lhs[0]) * rhs[1] + std::conj(lhs[2]) * rhs[3];
		dest[2] += std::conj(lhs[1]) * rhs[0] + std::conj(lhs[3]) * rhs[2];
		dest[3] += std::conj(lhs[1]) * rhs[1] + std::conj(lhs[3]) * rhs[3];
	}

	template<typename T>
	static bool Invert(T* matrix)
	{
		T d = ((matrix[0]*matrix[3]) - (matrix[1]*matrix[2]));
		if(d == 0.0)
			return false;
		T oneOverDeterminant = 1.0 / d;
		T temp;
		temp      = matrix[3] * oneOverDeterminant;
		matrix[1] = -matrix[1] * oneOverDeterminant;
		matrix[2] = -matrix[2] * oneOverDeterminant;
		matrix[3] = matrix[0] * oneOverDeterminant;
		matrix[0] = temp;
		return true;
	}

	static bool MultiplyWithInverse(std::complex<double>* lhs, const std::complex<double>* rhs)
	{
		std::complex<double> d = ((rhs[0]*rhs[3]) - (rhs[1]*rhs[2]));
		if(d == 0.0) return false;
		std::complex<double> oneOverDeterminant = 1.0 / d;
		std::complex<double> temp[4];
		temp[0] = rhs[3] * oneOverDeterminant;
		temp[1] = -rhs[1] * oneOverDeterminant;
		temp[2] = -rhs[2] * oneOverDeterminant;
		temp[3] = rhs[0] * oneOverDeterminant;
		
		std::complex<double> temp2 = lhs[0];
		lhs[0] = lhs[0] * temp[0] + lhs[1] * temp[2];
		lhs[1] =  temp2 * temp[1] + lhs[1] * temp[3];
		
		temp2 = lhs[2];
		lhs[2] = lhs[2] * temp[0] + lhs[3] * temp[2];
		lhs[3] = temp2 * temp[1] + lhs[3] * temp[3];
		return true;
	}

	static void SingularValues(const std::complex<double>* matrix, double &e1, double &e2)
	{
		// This is not the ultimate fastest method, since we
		// don't need to calculate the imaginary values of b,c at all.
		// Calculate M M^H
		std::complex<double> temp[4] = {
			matrix[0] * std::conj(matrix[0]) + matrix[1] * std::conj(matrix[1]),
			matrix[0] * std::conj(matrix[2]) + matrix[1] * std::conj(matrix[3]),
			matrix[2] * std::conj(matrix[0]) + matrix[3] * std::conj(matrix[1]),
			matrix[2] * std::conj(matrix[2]) + matrix[3] * std::conj(matrix[3])
		};
		// Use quadratic formula, with a=1.
		double
			b = -temp[0].real() - temp[3].real(),
			c = temp[0].real()*temp[3].real() - (temp[1]*temp[2]).real(),
			d = b*b - (4.0*1.0)*c,
			sqrtd = sqrt(d);

		e1 = sqrt((-b + sqrtd) * 0.5);
		e2 = sqrt((-b - sqrtd) * 0.5);
	}
	
	static void EigenValues(const double* matrix, double &e1, double &e2)
	{
		double tr = matrix[0] + matrix[3];
		double d = matrix[0]*matrix[3] - matrix[1]*matrix[2];
		double term = sqrt(tr*tr*0.25-d);
		double trHalf = tr*0.5;
		e1 = trHalf + term;
		e2 = trHalf - term;
	}
	
	static void EigenValues(const std::complex<double>* matrix, std::complex<double> &e1, std::complex<double> &e2)
	{
		std::complex<double> tr = matrix[0] + matrix[3];
		std::complex<double> d = matrix[0]*matrix[3] - matrix[1]*matrix[2];
		std::complex<double> term = sqrt(tr*tr*0.25-d);
		std::complex<double> trHalf = tr*0.5;
		e1 = trHalf + term;
		e2 = trHalf - term;
	}
	
	static void EigenValuesAndVectors(const double* matrix, double &e1, double &e2, double* vec1, double* vec2)
	{
		double tr = matrix[0] + matrix[3];
		double d = matrix[0]*matrix[3] - matrix[1]*matrix[2];
		double term = sqrt(tr*tr*0.25-d);
		double trHalf = tr*0.5;
		e1 = trHalf + term;
		e2 = trHalf - term;
		if(matrix[2] != 0.0)
		{
			vec1[0] = matrix[3] - e1;
			vec1[1] = -matrix[2];
			vec2[0] = matrix[3] - e2;
			vec2[1] = -matrix[2];
		}
		else if(matrix[1] != 0.0)
		{
			vec1[0] = -matrix[1];
			vec1[1] = matrix[0] - e1;
			vec2[0] = -matrix[1];
			vec2[1] = matrix[0] - e2;
		}
		else {
			vec1[0] = 1.0;
			vec1[1] = 0.0;
			vec2[0] = 0.0;
			vec2[1] = 1.0;
		}
	}
	
	static void SquareRoot(double* matrix)
	{
		double tr = matrix[0] + matrix[3];
		double d = matrix[0]*matrix[3] - matrix[1]*matrix[2];
		double s = /*+/-*/ sqrt(d);
		double t = /*+/-*/ sqrt(tr + 2.0*s);
		if(t != 0.0)
		{
			matrix[0] = (matrix[0]+s ) / t;
			matrix[1] = (matrix[1] / t);
			matrix[2] = (matrix[2] / t);
			matrix[3] = (matrix[3]+s) / t;
		}
		else {
			if(matrix[0] == 0.0 && matrix[1] == 0.0 &&
				matrix[2] == 0.0 && matrix[3] == 0.0)
			{
				// done: it's the zero matrix
			} else {
				for(size_t i=0; i!=4; ++i)
					matrix[i] = std::numeric_limits<double>::quiet_NaN();
			}
		}
	}
	
	template<typename T>
	static T RotationAngle(const std::complex<T>* matrix)
	{
		return std::atan2((matrix[2].real()-matrix[1].real())*0.5, (matrix[0].real()+matrix[3].real())*0.5);
	}
	
	template<typename T>
	static void RotationMatrix(std::complex<T>* matrix, double alpha)
	{
		T cosAlpha, sinAlpha;
		sincos(alpha, &sinAlpha, &cosAlpha);
		matrix[0] = cosAlpha; matrix[1] = -sinAlpha;
		matrix[2] = sinAlpha; matrix[3] = cosAlpha;
	}
};

class MC2x2
{
public:
	MC2x2() { }
	MC2x2(const MC2x2& source) { Matrix2x2::Assign(_values, source._values); }
	template<typename T>
	explicit MC2x2(const T source[4]) { Matrix2x2::Assign(_values, source); }
	MC2x2(double m00, double m01, double m10, double m11) {
		_values[0] = m00; _values[1] = m01;
		_values[2] = m10; _values[3] = m11;
	}
	MC2x2(std::complex<double> m00, std::complex<double> m01, std::complex<double> m10, std::complex<double> m11) {
		_values[0] = m00; _values[1] = m01;
		_values[2] = m10; _values[3] = m11;
	}
	MC2x2& operator=(const MC2x2& source)
	{
		Matrix2x2::Assign(_values, source._values);
		return *this;
	}
	MC2x2& operator+=(const MC2x2& rhs)
	{
		Matrix2x2::Add(_values, rhs._values);
		return *this;
	}
	MC2x2& operator*=(double rhs)
	{
		Matrix2x2::ScalarMultiply(_values, rhs);
		return *this;
	}
	MC2x2& operator/=(double rhs)
	{
		Matrix2x2::ScalarMultiply(_values, 1.0/rhs);
		return *this;
	}
	const std::complex<double>& operator[](size_t index) const { return _values[index]; }
	std::complex<double>& operator[](size_t index) { return _values[index]; }
	const double& IndexReal(size_t index) const { return reinterpret_cast<const double(&)[2]>(_values[index/2])[index%2]; }
	double& IndexReal(size_t index) { return reinterpret_cast<double(&)[2]>(_values[index/2])[index%2]; }
	static MC2x2 Zero()
	{
		return MC2x2(0.0, 0.0, 0.0, 0.0);
	}
	static MC2x2 Unity()
	{
		return MC2x2(1.0, 0.0, 0.0, 1.0);
	}
	static MC2x2 NaN()
	{
		return MC2x2(
			std::complex<double>(std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN()),
			std::complex<double>(std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN()),
			std::complex<double>(std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN()),
			std::complex<double>(std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN()));
	}
	std::complex<double>* Data() { return _values; }
	const std::complex<double>* Data() const { return _values; }
	MC2x2 Multiply(const MC2x2& rhs) const
	{
		MC2x2 dest;
		Matrix2x2::ATimesB(dest._values, _values, rhs._values);
		return dest;
	}
	MC2x2 MultiplyHerm(const MC2x2& rhs) const
	{
		MC2x2 dest;
		Matrix2x2::ATimesHermB(dest._values, _values, rhs._values);
		return dest;
	}
	MC2x2 HermThenMultiply(const MC2x2& rhs) const
	{
		MC2x2 dest;
		Matrix2x2::HermATimesB(dest._values, _values, rhs._values);
		return dest;
	}
	MC2x2 HermThenMultiplyHerm(const MC2x2& rhs) const
	{
		MC2x2 dest;
		Matrix2x2::HermATimesHermB(dest._values, _values, rhs._values);
		return dest;
	}
	void AddWithFactorAndAssign(const MC2x2& rhs, double factor)
	{
		Matrix2x2::MultiplyAdd(_values, rhs._values, factor);
	}
	bool Invert()
	{
		return Matrix2x2::Invert(_values);
	}
	static void ATimesB(MC2x2& dest, const MC2x2& lhs, const MC2x2& rhs)
	{
		Matrix2x2::ATimesB(dest._values, lhs._values, rhs._values);
	}
	static void ATimesHermB(MC2x2& dest, const MC2x2& lhs, const MC2x2& rhs)
	{
		Matrix2x2::ATimesHermB(dest._values, lhs._values, rhs._values);
	}
	static void HermATimesB(MC2x2& dest, const MC2x2& lhs, const MC2x2& rhs)
	{
		Matrix2x2::HermATimesB(dest._values, lhs._values, rhs._values);
	}
	static void HermATimesHermB(MC2x2& dest, const MC2x2& lhs, const MC2x2& rhs)
	{
		Matrix2x2::HermATimesHermB(dest._values, lhs._values, rhs._values);
	}
	std::string ToString() const
	{
		std::stringstream str;
		str << _values[0] << ", " << _values[1] << "; "
			<< _values[2] << ", " << _values[3];
		return str.str();
	}
	void CopyValues(std::complex<double>* values) const
	{
		Matrix2x2::Assign(values, _values);
	}
	void EigenValues(std::complex<double> &e1, std::complex<double> &e2) const
	{
		Matrix2x2::EigenValues(_values, e1, e2);
	}
	bool HasNaN() const
	{
		return !(
			std::isfinite(_values[0].real()) && std::isfinite(_values[0].imag()) &&
			std::isfinite(_values[1].real()) && std::isfinite(_values[1].imag()) &&
			std::isfinite(_values[2].real()) && std::isfinite(_values[2].imag()) &&
			std::isfinite(_values[3].real()) && std::isfinite(_values[3].imag())
		);
	}
private:
	std::complex<double> _values[4];
};

#endif
