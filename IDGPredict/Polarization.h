#ifndef POLARIZATION_H
#define POLARIZATION_H

#include <complex>
#include <stdexcept>
#include <set>
#include <vector>

/**
 * Class for various simple polarization related values.
 * 
 * The visibility relations for converting polarizations are:
 * 
 *   RR = I + V  ;   I = (RR + LL)/2
 *   RL = Q + iU ;   Q = (RL + LR)/2
 *   LR = Q - iU ;   U = -i (RL - LR)/2
 *   LL = I - V  ;   V = (RR - LL)/2
 *
 *   XX = I + Q  ;   I = (XX + YY)/2
 *   XY = U + iV ;   Q = (XX - YY)/2
 *   YX = U - iV ;   U = (XY + YX)/2
 *   YY = I - Q  ;   V = -i(XY - YX)/2
 * 
 * These definitions assume that 'X' and 'Y' are labelled as they are in
 * CASA measurement sets: X is North-South.
 * 
 * Note that uv-fits files have the polarizations labelled such that X
 * is East-West, confusingly. We have also observed that at present (2016)
 * CASA's importuvfits does not swap the polarizations, hence Measurement
 * Sets exist that have the polarizations labelled wrong.
 * The IEEE definition is to have X be N-S.
 */
class Polarization
{
public:
	enum PolarizationEnum
	{
		StokesI,
		StokesQ,
		StokesU,
		StokesV,
		RR,
		RL,
		LR,
		LL,
		XX,
		XY,
		YX,
		YY,
		/**
		 * Instrumental is a special value representing that four polarizations are
		 * stored, and these are the 'raw' measurement set polarizations. It is used
		 * as a special value that e.g. can be passed to a MSProvider, which would mean
		 * all values are gridded at once, as required for an a-term correcting gridder.
		 */
		Instrumental
	};
	
	static size_t TypeToIndex(enum PolarizationEnum polarization, size_t polCountInSet)
	{
		switch(polCountInSet)
		{
			case 1:
				if(polarization != StokesI)
					throw std::runtime_error("TypeTo4PolIndex(): can't convert given polarization to index");
				else
					return 0;
			case 2:
			switch(polarization)
			{
				case XX: return 0;
				case YY: return 1;
				default: throw std::runtime_error("TypeTo4PolIndex(): can't convert given polarization to index");
			}
			case 4:
			switch(polarization)
			{
				case XX: return 0;
				case XY: return 1;
				case YX: return 2;
				case YY: return 3;
				default: throw std::runtime_error("TypeTo4PolIndex(): can't convert given polarization to index");
			}
			default: throw std::runtime_error("TypeTo4PolIndex(): can't convert given polarization to index");
		}
	}
	
	static enum PolarizationEnum AipsIndexToEnum(int index)
	{
		switch(index)
		{
			case 1: return StokesI;
			case 2: return StokesQ;
			case 3: return StokesU;
			case 4: return StokesV;
			case 5: return RR;
			case 6: return RL;
			case 7: return LR;
			case 8: return LL;
			case 9: return XX;
			case 10: return XY;
			case 11: return YX;
			case 12: return YY;
			default: throw std::runtime_error("AipsIndexToEnum(): unknown aips polarization index");
		}
	}
	
	static bool TypeToIndex(enum PolarizationEnum polarization, const std::vector<PolarizationEnum>& polList, size_t& index)
	{
		for(size_t i=0; i!=polList.size(); ++i)
		{
			if(polList[i] == polarization)
			{
				index = i;
				return true;
			}
		}
		return false;
	}
	
	template<typename Range>
	static bool TypeToIndex(enum PolarizationEnum polarization, const Range& polList, size_t& index)
	{
		size_t curIndex = 0;
		for(typename Range::const_iterator i=polList.begin(); i!=polList.end(); ++i, ++curIndex)
		{
			if(*i == polarization)
			{
				index = curIndex;
				return true;
			}
		}
		return false;
	}
	
	static size_t StokesToIndex(enum PolarizationEnum polarization)
	{
		switch(polarization)
		{
			default:
			case StokesI: return 0;
			case StokesQ: return 1;
			case StokesU: return 2;
			case StokesV: return 3;
		}
	}
	
	static PolarizationEnum IndexToStokes(size_t index)
	{
		const static PolarizationEnum arr[4] = { StokesI, StokesQ, StokesU, StokesV };
		return arr[index];
	}
	
	static bool IsStokes(enum PolarizationEnum polarization)
	{
		return
			polarization == StokesI || polarization == StokesQ ||
			polarization == StokesU || polarization == StokesV;
	}
	
	static size_t TypeTo4PolIndex(enum PolarizationEnum polarization)
	{
		switch(polarization)
		{
			case XX: return 0;
			case XY: return 1;
			case YX: return 2;
			case YY: return 3;
			default: throw std::runtime_error("TypeTo4PolIndex(): can't convert given polarization to index");
		}
	}
	
	static std::string TypeToShortString(enum PolarizationEnum polarization)
	{
		switch(polarization)
		{
			case XX: return "XX";
			case XY: return "XY";
			case YX: return "YX";
			case YY: return "YY";
			case StokesI: return "I";
			case StokesQ: return "Q";
			case StokesU: return "U";
			case StokesV: return "V";
			case RR: return "RR";
			case RL: return "RL";
			case LR: return "LR";
			case LL: return "LL";
			case Instrumental: return "instr";
			default: return "";
		}
	}
	
	static std::string TypeToFullString(enum PolarizationEnum polarization)
	{
		switch(polarization)
		{
			case XX: return "XX";
			case XY: return "XY";
			case YX: return "YX";
			case YY: return "YY";
			case StokesI: return "Stokes I";
			case StokesQ: return "Stokes Q";
			case StokesU: return "Stokes U";
			case StokesV: return "Stokes V";
			case RR: return "RR";
			case RL: return "RL";
			case LR: return "LR";
			case LL: return "LL";
			case Instrumental: return "instrumental";
			default: return "Unknown polarization";
		}
	}
	
	static bool HasFullPolarization(const std::set<PolarizationEnum>& polarizations)
	{
		return
			HasFullLinearPolarization(polarizations) || HasFullStokesPolarization(polarizations)
			||
			(polarizations.count(RR)>0 && polarizations.count(RL)>0 &&
			 polarizations.count(LR)>0 && polarizations.count(LL)>0);
	}
	
	static bool HasFullLinearPolarization(const std::set<PolarizationEnum>& polarizations)
	{
		return
			(polarizations.count(XX)>0 && polarizations.count(XY)>0 &&
			 polarizations.count(YX)>0 && polarizations.count(YY)>0);
	}
	
	static bool HasFullCircularPolarization(const std::set<PolarizationEnum>& polarizations)
	{
		return
			(polarizations.count(LL)>0 && polarizations.count(LR)>0 &&
			 polarizations.count(RL)>0 && polarizations.count(RR)>0);
	}
	
	static bool HasFullStokesPolarization(const std::set<PolarizationEnum>& polarizations)
	{
		return
			(polarizations.count(StokesI)>0 && polarizations.count(StokesQ)>0 &&
			 polarizations.count(StokesU)>0 && polarizations.count(StokesV)>0);
	}
	
	static bool HasDualPolarization(const std::set<PolarizationEnum>& polarizations)
	{
		return HasDualLinearPolarization(polarizations) ||
			(polarizations.count(RR)>0 && polarizations.count(LL)>0);
	}
	
	static bool HasDualLinearPolarization(const std::set<PolarizationEnum>& polarizations)
	{
		return (polarizations.count(XX)>0 && polarizations.count(YY)>0);
	}
	
	static bool IsComplex(PolarizationEnum polarization)
	{
		return polarization==XY || polarization==YX;
	}
	
	static PolarizationEnum ParseString(const std::string& str)
	{
		if(str == "XX") return XX;
		else if(str == "XY") return XY;
		else if(str == "YX") return YX;
		else if(str == "YY") return YY;
		else if(str == "I") return StokesI;
		else if(str == "Q") return StokesQ;
		else if(str == "U") return StokesU;
		else if(str == "V") return StokesV;
		else if(str == "RR") return RR;
		else if(str == "RL") return RL;
		else if(str == "LR") return LR;
		else if(str == "LL") return LL;
		else throw std::runtime_error(std::string("Could not parse polarization string: ") + str);
	}
	
	static std::set<PolarizationEnum> ParseList(const std::string& listStr)
	{
		std::set<PolarizationEnum> list;
		enum { StartSt, GotXSt, GotYSt, GotLSt, GotRSt, GotSeperatorSt } state = StartSt;
		for(std::string::const_iterator i=listStr.begin(); i!=listStr.end(); ++i)
		{
			char c = (*i>='a' && *i<='z') ? *i-('a'-'A') : *i;
			switch(c)
			{
				case 'X':
					if(state == StartSt || state == GotSeperatorSt) state=GotXSt;
					else {
						if(state==GotXSt) list.insert(XX);
						else if(state==GotYSt) list.insert(YX);
						else throw std::runtime_error("Invalid polarization list: parse error near 'X'");
						state=StartSt;
					}
					break;
				case 'Y':
					if(state == StartSt || state == GotSeperatorSt) state=GotYSt;
					else {
						if(state==GotXSt) list.insert(XY);
						else if(state==GotYSt) list.insert(YY);
						else throw std::runtime_error("Invalid polarization list: parse error near 'Y'");
						state=StartSt;
					}
					break;
				case 'R':
					if(state == StartSt || state == GotSeperatorSt) state=GotRSt;
					else {
						if(state==GotRSt) list.insert(RR);
						else if(state==GotLSt) list.insert(LR);
						else throw std::runtime_error("Invalid polarization list: parse error near 'R'");
						state=StartSt;
					}
					break;
				case 'L':
					if(state == StartSt || state == GotSeperatorSt) state=GotLSt;
					else {
						if(state==GotRSt) list.insert(RL);
						else if(state==GotLSt) list.insert(LL);
						else throw std::runtime_error("Invalid polarization list: parse error near 'L'");
						state=StartSt;
					}
					break;
				case 'I':
					if(state == StartSt || state == GotSeperatorSt) list.insert(StokesI);
					else throw std::runtime_error("Invalid polarization list: parse error near 'I'");
					state = StartSt;
					break;
				case 'Q':
					if(state == StartSt || state == GotSeperatorSt) list.insert(StokesQ);
					else throw std::runtime_error("Invalid polarization list: parse error near 'Q'");
					state = StartSt;
					break;
				case 'U':
					if(state == StartSt || state == GotSeperatorSt) list.insert(StokesU);
					else throw std::runtime_error("Invalid polarization list: parse error near 'U'");
					state = StartSt;
					break;
				case 'V':
					if(state == StartSt || state == GotSeperatorSt) list.insert(StokesV);
					else throw std::runtime_error("Invalid polarization list: parse error near 'V'");
					state = StartSt;
					break;
				case ',':
				case ' ':
				case '/':
					if(state == StartSt) state = GotSeperatorSt;
					else throw std::runtime_error("Invalid polarization list: parse error near seperator");
			}
		}
		if(state!=StartSt)
			throw std::runtime_error("Invalid polarization list: parse error near string end");
		return list;
	}
	
	template<typename NumType>
	static void LinearToStokes(const std::complex<NumType> *linear, NumType* stokes)
	{
		// Note comments about X and Y at the top of this file!
		stokes[0] = 0.5 * (linear[0].real() + linear[3].real());
		stokes[1] = 0.5 * (linear[0].real() - linear[3].real());
		stokes[2] = 0.5 * (linear[1].real() + linear[2].real());
		stokes[3] = 0.5 * (linear[1].imag() - linear[2].imag());
	}
	
	template<typename NumType>
	static void StokesToLinear(const NumType* stokes, std::complex<NumType> *linear)
	{
		linear[0] = stokes[0] + stokes[1];
		linear[1] = std::complex<NumType>(stokes[2], stokes[3]);
		linear[2] = std::complex<NumType>(stokes[2], -stokes[3]);
		linear[3] = stokes[0] - stokes[1];
	}
	
	template<typename NumType>
	static void CircularToStokes(const std::complex<NumType>* circular, NumType* stokes)
	{
		stokes[0] = (circular[0] + circular[3]).real()*0.5;
		stokes[1] = (circular[1] + circular[2]).real()*0.5;
		stokes[2] = (circular[1] - circular[2]).imag()*0.5;
		stokes[3] = (circular[0] - circular[3]).real()*0.5;
	}
};

typedef Polarization::PolarizationEnum PolarizationEnum;

#endif
