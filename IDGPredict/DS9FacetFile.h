#ifndef DS9_FACET_FILE
#define DS9_FACET_FILE

#include <fstream>
#include <vector>

#include "FacetMap.h"

#include "../Common/ImageCoordinates.h"

class DS9FacetFile
{
public:
	enum TokenType { Empty, Word, Number, Symbol, Comment };
	
	DS9FacetFile(const std::string& filename) : _file(filename), _hasChar(false)
	{ 
		Skip();
	}
	
	void Read(FacetMap& map, double phaseCentreRA, double phaseCentreDec, double pxScaleX, double pxScaleY, size_t width, size_t height)
	{
    Facet* lastFacet = nullptr;
		while(Type() != DS9FacetFile::Empty)
		{
			std::string t = Token();
			if(t == "global" || t == "fk5")
				SkipLine();
			else {
				Skip();
				
				if(t == "polygon")
				{
					std::vector<double> vals = readNumList();
					if(vals.size()%2!=0)
						throw std::runtime_error("Polygon is expecting an even number of numbers in its list");
					std::vector<double>::const_iterator i = vals.begin();
					Facet& facet = map.AddFacet();
					while(i != vals.end())
					{
						double ra = *i*(M_PI/180.0);
						++i;
						double dec = *i*(M_PI/180.0);
						++i;
						double l, m;
						ImageCoordinates::RaDecToLM(ra, dec, phaseCentreRA, phaseCentreDec, l, m);
						int x, y;
						ImageCoordinates::LMToXY(l, m, pxScaleX, pxScaleY, width, height, x, y);
						facet.AddVertex(x, y);
					}
					lastFacet = &facet;
				}
				else if(t == "point")
				{
					std::vector<double> vals = readNumList();
					if(vals.size()!=2)
						throw std::runtime_error("Point is expecting exactly two numbers in its list");
          if(lastFacet)
          {
            lastFacet->SetRA(vals[0]*(M_PI/180.0));
            lastFacet->SetDec(vals[1]*(M_PI/180.0));
          }
				}
			}
		}
	}
	
	std::vector<double> readNumList()
	{
		std::vector<double> vals;
		if(Token()!="(")
			throw std::runtime_error("Excepting '(' after polygon keyword");
		Skip();
		while(Token()!=")") {
			if(Type() != Number)
				throw std::runtime_error("Excepted number or ')' after '(' ");
			vals.push_back(atof(Token().c_str()));
			Skip();
			if(Token()==",")
				Skip();
		}
		Skip();
		return vals;
	}
	
	std::string Token() const
	{
		return _token;
	}
	
	TokenType Type() const
	{
		return _type;
	}
	
	void SkipLine()
	{
		char c;
		while(nextChar(c))
		{
			if(c == '\n')
				break;
		}
		Skip();
	}
	
	void Skip()
	{
		bool cont = true;
		_type = Empty;
		_token = std::string();
		do {
			char c;
			if(nextChar(c))
			{
				switch(_type)
				{
				case Empty:
					if(isAlpha(c))
					{
						_type = Word;
						_token += c;
					}
					else if(isWhiteSpace(c))
					{ }
					else if(isNumeric(c))
					{
						_type = Number;
						_token += c;
					}
					else if(c == '(' || c ==')' || c == ',')
					{
						_type = Symbol;
						_token += c;
						cont = false;
					}
					else if(c == '#')
					{
						_type = Comment;
					}
					break;
				case Word:
					if(isAlpha(c) || (c>='0' && c<='9'))
					{
						_token += c;
					}
					else {
						cont = false;
						pushChar(c);
					}
					break;
				case Number:
					if(isNumeric(c))
					{
						_token += c;
					}
					else {
						cont = false;
						pushChar(c);
					}
					break;
				case Symbol:
					pushChar(c);
					cont = false;
					break;
				case Comment:
					if(c == '\n')
						_type = Empty;
					break;
				}
			}
			else {
				cont = false;
			}
		} while(cont);
	}
	
private:
	bool nextChar(char& c)
	{
		if(_hasChar)
		{
			_hasChar = false;
			c = _char;
			return true;
		}
		else {
			_file.read(&c, 1);
			return _file.good();
		}
	}
	void pushChar(char c)
	{
		_hasChar = true;
		_char = c;
	}
	
	std::ifstream _file;
	std::string _token;
	TokenType _type;
	bool _hasChar;
	char _char;
	
	static bool isAlpha(char c) {
		return ((c>='A' && c<='Z') || (c>='a' && c<='z') || c=='_');
	}
	static bool isWhiteSpace(char c) {
		return c==' ' || c=='\n' || c=='\r' || c=='\t';
	}
	static bool isNumeric(char c) {
		return (c>='0' && c<='9') || c=='-' || c=='.';
	}
};

#endif
