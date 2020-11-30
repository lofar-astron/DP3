#include "../../StringUtil.h"

#include <boost/test/unit_test.hpp>

#include <iomanip>
#include <limits>
#include <type_traits>

using StringList = std::vector<std::string>;

BOOST_AUTO_TEST_SUITE(stringutil)

BOOST_AUTO_TEST_CASE(tokenize) {
  using DP3::StringUtil::tokenize;

  BOOST_TEST(tokenize("", "") == StringList());
  BOOST_TEST(tokenize("foo", "!") == StringList({"foo"}));
  BOOST_TEST(tokenize("foo!bar", "!?") == StringList({"foo", "bar"}));
  BOOST_TEST(tokenize("!!foo!?bar??", "!?") == StringList({"foo", "bar"}));
}

BOOST_AUTO_TEST_CASE(compare) {
  using DP3::StringUtil::Compare;
  const Compare normal;
  const Compare nocase(Compare::NOCASE);

  BOOST_TEST(normal("a", "b"));
  BOOST_TEST(!normal("a", "a"));
  BOOST_TEST(!normal("b", "a"));

  BOOST_TEST(nocase("a", "b"));
  BOOST_TEST(!nocase("a", "a"));
  BOOST_TEST(!nocase("b", "a"));

  BOOST_TEST(normal("B", "a"));
  BOOST_TEST(!nocase("B", "a"));
}

BOOST_AUTO_TEST_CASE(formatString) {
  BOOST_TEST(
      DP3::formatString("The %s to life and everything is %d!", "answer", 42) ==
      "The answer to life and everything is 42!");
}

BOOST_AUTO_TEST_CASE(skipws) {
  BOOST_TEST(DP3::lskipws("", 0, 0) == 0u);
  BOOST_TEST(DP3::lskipws("foo", 0, 3) == 0u);
  BOOST_TEST(DP3::lskipws("  foo", 0, 5) == 2u);
  BOOST_TEST(DP3::lskipws("  foo", 0, 1) == 1u);

  BOOST_TEST(DP3::rskipws("", 0, 0) == 0u);
  BOOST_TEST(DP3::rskipws("foo", 0, 3) == 3u);
  BOOST_TEST(DP3::rskipws("foo  ", 0, 5) == 3u);
  BOOST_TEST(DP3::rskipws("foo  ", 4, 5) == 4u);
}

BOOST_AUTO_TEST_CASE(skipQuoted) {
  BOOST_TEST(DP3::skipQuoted("'foo'", 0) == 5u);
  BOOST_TEST(DP3::skipQuoted("  \"foo\" bar", 2) == 7u);
  BOOST_CHECK_THROW(DP3::skipQuoted("'foo", 0), std::runtime_error);
}

BOOST_AUTO_TEST_CASE(skipBalanced) {
  BOOST_TEST(DP3::skipBalanced("[]", 0, 2, ']') == 2u);
  BOOST_TEST(DP3::skipBalanced("[1[23]45]]]", 0, 11, ']') == 9u);
  BOOST_TEST(DP3::skipBalanced("x1234x", 0, 6, 'x') == 6u);
  BOOST_TEST(DP3::skipBalanced("  x1234x  ", 2, 10, 'x') == 8u);

  // Quoted substrings are ignored:
  BOOST_TEST(DP3::skipBalanced("x 'xx' x", 0, 8, 'x') == 8u);
  BOOST_TEST(DP3::skipBalanced("[ \"]]\" ]", 0, 8, ']') == 8u);

  BOOST_CHECK_THROW(DP3::skipBalanced("x", 0, 1, 'x'), std::runtime_error);
  BOOST_CHECK_THROW(DP3::skipBalanced("[ [][]", 0, 6, ']'), std::runtime_error);
}

BOOST_AUTO_TEST_CASE(strToBool) {
  BOOST_TEST(DP3::strToBool("false") == false);
  BOOST_TEST(DP3::strToBool("0") == false);
  BOOST_TEST(DP3::strToBool("Njet") == false);
  BOOST_TEST(DP3::strToBool("TRUE") == true);
  BOOST_TEST(DP3::strToBool("111") == true);
  BOOST_TEST(DP3::strToBool("yep") == true);
  BOOST_CHECK_THROW(DP3::strToBool(""), std::runtime_error);
  BOOST_CHECK_THROW(DP3::strToBool("bar"), std::runtime_error);
}

namespace {

// std::to_string does not suffice for floating point values.
// -> Use std::stringstream indead, with increased precision.
template <typename T>
std::string ToString(T t) {
  std::stringstream ss;
  ss << std::setprecision(20) << t;
  return ss.str();
}

template <typename T>
void testMaxStrTo(T (*f)(const std::string&)) {
  const T max = std::numeric_limits<T>::max();
  std::string max_str = ToString(max);
  BOOST_TEST(f(max_str) == max);

  // Raise last digit, so the value become too large.
  // (For floats or doubles, raise the digit of the exponent.)
  char& last_digit = *(max_str.end() - 1);
  BOOST_TEST_REQUIRE(last_digit != '9');
  ++last_digit;
  BOOST_CHECK_THROW(f(max_str), std::runtime_error);
}

template <typename T>
void testMinStrTo(T (*f)(const std::string&)) {
  const T min = std::numeric_limits<T>::min();
  std::string min_str = ToString(min);
  BOOST_TEST(f(min_str) == min);

  // Raise last digit, so the value become too small.
  // (For floats or doubles, raise the digit of the exponent.)
  char& last_digit = *(min_str.end() - 1);
  BOOST_TEST_REQUIRE(last_digit != '9');
  ++last_digit;
  BOOST_CHECK_THROW(f(min_str), std::runtime_error);
}

template <typename T>
void testStrToInt(T (*f)(const std::string&)) {
  // Basic tests:
  BOOST_TEST(f("42") == T(42));
  if (std::is_unsigned<T>()) {
    BOOST_CHECK_THROW(f("-42"), std::runtime_error);
  } else {
    BOOST_TEST(f("-42") == T(-42));
  }
  BOOST_TEST(f(" 0x2A ") == T(42));
  BOOST_CHECK_THROW(f("  foo  "), std::runtime_error);
  BOOST_CHECK_THROW(f("-0x0"), std::runtime_error);
  BOOST_CHECK_THROW(f("42bar"), std::runtime_error);
  BOOST_CHECK_THROW(f("bar42"), std::runtime_error);

  // Test limit values:
  testMaxStrTo(f);
  if (!std::is_unsigned<T>()) testMinStrTo(f);
}

template <typename T>
void testStrToFloat(T (*f)(const std::string&)) {
  BOOST_TEST(f("42") == T(42.0));
  BOOST_TEST(f(" -42.0 ") == T(-42.0));
  BOOST_TEST(f("  42.25") == T(42.25));
  BOOST_TEST(f("-1e4  ") == T(-1e4));
  BOOST_CHECK_THROW(f("  foo  "), std::runtime_error);
  BOOST_CHECK_THROW(f("42bar"), std::runtime_error);
  BOOST_CHECK_THROW(f("bar42"), std::runtime_error);

  testMaxStrTo(f);
  testMinStrTo(f);
}

}  // namespace

BOOST_AUTO_TEST_CASE(strToLong) { testStrToInt(DP3::strToLong); }
BOOST_AUTO_TEST_CASE(strToInt) { testStrToInt(DP3::strToInt); }
BOOST_AUTO_TEST_CASE(strToInt16) { testStrToInt(DP3::strToInt16); }
BOOST_AUTO_TEST_CASE(strToInt32) { testStrToInt(DP3::strToInt32); }
BOOST_AUTO_TEST_CASE(strToInt64) { testStrToInt(DP3::strToInt64); }
BOOST_AUTO_TEST_CASE(strToUlong) { testStrToInt(DP3::strToUlong); }
BOOST_AUTO_TEST_CASE(strToUint) { testStrToInt(DP3::strToUint); }
BOOST_AUTO_TEST_CASE(strToUint16) { testStrToInt(DP3::strToUint16); }
BOOST_AUTO_TEST_CASE(strToUint32) { testStrToInt(DP3::strToUint32); }
BOOST_AUTO_TEST_CASE(strToUint64) { testStrToInt(DP3::strToUint64); }
BOOST_AUTO_TEST_CASE(strToFloat) { testStrToFloat(DP3::strToFloat); }
BOOST_AUTO_TEST_CASE(strToDouble) { testStrToFloat(DP3::strToDouble); }

BOOST_AUTO_TEST_CASE(expandRangeString) {
  BOOST_TEST(DP3::expandRangeString("") == "");
  BOOST_TEST(DP3::expandRangeString("foo, bar") == "foo, bar");
  BOOST_TEST(DP3::expandRangeString("1..4") == "1,2,3,4");
  BOOST_TEST(DP3::expandRangeString("'1..4'") == "'1..4'");
  BOOST_CHECK_THROW(DP3::expandMultString("' unbalanced quotes"),
                    std::runtime_error);
  BOOST_TEST(DP3::expandRangeString("[1..4]") == "[1,2,3,4]");
  BOOST_TEST(DP3::expandRangeString("{1..4}") == "1,2,3,4");
  BOOST_TEST(DP3::expandRangeString("001..001") == "001");
  BOOST_TEST(DP3::expandRangeString("xx042..044yy, ,,'1..4',foo1..2bar") ==
             "xx042yy,xx043yy,xx044yy, ,,'1..4',foo1bar,foo2bar");
}

BOOST_AUTO_TEST_CASE(expandMultString) {
  BOOST_TEST(DP3::expandMultString("") == "");
  BOOST_TEST(DP3::expandMultString("1*42") == "42");
  BOOST_TEST(DP3::expandMultString("11*x") == "x,x,x,x,x,x,x,x,x,x,x");
  BOOST_TEST(DP3::expandMultString("0*foo") == "");
  BOOST_TEST(DP3::expandMultString("'1*42'") == "'1*42'");
  BOOST_CHECK_THROW(DP3::expandMultString("' unbalanced quotes"),
                    std::runtime_error);
  BOOST_TEST(DP3::expandMultString("[ 3*bar]") == "[ bar,bar,bar]");
  BOOST_TEST(DP3::expandMultString("{3*bar}") == "{3*bar}");
  BOOST_TEST(DP3::expandMultString(" ( 3*bar ) ") == " ( bar ,bar ,bar ) ");
}

BOOST_AUTO_TEST_CASE(expandArrayString) {
  BOOST_TEST(DP3::expandArrayString("") == "");
  BOOST_TEST(DP3::expandArrayString("3*5") == "3*5");
  BOOST_TEST(DP3::expandArrayString("[1..4]") == "[1,2,3,4]");
  BOOST_TEST(DP3::expandArrayString("[3*5]") == "[5,5,5]");
  BOOST_TEST(DP3::expandArrayString("[1..4,3*5]") == "[1,2,3,4,5,5,5]");
  BOOST_TEST(DP3::expandArrayString("[3*5,1..4]") == "[5,5,5,1,2,3,4]");
}

BOOST_AUTO_TEST_SUITE_END()
