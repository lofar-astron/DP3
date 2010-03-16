#include <stack>
#include <vector>
#include <iostream>
#include <Common/StringUtil.h>
#include <Common/StreamUtil.h>
#include <Common/LofarLogger.h>

using namespace LOFAR;
using namespace std;

    void parse (const string& origExpr)
    {
      string expr = toUpper(origExpr);
      stack<int> tokens;
      vector<int> rpn;
      vector<string> names;
      unsigned i=0;
      bool hadName = false;    // the last token was a name.
      while (i < expr.size()) {
        int oper = 1;
        if (expr[i] == ' '  ||  expr[i] == '\t') {
          i++;
        } else if (expr[i] == '(') {
          ASSERTSTR (!hadName, "no operator before opening parenthesis at pos. "
                     << i << " in expression: " << origExpr);
          oper = 0;
          tokens.push (oper);
          i++;
        } else if (expr[i] == '|') {
          oper = -1;
          i++;
          if (i < expr.size()  &&  expr[i] == '|') i++;
        } else if (expr[i] == ',') {
          oper = -1;
          i++;
        } else if (expr[i] == '&') {
          oper = -2;
          i++;
          if (i < expr.size()  &&  expr[i] == '&') i++;
        } else if (expr[i] == '!') {
          oper = -3;
          i++;
        } else if (expr.size()-i >= 3  &&  (expr.substr(i,3) == "OR "  ||
                                            expr.substr(i,3) == "OR\t" ||
                                            expr.substr(i,3) == "OR!"  ||
                                            expr.substr(i,3) == "OR(")) {
          oper = -1;
          i+=2;
        } else if (expr.size()-i >= 4  &&  (expr.substr(i,4) == "AND "  ||
                                            expr.substr(i,4) == "AND\t" ||
                                            expr.substr(i,4) == "AND!"  ||
                                            expr.substr(i,4) == "AND(")) {
          oper = -2;
          i+=3;
        } else if (expr.size()-i >= 4  &&  (expr.substr(i,4) == "NOT "  ||
                                            expr.substr(i,4) == "NOT\t" ||
                                            expr.substr(i,4) == "NOT(")) {
          oper = -3;
          i+=3;
        } else if (expr[i] == ')') {
          ASSERTSTR (hadName, "no set name before closing parenthesis at pos. "
                     << i << " in expression: " << origExpr);
          while (true) {
            ASSERTSTR (!tokens.empty(), "mismatched parentheses at pos. "
                       << i << " in expression: " << origExpr);
            if (tokens.top() == 0) {
              tokens.pop();
              break;
            }
            rpn.push_back (tokens.top());
            tokens.pop();
          }
          i++;
        } else {
          int st=i;
          while (i < expr.size() && 
                 expr[i] != ' ' && expr[i] != '\t' &&
                 expr[i] != '(' && expr[i] != ')' && expr[i] != '!' &&
                 expr[i] != ',' && expr[i] != '&' && expr[i] != '|') {
            i++;
          }
          ASSERTSTR (!hadName, "No operator between set names at pos. "
                     << i << " in expression: " << origExpr);
          hadName = true;
          rpn.push_back (names.size());
          names.push_back (origExpr.substr(st, i-st));
        }
        if (oper < 0) {
          if (oper == -3) {
            ASSERTSTR (!hadName, "No set name before operator ! at pos. "
                       << i << " in expression: " << origExpr);
          } else {
            ASSERTSTR (hadName, "No set name before operator at pos. "
                       << i << " in expression: " << origExpr);
          }
          hadName = false;
          while (!tokens.empty()  &&  tokens.top() <= oper) {
            rpn.push_back (tokens.top());
            tokens.pop();
          }
          tokens.push (oper);
        }
      }
      ASSERTSTR (hadName, "no set name after last operator in expression: "
                 << origExpr);
      while (!tokens.empty()) {
        ASSERTSTR (tokens.top()<0, "mismatched parentheses in expression: "
                   << origExpr);
        rpn.push_back (tokens.top());
        tokens.pop();
      }

      cout << rpn << endl;
      cout << names << endl;
    }

int main(int argc, char* argv[])
{
  try {
    if (argc > 1) {
      parse (argv[1]);
    }
  } catch (std::exception& x) {
    cout << x.what() << endl;
  }
}
