#include "stdafx.h"
#include "Express.h"

/* 初始化静态数据 */
string CExpress::ptrList[] = { "f(", ",", "+", "-", "*", "/", "%", "^", "(", ")", "#" };
int CExpress::ptrArgCnt[] = { 0, 1, 2, 2, 2, 2, 2, 2, 0, 0, 0 };

string CExpress::funList[] = {
	"sin(",  "cos(",  "tan(",
	"pow(", "pow2(", "pow3(",
	"max(",  "min(", "sqrt(",
	"log(",   "ln(",  "exp(",
};
int CExpress::funArgCnt[] = {
	1, 1, 1,
	2, 1, 1,
	2, 2, 1,
	1, 1, 1,
};

int CExpress::preceMap[][EXP_PTR_NUM] = {
  //{f(,  ,,   +,  -,  *,  /,  %,  ^,  (,  ),  #}
	{ -1, -1, -1, -1, -1, -1, -1, -1, -1,  1,  1 },  //f(
	{  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1 },  //,
	{ -1,  1,  1,  1, -1, -1, -1, -1, -1,  1,  1 },  //+
	{ -1,  1,  1,  1, -1, -1, -1, -1, -1,  1,  1 },  //-
	{ -1,  1,  1,  1,  1,  1,  1,  1, -1,  1,  1 },  //*
	{ -1,  1,  1,  1,  1,  1,  1,  1, -1,  1,  1 },  ///
	{ -1,  1,  1,  1,  1,  1,  1,  1, -1,  1,  1 },  //%
	{ -1,  1,  1,  1,  1,  1,  1,  1, -1,  1,  1 },  //^
	{ -1,  1, -1, -1, -1, -1, -1, -1, -1,  0,  1 },  //(
	{  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1 },  //)
	{ -1,  1, -1, -1, -1, -1, -1, -1, -1, -1,  0 },  //#
};

CExpress::CExpress(const string &expression)
{
	this->expression = expression + "#";
	this->token = "";
	this->pos = 0;
	this->length = 1 + expression.length();
}

CExpress::~CExpress()
{
}

int CExpress::getPtrIndex(const string &ptr)
{
	for (int i = 0; i < EXP_PTR_NUM; i++) 
	{
		if (ptrList[i] == ptr)
		{
			return i;
		}
	}
	return -1;
}

int CExpress::getFunIndex(const string &fun)
{
	for (int i = 0; i < EXP_FUN_NUM; i++)
	{
		if (funList[i] == fun) 
		{
			return i;
		}
	}
	return -1;
}

int CExpress::comaprePrece(const string &ptr1, const string &ptr2)
{
	int m = getPtrIndex(ptr1), n = getPtrIndex(ptr2);
	if (m == -1) {
		m = 0;
	}
	if (n == -1) {
		n = 0;
	}
	return preceMap[m][n];
}

double CExpress::calculate(const string &ptr, double arg[])
{
	switch (getPtrIndex(ptr))
	{
	case 1:
		return arg[0];
	case 2:
		return arg[0] + arg[1];
	case 3:
		return arg[0] - arg[1];
	case 4:
		return arg[0] * arg[1];
	case 5:
		return arg[0] / arg[1];
	case 6:
		return (int(arg[0]) % int(arg[1]));
	case 7:
		return pow(arg[0],arg[1]);
	}

	return 0;
}

double CExpress::callFun(const string &fun, double arg[])
{
	switch (getFunIndex(fun))
	{
	case 0:
		return sin(arg[0]);
	case 1:
		return cos(arg[0]);
	case 2:
		return tan(arg[0]);
	case 3:
		return pow(arg[0], arg[1]);
	case 4:
		return arg[0] * arg[0];
	case 5:
		return arg[0] * arg[0] * arg[0];
	case 6:
		return arg[0] > arg[1] ? arg[0] : arg[1];
	case 7:
		return arg[0] < arg[1] ? arg[0] : arg[1];
	case 8:
		return sqrt(arg[0]);
	case 9:
		return log10(arg[0]);
	case 10:
		return log(arg[0]);
	case 11:
		return exp(arg[0]);
	}
	return 0;
}

void CExpress::readToken()
{
	while (pos < length && expression[pos] == ' ')  //去掉前空格 
	{  
		++pos;
	}
	if (pos >= length) 
	{
		tkType = TKT_ENDSIGN;
		return;
	}
	int pos_t = pos;
	char ch = expression[pos_t++];
	char ch_n = pos_t < length ? expression[pos_t] : 0;

	if (isdigit(ch) || (ch == '-' && isdigit(ch_n) && tkType != TKT_NUMBER)) //判断数字 
	{  
		if (ch == '-')
		{
			++pos_t;
		}
		while (pos_t < length && isdigit(ch = expression[pos_t]))
		{
			++pos_t;
		}
		if (ch == '.') 
		{
			++pos_t;
			while (pos_t < length && isdigit(expression[pos_t]))
			{
				++pos_t;
			}
		}
		tkType = TKT_NUMBER;
	}
	else if (-1 != getPtrIndex(string(1, ch))) //判断运算符 
	{  
		tkType = TKT_OPERATOR;
	}
	else if (isalpha(ch)) //判断函数 
	{  
		while (pos_t < length && (isalnum(ch) || ch == '_'))
		{
			ch = expression[++pos_t];
		}
		if (ch == '(') 
		{
			++pos_t;
			if (-1 != getFunIndex(expression.substr(pos, pos_t - pos)))
			{
				tkType = TKT_FUNCTION;
			}
			else 
			{
				tkType = TKT_UNKNOW;
			}
		}
		else 
		{
			tkType = TKT_UNKNOW;
		}
	}
	else 
	{
		tkType = TKT_UNKNOW;
	}
	token = expression.substr(pos, pos_t - pos);
	pos = pos_t;
}

bool CExpress::getArg(stack<double> &opnd, double arg[], int n)
{
	if (int(opnd.size()) < n) 
	{
		return false;
	}
	for (int i = n - 1; i >= 0; --i) 
	{
		arg[i] = opnd.top();
		opnd.pop();
	}
	return true;
}

int CExpress::getVal(double &res)
{
	stack<string> optr;  //算符栈 
	stack<double> opnd;  //算数栈 

	optr.push("#");
	pos = 0;
	readToken();

	while (tkType != TKT_ENDSIGN || !optr.empty())
	{
#ifdef EXP_DEBUG
		std::cout << "TkT = " << tkType << ", ";
		std::cout << "Pos = " << pos << "/" << length << ", ";
		std::cout << "Token = '" << token << "'" << std::endl;
#endif

		if (tkType == TKT_UNKNOW) 
		{
			return -1; //未知符号 
		}

		if (tkType == TKT_NUMBER) 
		{
			opnd.push(atof(token.c_str()));
			readToken();
		}
		else 
		{
			int comRes = comaprePrece(optr.top(), token);
#ifdef EXP_DEBUG
			std::cout << "compare('" << optr.top() << "', '" << token << "') = " << comRes << std::endl;
#endif
			switch (comRes) 
			{
			case -1:
				{
					optr.push(token);
					readToken();
					break;
				}
			case 1:
				{
					std::string ptr = optr.top();
					optr.pop();

					int idx = getPtrIndex(ptr), argCnt;
					double arg[10], res;
					if (-1 != idx) 
					{
						argCnt = ptrArgCnt[idx];
						if (argCnt) 
						{
							if (!getArg(opnd, arg, argCnt)) 
							{
								return -2;//表达式错误 
							}
							res = calculate(ptr, arg);
							opnd.push(res);
						}
					}
					else 
					{
						idx = getFunIndex(ptr);
						argCnt = funArgCnt[idx];
						if (!getArg(opnd, arg, argCnt))
						{
							return -2;//表达式错误 
						}
						res = callFun(ptr, arg);
						opnd.push(res);
						readToken();
					}
#ifdef EXP_DEBUG
					if (argCnt) {
						std::cout << "('" << ptr << "', ";
						for (int i = 0; i < argCnt; ++i) {
							std::cout << arg[i];
							if (i + 1 < argCnt) {
								std::cout << ", ";
							}
						}
						std::cout << ") = " << res << std::endl;
					}
#endif
					break;
				}
			case 0:
				{	
					optr.pop();
					readToken();
					break;
				}
			}
		}
	}
	res = opnd.top();
	return 0;
}
