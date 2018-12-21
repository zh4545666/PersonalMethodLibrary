#pragma once

#include "stdafx.h"

#define EXP_PTR_NUM 11  //运算符个数
#define EXP_FUN_NUM 12  //函数个数

#ifdef AFX_CLASS  
#define AFX_EX_CLASS _declspec(dllexport)  
#else  
#define AFX_EX_CLASS _declspec(dllimport)  
#endif 

/* 单词符号类型 */
typedef enum{
	TKT_NUMBER,    //数字 
	TKT_OPERATOR,  //运算符 
	TKT_FUNCTION,  //函数 

	TKT_ENDSIGN,   //结束符 
	TKT_UNKNOW     //未知符号 
}TokenType_Express;


class AFX_EX_CLASS CExpress
{
public:
	CExpress(const string &expression);
	CExpress();
	~CExpress();

public:
	int getVal(double &res, const string &expression);
	int getVal(double &res);

private:
	string expression;        //表达式 
	string token;             //每次读取的单词 
	TokenType_Express tkType; //读取的单词类型 
	int pos, length;          //读取的位置和长度 

	/* 静态数据 */
	static string ptrList[];   //支持的运算符
	static int ptrArgCnt[];    //运算符参数个数

	static string funList[];   //支持的函数列表 
	static int funArgCnt[];    //运算符参数个数 
	static int preceMap[][EXP_PTR_NUM];      //运算优先级表 

	/* 调试输出 */
	void debug();

	/* 读取下一个单词 */
	void readToken();

	/* 比较两个运算符的优先 */
	int comparePrece(const string &ptr1, const std::string &ptr2);

	/* 单步运算符计算 */
	double calculate(const string &ptr, double arg[]);

	/* 单步函数计算 */
	double callFun(const string &fun, double arg[]);

	/* 获取运算符序号 */
	int getPtrIndex(const string &ptr);

	/* 获取函数序号 */
	int getFunIndex(const string &fun);

	/* 从操作数栈获取n个参数 */
	bool getArg(stack<double> &opnd, double arg[], int n);
};

