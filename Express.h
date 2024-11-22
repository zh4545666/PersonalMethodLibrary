#pragma once

#include <string>
#include <stack>

#ifdef AFX_CLASS  
#define AFX_EX_CLASS _declspec(dllexport)  
#else  
#define AFX_EX_CLASS _declspec(dllimport)  
#endif 

namespace PersonalMethod {

#define EXP_PTR_NUM 11  //运算符个数
#define EXP_FUN_NUM 12  //函数个数

	/* 单词符号类型 */
	typedef enum {
		TKT_NUMBER,    //数字 
		TKT_OPERATOR,  //运算符 
		TKT_FUNCTION,  //函数 

		TKT_ENDSIGN,   //结束符 
		TKT_UNKNOW     //未知符号 
	}TokenType_Express;

	class AFX_EX_CLASS CExpress
	{
	public:
		CExpress(const std::string& expression);
		CExpress();
		~CExpress();

	public:
		int getVal(double& res, const std::string& expression);
		int getVal(double& res);

	private:
		std::string expression;        //表达式 
		std::string token;             //每次读取的单词 
		TokenType_Express tkType; //读取的单词类型 
		int pos, length;          //读取的位置和长度 

		/* 静态数据 */
		static std::string ptrList[];   //支持的运算符
		static int ptrArgCnt[];    //运算符参数个数

		static std::string funList[];   //支持的函数列表 
		static int funArgCnt[];    //运算符参数个数 
		static int preceMap[][EXP_PTR_NUM];      //运算优先级表 

		/* 调试输出 */
		void debug();

		/* 读取下一个单词 */
		void readToken();

		/* 比较两个运算符的优先 */
		int comparePrece(const std::string& ptr1, const std::string& ptr2);

		/* 单步运算符计算 */
		double calculate(const std::string& ptr, double arg[]);

		/* 单步函数计算 */
		double callFun(const std::string& fun, double arg[]);

		/* 获取运算符序号 */
		int getPtrIndex(const std::string& ptr);

		/* 获取函数序号 */
		int getFunIndex(const std::string& fun);

		/* 从操作数栈获取n个参数 */
		bool getArg(std::stack<double>& opnd, double arg[], int n);
	};

}