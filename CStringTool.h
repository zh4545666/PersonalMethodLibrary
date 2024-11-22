#pragma once

#include <string>

#ifdef AFX_CLASS  
#define AFX_EX_CLASS _declspec(dllexport)  
#else  
#define AFX_EX_CLASS _declspec(dllimport)  
#endif

namespace PersonalMethod {

	class AFX_EX_CLASS CStringTool
	{
	public:
		CStringTool();
		~CStringTool();

	public:
		static std::string GbkToUtf8(const char* src_str);

		static std::string Utf8ToGbk(const char* src_str);
	};
}

