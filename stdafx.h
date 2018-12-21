// stdafx.h : 标准系统包含文件的包含文件，
// 或是经常使用但不常更改的
// 特定于项目的包含文件
//

#pragma once

#include "targetver.h"

#define WIN32_LEAN_AND_MEAN             //  从 Windows 头文件中排除极少使用的信息
// Windows 头文件: 
#include <windows.h>

//宏定义
#define AFX_CLASS

#define InitRandom srand((unsigned)time(NULL))
#define GetRandom( min, max ) ((rand() % (int)(((max)+1) - (min))) + (min))

// TODO:  在此处引用程序需要的其他头文件
#include <stdio.h>
#include <stdlib.h>

#include <string>
#include <stack>
#include <vector>
#include <sstream>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <iomanip> 
using namespace std;

#include "math.h"
#include "assert.h"

template class _declspec(dllexport) std::basic_string<char, std::char_traits<char>, std::allocator<char>>;
//template class _declspec(dllexport) std::allocator<double>;
//template class _declspec(dllexport) std::vector<double, std::allocator<double> >;

