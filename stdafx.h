// stdafx.h : ��׼ϵͳ�����ļ��İ����ļ���
// ���Ǿ���ʹ�õ��������ĵ�
// �ض�����Ŀ�İ����ļ�
//

#pragma once

#include "targetver.h"

#define WIN32_LEAN_AND_MEAN             //  �� Windows ͷ�ļ����ų�����ʹ�õ���Ϣ
// Windows ͷ�ļ�: 
#include <windows.h>

//�궨��
#define AFX_CLASS

#define InitRandom srand((unsigned)time(NULL))
#define GetRandom( min, max ) ((rand() % (int)(((max)+1) - (min))) + (min))

// TODO:  �ڴ˴����ó�����Ҫ������ͷ�ļ�
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

