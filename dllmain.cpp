// dllmain.cpp : ���� DLL Ӧ�ó������ڵ㡣
#include "stdafx.h"


bool APIENTRY DllMain( HMODULE hModule,
                       DWORD  ul_reason_for_call,
                       LPVOID lpReserved
					 )
{
	switch (ul_reason_for_call)
	{
	case DLL_PROCESS_ATTACH:
	{
		bool bflag = false;
		wchar_t path[MAX_PATH];
		//NULLΪ���ó����ַ�� hModuleΪ����̬�����
		GetModuleFileName(hModule, path, MAX_PATH);
		wstring strPath(path);

		int nPos = strPath.rfind(L"\\");
		if (nPos > 0)
		{
			strPath = strPath.substr(0, nPos+1) + L"LICENSE";
			ifstream in;
			in.open(strPath, ios::in | ios::binary);
			if (in.is_open())
			{
				bflag = true;
			}			
		}

		if (!bflag)
			assert(1 == 2);
	}break;
	case DLL_THREAD_ATTACH:
	case DLL_THREAD_DETACH:
	case DLL_PROCESS_DETACH:
		break;
	}
	return true;
}
