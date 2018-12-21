///////////////////////////////////////////////////////////////////////////// 
// �ļ���: 
// ������: 
// ��������: 2009-09-27 ���� 04:51:46 
// 
// ˵��:ѹ����ѹ����ͼ�ļ��� 
/////////////////////////////////////////////////////////////////////////////
#pragma once
#include "zip.h" 
#include "unzip.h"
class CZipImplement 
{ 
public: 
    CZipImplement(void); 
    ~CZipImplement(void);
private: 
    HZIP hz;          //Zip�ļ���� 
    ZRESULT zr;    //��������ֵ 
    ZIPENTRY ze;  //Zip�ļ����
    CString m_FolderPath;     //folder·�� 
    CString  m_FolderName;  //folder��Ҫ��ѹ�����ļ�����
private: 
    //ʵ�ֱ����ļ��� 
    void BrowseFile(CString &strFile);
    //��ȡ���·�� 
    void GetRelativePath(CString& pFullPath, CString& pSubString);
    //����·�� 
//    BOOL CreatedMultipleDirectory(wchar_t* direct); 
    ///* 
    //*********************************************************************** 
    //* ������ TransCStringToTCHAR 
    //* ��������CString ת��Ϊ TCHAR* 
    //*********************************************************************** 
    //*/ 
    //TCHAR* CString2TCHAR(CString &str) 
    //{ 
    //    int iLen = str.GetLength(); 
    //    TCHAR* szRs = new TCHAR[iLen]; 
    //    lstrcpy(szRs, str.GetBuffer(iLen)); 
    //    str.ReleaseBuffer(); 
    //    return szRs; 
    //}
public: 
    //ѹ���ļ��нӿ� 
    BOOL Zip_PackFiles(CString& pFilePath, CString& mZipFileFullPath);
    //��ѹ���ļ��нӿ� 
    BOOL Zip_UnPackFiles(CString &mZipFileFullPath, CString& mUnPackPath);
public: 
    //��̬�����ṩ�ļ���·����� 
    static BOOL FolderExist(CString& strPath); 
};