///////////////////////////////////////////////////////////////////////////// 
// �ļ���: 
// ������: 
// ��������: 2009-09-27 ���� 04:51:46 
// 
// ˵��:ѹ����ѹ����ͼ�ļ��� 
/////////////////////////////////////////////////////////////////////////////
#include "StdAfx.h" 
#include "zipimplement.h" 
#include <direct.h>
#include <vector>
#include <xstring>
CZipImplement::CZipImplement(void) 
{ 
}
CZipImplement::~CZipImplement(void) 
{ 
}
///////////////////////////////////////////////////////////////////////////// 
// ����˵��: ʵ��ѹ���ļ��в��� 
// ����˵��: [in]�� pFilePath         Ҫ��ѹ�����ļ��� 
//                         mZipFileFullPath  ѹ������ļ�����·�� 
// ����ֵ: �������������·���FALSE,ѹ���ɹ��󷵻�TRUE 
// ��������: 
// ��������: 2009-09-27 ���� 04:58:52 
///////////////////////////////////////////////////////////////////////////// 
BOOL CZipImplement::Zip_PackFiles(CString& pFilePath, CString& mZipFileFullPath) 
{ 
    //�������� 
    if ((pFilePath == L"") || (mZipFileFullPath == L"")) 
    { 
        //·���쳣���� 
        return FALSE;
    }
    if(!CZipImplement::FolderExist(pFilePath)) 
    { 
        //Ҫ��ѹ�����ļ��в����� 
        return FALSE;
    }
    CString tZipFilePath = mZipFileFullPath.Left(mZipFileFullPath.ReverseFind('\\') + 1); 
//     if(!CZipImplement::FolderExist(tZipFilePath)) 
//     { 
//         //ZIP�ļ���ŵ��ļ��в����ڴ����� 
//         wchar_t* temp=tZipFilePath.GetBuffer(tZipFilePath.GetLength()); 
//         if (FALSE == CreatedMultipleDirectory(temp)) 
//         { 
//             //����Ŀ¼ʧ�� 
//             return FALSE; 
//         } 
//     }
    //����ļ��е����� 
    if(pFilePath.Right(1) == L"\\") 
    { 
        this->m_FolderPath = pFilePath.Left(pFilePath.GetLength() - 1); 
        m_FolderName = m_FolderPath.Right(m_FolderPath.GetLength() - m_FolderPath.ReverseFind('\\') - 1); 
    } 
    else 
    { 
        this->m_FolderPath = pFilePath; 
        m_FolderName = pFilePath.Right(pFilePath.GetLength() - pFilePath.ReverseFind('\\') - 1); 
    }
    /************************************************************************/
    //����ZIP�ļ� 
    hz=CreateZip(mZipFileFullPath,0); 
    if(hz == 0) 
    { 
        //����Zip�ļ�ʧ�� 
        return FALSE; 
    }
    //�ݹ��ļ���,����ȡ���ʼۼ���ZIP�ļ� 
    BrowseFile(pFilePath); 
    //�ر�ZIP�ļ� 
    CloseZip(hz);
    /************************************************************************/
    CFileFind tFFind; 
    if (!tFFind.FindFile(mZipFileFullPath)) 
    { 
        //ѹ��ʧ��(δ����ѹ������ļ�) 
        return FALSE; 
    }
    return TRUE; 
}
///////////////////////////////////////////////////////////////////////////// 
// ����˵��: ��ѹ���ļ��� 
// ����˵��: [in]�� mUnPackPath     ��ѹ���ļ���ŵ�·�� 
//                         mZipFileFullPath  ZIP�ļ���·�� 
// ����ֵ: 
// ��������: 
// ��������: 2009-09-27 ���� 11:04:28 
///////////////////////////////////////////////////////////////////////////// 
BOOL CZipImplement::Zip_UnPackFiles(CString &mZipFileFullPath, CString& mUnPackPath) 
{ 
    //�������� 
    if ((mUnPackPath == L"") || (mZipFileFullPath == L"")) 
    { 
        //·���쳣���� 
        return FALSE ; 
    }
    CFileFind tFFind; 
    if (!tFFind.FindFile(mZipFileFullPath)) 
    { 
        //ѹ��ʧ��(δ����ѹ���ļ�) 
        return FALSE; 
    }
    //�����ѹ����·�������� ��ͼ������ 
    CString tZipFilePath = mUnPackPath; 
//     if(!CZipImplement::FolderExist(tZipFilePath)) 
//     { 
//         //��ѹ���ŵ��ļ��в����� ������ 
//         wchar_t* temp=tZipFilePath.GetBuffer(tZipFilePath.GetLength()); 
//         if (FALSE == CreatedMultipleDirectory(temp)) 
//         { 
//             //����Ŀ¼ʧ�� 
//             return FALSE; 
//         } 
//     } 
    /************************************************************************/ 
    //��ZIP�ļ� 
    hz=OpenZip(mZipFileFullPath,0); 
    if(hz == 0) 
    { 
        //��Zip�ļ�ʧ�� 
        return FALSE; 
    }
    zr=SetUnzipBaseDir(hz,mUnPackPath); 
    if(zr != ZR_OK) 
    { 
        //��Zip�ļ�ʧ�� 
        CloseZip(hz); 
        return FALSE;       
    }
    zr=GetZipItem(hz,-1,&ze); 
    if(zr != ZR_OK) 
    { 
        //��ȡZip�ļ�����ʧ�� 
        CloseZip(hz); 
        return FALSE;       
    }
    int numitems=ze.index; 
    for (int i=0; i<numitems; i++)   
	{ 
        zr=GetZipItem(hz,i,&ze); 
        zr=UnzipItem(hz,i,ze.name);
        if(zr != ZR_OK) 
        { 
            //��ȡZip�ļ�����ʧ�� 
            CloseZip(hz); 
            return FALSE;       
        } 
    }
    CloseZip(hz); 
    return TRUE; 
}
///////////////////////////////////////////////////////////////////////////// 
// ����˵��: ���ָ�����ļ����Ƿ���� 
// ����˵��: [in]��strPath �����ļ��� (�˷�����������·��ĩβ���*.*) 
// ����ֵ:BOOL����,���ڷ���TRUE,����ΪFALSE 
// ��������: 
// ��������: 2009-09-27 ���� 02:16:36 
///////////////////////////////////////////////////////////////////////////// 
BOOL CZipImplement::FolderExist(CString& strPath) 
{ 
    CString sCheckPath = strPath;
    if(sCheckPath.Right(1) != L"\\") 
        sCheckPath += L"\\";
    sCheckPath += L"*.*";
    WIN32_FIND_DATA wfd; 
    BOOL rValue = FALSE;
    HANDLE hFind = FindFirstFile(sCheckPath, &wfd);
    if ((hFind!=INVALID_HANDLE_VALUE) && 
        (wfd.dwFileAttributes&FILE_ATTRIBUTE_DIRECTORY) || (wfd.dwFileAttributes&FILE_ATTRIBUTE_ARCHIVE)) 
    { 
        //������ڲ��������ļ��� 
        rValue = TRUE; 
    }
    FindClose(hFind); 
    return rValue; 
}
///////////////////////////////////////////////////////////////////////////// 
// ����˵��: �����ļ��� 
// ����˵��: [in]��strFile �������ļ���(�˷�����������·��ĩβ���*.*) 
// ����ֵ:BOOL����,���ڷ���TRUE,����ΪFALSE 
// ��������: 
// ��������: 2009-09-27 ���� 02:16:36 
///////////////////////////////////////////////////////////////////////////// 
void CZipImplement::BrowseFile(CString &strFile) 
{ 
    CFileFind ff; 
    CString szDir = strFile;
    if(szDir.Right(1) != L"\\") 
        szDir += L"\\";
    szDir += L"*.*";
    BOOL res = ff.FindFile(szDir); 
    while(res) 
    { 
        res = ff.FindNextFile(); 
        if(ff.IsDirectory() && !ff.IsDots()) 
        { 
            //�����һ����Ŀ¼���õݹ��������һ���� 
            CString strPath = ff.GetFilePath();
            CString subPath; 
            GetRelativePath(strPath,subPath); 
            //���ļ���ӵ�ZIP�ļ� 
            ZipAddFolder(hz,subPath); 
            BrowseFile(strPath); 
        } 
        else if(!ff.IsDirectory() && !ff.IsDots()) 
        { 
            //��ʾ��ǰ���ʵ��ļ�(����·��) 
            CString strPath = ff.GetFilePath(); 
            CString subPath;
            GetRelativePath(strPath,subPath); 
            //���ļ���ӵ�ZIP�ļ� 
            ZipAdd(hz,subPath,strPath); 
        } 
    }
    //�ر� 
    ff.Close(); 
}
///////////////////////////////////////////////////////////////////////////// 
// ����˵��: ��ȡ���·�� 
// ����˵��: [in]��pFullPath ��ǰ�ļ�������·�� [out] : ����������·�� 
// ��������: 
// ��������: 2009-9-28 ���� 11:17:21 
///////////////////////////////////////////////////////////////////////////// 
void CZipImplement::GetRelativePath(CString& pFullPath,CString& pSubString) 
{ 
    pSubString = pFullPath.Right(pFullPath.GetLength() - this->m_FolderPath.GetLength() + this->m_FolderName.GetLength()); 
}
/*
///////////////////////////////////////////////////////////////////////////// 
// ����˵��: �����༶Ŀ¼ 
// ����˵��: [in]�� ·���ַ��� 
// ����ֵ: BOOL �ɹ�True ʧ��False 
// ��������: 
// ��������: 2009-9-28 ���� 04:53:20 
///////////////////////////////////////////////////////////////////////////// 
BOOL CZipImplement::CreatedMultipleDirectory(wchar_t* direct) 
{ 
    std::wstring Directoryname = direct;
    if (Directoryname[Directoryname.length() - 1] !=  '\\') 
    { 
        Directoryname.append(1, '\\'); 
    } 
    std::vector< std::wstring> vpath; 
    std::wstring strtemp; 
    BOOL  bSuccess = FALSE; 
    for (int i = 0; i < Directoryname.length(); i++) 
    { 
        if ( Directoryname[i] != '\\') 
        { 
            strtemp.append(1,Directoryname[i]);    
        } 
        else 
        { 
            vpath.push_back(strtemp); 
            strtemp.append(1, '\\'); 
        } 
    } 
    std::vector<std::wstring>:: const_iterator vIter; 
    for (vIter = vpath.begin();vIter != vpath.end(); vIter++) 
    { 
        bSuccess = CreateDirectory(vIter->c_str(), NULL) ? TRUE :FALSE; 
    }
    return bSuccess; 
}
*/