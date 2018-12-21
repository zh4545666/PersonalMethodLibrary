// PersonalMethodLibrary.cpp : ���� DLL Ӧ�ó���ĵ���������
//

#include "stdafx.h"

/**************************************************
*
* Description:��ģ��˵��
*
**************************************************/
//AES AES����
//MD5 MD5ֵ����
//Fit ���
//Complex ����
//Integral ��ֵ����
//Interpolate ��ֵ
//Matrix ����
//LEquations ������Է�����
//NLEquations �������Է�����
//IntegralTransform ���ֱ任

//DataAnalysis ��ֵͳ�Ʒ���
//SortAlgorithm �����㷨
//StringBuffer �ַ�������
//Express ��ʽ����

#ifdef AFX_CLASS  
#define AFX_EX_CLASS _declspec(dllexport)  
#else  
#define AFX_EX_CLASS _declspec(dllimport)  
#endif 


/**************************************************
*
*function�����DLL����
*
**************************************************/
extern "C" AFX_EX_CLASS bool GetDllName(char * src,size_t len)
{
	char str[] = { "HZA-2018-1.0" };
	if (len < strlen(str) || !src) 
		return false;
	memset(src, 0, sizeof(char)*len);
	strcpy_s(src, strlen(str), str);
	return true;
}

/**************************************************
*
*function��N�ʺ�����-����㷨1
*
**************************************************/
//#include<iostream>
//#define N 8
//using namespace std;
//static int gEightQueen[N] = { 0 }, gCount = 0;
//void print()//���ÿһ������������лʺ�İڷ����
//{
//	cout << gCount <<endl;
//	return;
//	for (int i = 0; i < N; i++)
//	{
//		int inner;
//		for (inner = 0; inner < gEightQueen[i]; inner++)
//			cout << "0 ";
//		cout << "# ";
//		for (inner = gEightQueen[i] + 1; inner < N; inner++)
//			cout << "0 ";
//		cout << endl;
//	}
//	cout << "==========================\n";
//}
//int check_pos_valid(int loop, int value)//����Ƿ�����ж���ʺ���ͬһ��/��/�Խ��ߵ����
//{
//	int index;
//	int data;
//	for (index = 0; index < loop; index++)
//	{
//		data = gEightQueen[index];
//		if (value == data)
//			return 0;
//		if ((index + data) == (loop + value))
//			return 0;
//		if ((index - data) == (loop - value))
//			return 0;
//	}
//	return 1;
//}
//void eight_queen(int index)
//{
//	int loop;
//	for (loop = 0; loop < N; loop++)
//	{
//		if (check_pos_valid(index, loop))
//		{
//			gEightQueen[index] = loop;
//			if (N-1 == index)
//			{
//				gCount++, print();
//				gEightQueen[index] = 0;
//				return;
//			}
//			eight_queen(index + 1);
//			gEightQueen[index] = 0;
//		}
//	}
//}
//int main(int argc, char*argv[])
//{
//	eight_queen(0);
//	cout << "total=" << gCount << endl;
//
//	system("pause");
//	return 0;
//}


/**************************************************
*
*function��N�ʺ�����-����㷨2
* ��̽-�����㷨���ݹ�ʵ��
*
**************************************************/
//#include "iostream"  
//using namespace std;
//#include "time.h"  
//
//// sum������¼�ʺ���óɹ��Ĳ�ͬ��������upperlim������������ж��Ѿ����ú��˻ʺ�  
//long sum = 0, upperlim = 1;
//
//// ��̽�㷨�����ұߵ��п�ʼ��  
//void test(long row, long ld, long rd)
//{
//	if (row != upperlim)
//	{
//		// row��ld��rd���С������㣬������п��Է��ûʺ����,��ӦλΪ0��  
//		// Ȼ����ȡ�����롱��ȫ1����������õ�ǰ���п��Է��ûʺ��λ�ã���Ӧ�и�Ϊ1  
//		// Ҳ������ȡ��ǰ��Щ�п��Է��ûʺ�  
//		long pos = upperlim & ~(row | ld | rd);
//		while (pos)    // 0 -- �ʺ�û�еط��ɷţ�����  
//		{
//			// ����pos���ұ�Ϊ1��bit������bit��0  
//			// Ҳ����ȡ�ÿ��ԷŻʺ�����ұߵ���  
//			long p = pos & -pos;
//
//			// ��pos���ұ�Ϊ1��bit����  
//			// Ҳ����Ϊ��ȡ��һ�ε����ҿ�����ʹ����׼����  
//			// ����������ݵ����λ�ü�����̽  
//			pos -= p;
//
//			// row + p������ǰ����1����ʾ��¼��λʺ���õ��С�  
//			// (ld + p) << 1����ǵ�ǰ�ʺ�������ڵ��в�������һ���ʺ���á�  
//			// (ld + p) >> 1����ǵ�ǰ�ʺ��ұ����ڵ��в�������һ���ʺ���á�  
//			// �˴�����λ����ʵ�����Ǽ�¼�Խ����ϵ����ƣ�ֻ����Ϊ���ⶼ����  
//			// ��һ������������������Ա�ʾΪ�е����ƾͿ����ˡ���Ȼ��������λ  
//			// ��ÿ��ѡ����֮ǰ���У�ԭ��N��N������ĳ���ѷ��õĻʺ������Խ���  
//			// �ϲ��������ƶ�����¼������  
//			test(row + p, (ld + p) << 1, (rd + p) >> 1);
//		}
//	}
//	else
//	{
//		// row������λ��Ϊ1�����ҵ���һ���ɹ��Ĳ��֣�����  
//		sum++;
//	}
//}
//
//int main(int argc, char *argv[])
//{
//	time_t tm;
//	int n = 16;
//	cout << "����ʺ�����";
//	cin >> n;
//	cout << endl;
//
//	if (argc != 1)
//		n = atoi(argv[1]);
//	tm = time(0);
//
//	// ��Ϊ�����������ƣ����ֻ��32λ��  
//	// ����봦��N����32�Ļʺ����⣬��Ҫ  
//	// ��bitset���ݽṹ���д洢  
//	if ((n < 1) || (n > 32))
//	{
//		printf(" ֻ�ܼ���1-32֮��\n");
//		exit(-1);
//	}
//	printf("%d �ʺ�\n", n);
//
//	// N���ʺ�ֻ��Nλ�洢��N����ĳ���лʺ����Ӧbit��1��  
//	upperlim = (upperlim << n) - 1;
//
//	test(0, 0, 0);
//	printf("����%ld������, ����ʱ��%d�� \n", sum, (int)(time(0) - tm));
//	system("pause");
//	return 0;
//}


/**************************************************
*
*function�������ⷨ
*
**************************************************/
//bool sign = false;
//
// //������������ 
//int num[9][9];
//
// //�������� 
//void Input();
//void Output();
//bool Check(int n, int key);
//void DFS(int n);
//
// //������ 
//int main()
//{
//	cout << "������һ��9*9���������󣬿�λ��0��ʾ:" << endl;
//	Input();
//	DFS(0);
//	Output();
//	system("pause");
//}
//
// //������������ 
//void Input()
//{
//	char temp[9][9];
//	for (int i = 0; i < 9; i++)
//	{
//		for (int j = 0; j < 9; j++)
//		{
//			cin >> temp[i][j];
//			num[i][j] = temp[i][j] - '0';
//		}
//	}
//}
//
// //����������� 
//void Output()
//{
//	cout << endl;
//	for (int i = 0; i < 9; i++)
//	{
//		for (int j = 0; j < 9; j++)
//		{
//			cout << num[i][j] << " ";
//			if (j % 3 == 2)
//			{
//				cout << "   ";
//			}
//		}
//		cout << endl;
//		if (i % 3 == 2)
//		{
//			cout << endl;
//		}
//	}
//}
//
// //�ж�key����nʱ�Ƿ��������� 
//bool Check(int n, int key)
//{
//	 //�ж�n���ں����Ƿ�Ϸ� 
//	for (int i = 0; i < 9; i++)
//	{
//		 //jΪn������ 
//		int j = n / 9;
//		if (num[j][i] == key) return false;
//	}
//
//	 //�ж�n���������Ƿ�Ϸ� 
//	for (int i = 0; i < 9; i++)
//	{
//		 //jΪn������ 
//		int j = n % 9;
//		if (num[i][j] == key) return false;
//	}
//
//	 //xΪn���ڵ�С�Ź����󶥵������� 
//	int x = n / 9 / 3 * 3;
//
//	 //yΪn���ڵ�С�Ź����󶥵������ 
//	int y = n % 9 / 3 * 3;
//
//	 //�ж�n���ڵ�С�Ź����Ƿ�Ϸ� 
//	for (int i = x; i < x + 3; i++)
//	{
//		for (int j = y; j < y + 3; j++)
//		{
//			if (num[i][j] == key) return false;
//		}
//	}
//
//	 //ȫ���Ϸ���������ȷ 
//	return true;
//}
//
// //���ѹ������� 
//void DFS(int n)
//{
//	 //���еĶ����ϣ��˳��ݹ� 
//	if (n > 80)
//	{
//		sign = true;
//		return;
//	}
//	 //��ǰλ��Ϊ��ʱ���� 
//	if (num[n / 9][n % 9] != 0)
//	{
//		DFS(n + 1);
//	}
//	else
//	{
//		 //����Ե�ǰλ����ö�ٲ��� 
//		for (int i = 1; i <= 9; i++)
//		{
//			 //��������ʱ�������� 
//			if (Check(n, i) == true)
//			{
//				num[n / 9][n % 9] = i;
//				 //�������� 
//				DFS(n + 1);
//				 //����ʱ�������ɹ�����ֱ���˳� 
//				if (sign == true) return;
//				 //������첻�ɹ�����ԭ��ǰλ 
//				num[n / 9][n % 9] = 0;
//			}
//		}
//	}
//}


/**************************************************
*
*function����ɶ˿�ʵ��
*
**************************************************/
//#include <Winsock2.h> 
//#pragma comment(lib, "ws2_32.lib ")
//#include <windows.h>
//
//#include <string>
//#include <list>
//using namespace std;
//
//#define DefaultPort 20000  
//#define DataBuffSize 8 * 1024
//
////�������ݴ����
//class CSocketData
//{
//public:
//
//	CSocketData()
//	{
//		m_str = "";
//	}
//	string m_str;
//
//
//	void add(const char* pch)
//	{
//		m_str += pch;
//	}
//
//	void remove(size_t len)
//	{
//		if (len >= m_str.length())
//			m_str = "";
//		else
//		{
//			for (int i = len; i >= 0; i--)
//				m_str.erase(m_str.begin() + i);
//		}
//	}
//protected:
//private:
//};
//typedef CSocketData* PSOCKETDATA;
////��ɶ˿������ж�
//enum TYPE_SOCKET
//{
//	SEND_SOCKET, RECV_SOCKET
//};
//
//typedef struct
//{
//	OVERLAPPED overlapped;
//	CHAR buffer[DataBuffSize];
//	TYPE_SOCKET Type;
//	PSOCKETDATA pSocketData;
//}PER_IO_OPERATEION_DATA, *LPPER_IO_OPERATION_DATA;
//
//typedef struct
//{
//	SOCKET socket;
//}PER_HANDLE_DATA, *LPPER_HANDLE_DATA;
//
//DWORD WINAPI ServerWorkThread(LPVOID CompletionPortID);
//DWORD WINAPI RecvProc(LPVOID CompletionPortID);
//
//typedef struct SOURCEDATA
//{
//	LPPER_IO_OPERATION_DATA pIoData;
//	LPPER_HANDLE_DATA pHandleData;
//
//	bool operator ==(const SOURCEDATA& other) const
//	{
//		return ((this->pIoData == other.pIoData) && (this->pHandleData == other.pHandleData));
//	}
//}IOCPDATA;
//typedef list<IOCPDATA> IOCPLIST;
//
////��Դ����(��֤�رճ���ʱ��������Դ�����ͷ�)
//IOCPLIST IOCPList;
////������
//HANDLE hmutex;
//
//void main()
//{
//	//����������
//	hmutex = CreateMutex(NULL, TRUE, L"ICOPLIST");
//	ReleaseMutex(hmutex);
//
//	//��ɶ˿ھ��
//	HANDLE completionPort;
//
//	//��ʼ���׽���
//	WSADATA wsaData;
//	DWORD ret;
//	if (ret = WSAStartup(0x0202, &wsaData) != 0)
//	{
//		std::cout << "WSAStartup failed. Error:" << ret << std::endl;
//		return;
//	}
//	//������ɶ˿�
//	completionPort = CreateIoCompletionPort(INVALID_HANDLE_VALUE, NULL, 0, 0);
//	if (completionPort == NULL)
//	{
//		std::cout << "CreateIoCompletionPort failed. Error:" << GetLastError() << std::endl;
//		return;
//	}
//	//���ϵͳ��Ϣ
//	SYSTEM_INFO mySysInfo;
//	GetSystemInfo(&mySysInfo);
//	// ���� 2 * CPU���� + 1 ���߳�  
//	DWORD threadID;
//	for (DWORD i = 0; i < (mySysInfo.dwNumberOfProcessors * 2 + 1); ++i)
//	{
//		HANDLE threadHandle;
//		threadHandle = CreateThread(NULL, 0, ServerWorkThread, completionPort, 0, &threadID);
//		if (threadHandle == NULL)
//		{
//			std::cout << "CreateThread failed. Error:" << GetLastError() << std::endl;
//			return;
//		}
//
//		CloseHandle(threadHandle);
//	}
//	//���������߳�
//	{
//		HANDLE threadHandle = CreateThread(NULL, 0, RecvProc, completionPort, 0, &threadID);
//		if (threadHandle == NULL)
//		{
//			std::cout << "CreateThread failed. Error:" << GetLastError() << std::endl;
//			if (completionPort)
//				CloseHandle(completionPort);
//			return;
//		}
//		CloseHandle(threadHandle);
//	}
//
//	system("pause");
//
//	//�ر���ɶ˿ڣ����˳������߳�
//	if (completionPort)
//		CloseHandle(completionPort);
//	//ȷ��������Դ�ͷ�
//	for (IOCPLIST::iterator iter = IOCPList.begin(); iter != IOCPList.end();)
//	{
//		delete iter->pIoData->pSocketData;
//		GlobalFree(iter->pIoData);
//		if (iter->pHandleData)
//			closesocket(iter->pHandleData->socket);
//		GlobalFree(iter->pHandleData);
//		iter = IOCPList.erase(iter);
//	}
//
//}
//
//DWORD WINAPI RecvProc(LPVOID CompletionPortID)
//{
//	HANDLE complationPort = (HANDLE)CompletionPortID;
//
//	SOCKET acceptSocket;
//	LPPER_HANDLE_DATA pHandleData;
//	LPPER_IO_OPERATION_DATA pIoData;
//	DWORD recvBytes;
//	DWORD flags;
//
//	// ����һ������socket  
//	SOCKET listenSocket = WSASocket(AF_INET, SOCK_STREAM, 0, NULL, 0, WSA_FLAG_OVERLAPPED);
//	if (listenSocket == INVALID_SOCKET)
//	{
//		std::cout << " WSASocket( listenSocket ) failed. Error:" << GetLastError() << std::endl;
//		return 0;
//	}
//
//	SOCKADDR_IN internetAddr;
//	internetAddr.sin_family = AF_INET;
//	internetAddr.sin_addr.s_addr = htonl(INADDR_ANY);
//	internetAddr.sin_port = htons(DefaultPort);
//
//	// �󶨼����˿�  
//	if (bind(listenSocket, (PSOCKADDR)&internetAddr, sizeof(internetAddr)) == SOCKET_ERROR)
//	{
//		std::cout << "Bind failed. Error:" << GetLastError() << std::endl;
//		return 0;
//	}
//
//	if (listen(listenSocket, 5) == SOCKET_ERROR)
//	{
//		std::cout << "listen failed. Error:" << GetLastError() << std::endl;
//		return 0;
//	}
//
//	// ��ʼ��ѭ������������  
//	while (1)
//	{
//		sockaddr_in  cliaddr;
//		int len = sizeof(cliaddr);
//
//		acceptSocket = WSAAccept(listenSocket, (sockaddr*)&cliaddr, &len, NULL, 0);
//		if (acceptSocket == SOCKET_ERROR)
//		{
//			std::cout << "WSAAccept failed. Error:" << GetLastError() << std::endl;
//			continue;
//		}
//		//�ͻ���IP��ַ������ָ���м�·������ַ��
//		cout << inet_ntoa(cliaddr.sin_addr) << endl;
//		//�ͻ���IP�˿�
//		cout << cliaddr.sin_port << endl;
//
//		pHandleData = (LPPER_HANDLE_DATA)GlobalAlloc(GPTR, sizeof(PER_HANDLE_DATA));
//		if (pHandleData == NULL)
//		{
//			std::cout << "GlobalAlloc( HandleData ) failed. Error:" << GetLastError() << std::endl;
//			continue;
//		}
//		//���׽��ֺ���ɶ˿ڰ�
//		pHandleData->socket = acceptSocket;
//		if (CreateIoCompletionPort((HANDLE)acceptSocket, complationPort, (DWORD)pHandleData, 0) == NULL)
//		{
//			std::cout << "CreateIoCompletionPort failed. Error:" << GetLastError() << std::endl;
//			continue;
//		}
//
//		pIoData = (LPPER_IO_OPERATION_DATA)GlobalAlloc(GPTR, sizeof(PER_IO_OPERATEION_DATA));
//		if (pIoData == NULL)
//		{
//			std::cout << "GlobalAlloc( IoData ) failed. Error:" << GetLastError() << std::endl;
//			continue;
//		}
//
//		ZeroMemory(&(pIoData->overlapped), sizeof(pIoData->overlapped));
//		ZeroMemory(pIoData->buffer, sizeof(CHAR)*DataBuffSize);
//		pIoData->Type = RECV_SOCKET;
//		pIoData->pSocketData = new CSocketData;
//		WSABUF buf;
//		buf.len = DataBuffSize;
//		buf.buf = pIoData->buffer;
//
//		flags = 0;
//		if (WSARecv(acceptSocket, &buf, 1, &recvBytes, &flags, &(pIoData->overlapped), NULL) == SOCKET_ERROR)
//		{
//			if (WSAGetLastError() != ERROR_IO_PENDING)
//			{
//				std::cout << "WSARecv() failed. Error:" << GetLastError() << std::endl;
//				continue;
//			}
//			else
//			{
//// 				std::cout << "WSARecv() io pending" << std::endl;  
//// 				return;
//			}
//		}
//
//		IOCPDATA Info;
//		Info.pHandleData = pHandleData;
//		Info.pIoData = pIoData;
//		WaitForSingleObject(hmutex, INFINITE);
//		IOCPList.push_back(Info);
//		ReleaseMutex(hmutex);
//	}
//}
//
//DWORD WINAPI ServerWorkThread(LPVOID CompletionPortID)
//{
//	HANDLE complationPort = (HANDLE)CompletionPortID;
//
//	DWORD bytesTransferred;
//	LPPER_HANDLE_DATA pHandleData = NULL;
//	LPPER_IO_OPERATION_DATA pIoData = NULL;
//	DWORD sendBytes = 0;
//	DWORD recvBytes = 0;
//	DWORD flags = 0;
//
//	while (1)
//	{
//		if (GetQueuedCompletionStatus(complationPort, &bytesTransferred, (LPDWORD)&pHandleData, (LPOVERLAPPED *)&pIoData, INFINITE) == 0)
//		{
//			std::cout << "GetQueuedCompletionStatus failed. Error:" << GetLastError() << std::endl;
//			break;
//		}
//
//		// ��������Ƿ��Ѿ���������  
//		if (bytesTransferred == 0)
//		{
//			//�˳��߳�
//			if (!pHandleData && !pIoData)
//				break;
//
//			//�ͻ��˶Ͽ�	//����б�����ж���
//			IOCPDATA Info;
//			Info.pHandleData = pHandleData;
//			Info.pIoData = pIoData;
//			WaitForSingleObject(hmutex, INFINITE);
//			IOCPLIST::iterator iter = find(IOCPList.begin(), IOCPList.end(), Info);
//			if (iter != IOCPList.end())
//				IOCPList.erase(iter);
//			ReleaseMutex(hmutex);
//
//			//�ͷ���Դ
//			std::cout << " Start closing socket..." << std::endl;
//			if (pHandleData)
//			{
//				if (CloseHandle((HANDLE)pHandleData->socket) == SOCKET_ERROR)
//				{
//					std::cout << "Close socket failed. Error:" << GetLastError() << std::endl;
//					continue;
//				}
//				GlobalFree(pHandleData);
//			}
//			if (pIoData)
//			{
//				delete pIoData->pSocketData;
//				GlobalFree(pIoData);
//			}
//			continue;
//		}
//
//		if (RECV_SOCKET == pIoData->Type)
//		{
//			if (bytesTransferred == DataBuffSize)
//			{
//				//����δ�����ꡢ�����ѽ�������
//				if (pIoData->pSocketData)
//				{
//					pIoData->pSocketData->add(pIoData->buffer);
//				}
//
//				//����ص��ṹ
//				ZeroMemory(&(pIoData->overlapped), sizeof(OVERLAPPED));
//				//����Ͷ�Ž���I/O����
//				ZeroMemory(pIoData->buffer, sizeof(CHAR)*DataBuffSize);
//				WSABUF buf;
//				buf.len = DataBuffSize;
//				buf.buf = pIoData->buffer;
//				if (WSARecv(pHandleData->socket, &buf, 1, &recvBytes, &flags, &(pIoData->overlapped), NULL) == SOCKET_ERROR)
//				{
//					if (WSAGetLastError() != ERROR_IO_PENDING)
//					{
//						std::cout << "WSARecv() failed. Error:" << GetLastError() << std::endl;
//					}
//				}
//			}
//			else
//			{
//				if (pIoData->pSocketData)
//				{
//					pIoData->pSocketData->add(pIoData->buffer);
//				}
//
//				//����ص��ṹ			
//				ZeroMemory(&(pIoData->overlapped), sizeof(OVERLAPPED));
//				//������յ�������
//				string str = "Send : \t" + pIoData->pSocketData->m_str + "\r\n";
//				std::cout << str;
//				//Ͷ�ŷ���I/O����
//				pIoData->Type = SEND_SOCKET;
//				WSABUF buf;
//				buf.len = pIoData->pSocketData->m_str.length();
//				buf.buf = const_cast<char*>(pIoData->pSocketData->m_str.c_str());
//				if (WSASend(pHandleData->socket, &buf, 1, &sendBytes, 0, &(pIoData->overlapped), NULL) == SOCKET_ERROR)
//				{
//					if (WSAGetLastError() != ERROR_IO_PENDING)
//						std::cout << "WSASend() failed. Error:" << GetLastError() << std::endl;
//				}
//			}
//
//		}
//		else if (SEND_SOCKET == pIoData->Type)
//		{
//			if (bytesTransferred != pIoData->pSocketData->m_str.length())
//			{
//				//����δ�����꣬����Ͷ�ŷ���I/O����
//				pIoData->pSocketData->remove(bytesTransferred);
//
//				WSABUF buf;
//				buf.len = pIoData->pSocketData->m_str.length();
//				buf.buf = const_cast<char*>(pIoData->pSocketData->m_str.c_str());
//				if (WSASend(pHandleData->socket, &buf, 1, &sendBytes, 0, &(pIoData->overlapped), NULL) == SOCKET_ERROR)
//				{
//					if (WSAGetLastError() != ERROR_IO_PENDING)
//						std::cout << "WSASend() failed. Error:" << GetLastError() << std::endl;
//				}
//			}
//			else
//			{
//				pIoData->pSocketData->m_str = "";
//
//				//����ص��ṹ			
//				ZeroMemory(&(pIoData->overlapped), sizeof(OVERLAPPED));
//				//Ͷ�Ž���I/O����
//				ZeroMemory(pIoData->buffer, sizeof(CHAR)*DataBuffSize);
//				pIoData->Type = RECV_SOCKET;
//				WSABUF buf;
//				buf.len = DataBuffSize;
//				buf.buf = pIoData->buffer;
//				if (WSARecv(pHandleData->socket, &buf, 1, &recvBytes, &flags, &(pIoData->overlapped), NULL) == SOCKET_ERROR)
//				{
//					if (WSAGetLastError() != ERROR_IO_PENDING)
//					{
//						std::cout << "WSARecv() failed. Error:" << GetLastError() << std::endl;
//					}
//				}
//			}
//		}
//	}
//	return 0;
//}


/**************************************************
*
*function��ʮ�����Ƹ�ʽ���
*
**************************************************/
//void printHex(std::ifstream& ifs, std::ostream& ostream){
//	using namespace std;
//	ostream << setfill('0') << hex << uppercase;
//
//	unsigned char byte;
//	unsigned long count = 0;
//	while (true){
//		ostream << setw(8) << count << "    ";
//		for (int i = 0; i<8; ++i){
//			if (ifs.read((char*)&byte, 1))
//				ostream << setw(2) << (int)byte << " ";
//			else
//				goto endfile;
//		}
//		ostream << " ";
//		for (int i = 0; i<8; ++i){
//			if (ifs.read((char*)&byte, 1))
//				ostream << setw(2) << (int)byte << " ";
//			else
//				goto endfile;
//		}
//		ostream << endl;
//		count += 16;
//	}
//
//endfile:
//	ostream << setfill(' ') << dec;
//}


/**************************************************
*
*function��wchar_t��stringת��
*
**************************************************/
////��Ҫ������ʹ����wchar_t*��delete[]�ͷ��ڴ�
//wchar_t *multiByteToWideChar(string& pKey)
//{
//	//��һ�ε��÷���ת������ַ������ȣ�����ȷ��Ϊwchar_t*���ٶ����ڴ�ռ�
//	int pSize = MultiByteToWideChar(CP_OEMCP, 0, pKey.c_str(), strlen(pKey.c_str()) + 1, NULL, 0);
//	wchar_t *pWCStrKey = new wchar_t[pSize];
//	//�ڶ��ε��ý����ֽ��ַ���ת����˫�ֽ��ַ���
//	MultiByteToWideChar(CP_OEMCP, 0, pKey.c_str(), strlen(pKey.c_str()) + 1, pWCStrKey, pSize);
//	return pWCStrKey;
//}
////��Ҫ����ʹ����char*��delete[]�ͷ��ڴ�
//char* wideCharToMultiByte(wchar_t* pWCStrKey)
//{
//	//��һ�ε���ȷ��ת�����ֽ��ַ����ĳ��ȣ����ڿ��ٿռ�
//	int pSize = WideCharToMultiByte(CP_OEMCP, 0, pWCStrKey, wcslen(pWCStrKey), NULL, 0, NULL, NULL);
//	char* pCStrKey = new char[pSize + 1];
//	//�ڶ��ε��ý�˫�ֽ��ַ���ת���ɵ��ֽ��ַ���
//	WideCharToMultiByte(CP_OEMCP, 0, pWCStrKey, wcslen(pWCStrKey), pCStrKey, pSize, NULL, NULL);
//	pCStrKey[pSize] = '\0';
//	return pCStrKey;
//
//	//�����Ҫת����string��ֱ�Ӹ�ֵ����
//	//string pKey = pCStrKey;
//}