// PersonalMethodLibrary.cpp : 定义 DLL 应用程序的导出函数。
//

#include "stdafx.h"

/**************************************************
*
* Description:各模块说明
*
**************************************************/
//AES AES加密
//MD5 MD5值计算
//Fit 拟合
//Complex 复数
//Integral 数值积分
//Interpolate 插值
//Matrix 矩阵
//LEquations 求解线性方程组
//NLEquations 求解非线性方程组
//IntegralTransform 积分变换

//DataAnalysis 数值统计分析
//SortAlgorithm 排序算法
//StringBuffer 字符串处理
//Express 公式计算

#ifdef AFX_CLASS  
#define AFX_EX_CLASS _declspec(dllexport)  
#else  
#define AFX_EX_CLASS _declspec(dllimport)  
#endif 


/**************************************************
*
*function：获得DLL名称
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
*function：N皇后问题-求解算法1
*
**************************************************/
//#include<iostream>
//#define N 8
//using namespace std;
//static int gEightQueen[N] = { 0 }, gCount = 0;
//void print()//输出每一种情况下棋盘中皇后的摆放情况
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
//int check_pos_valid(int loop, int value)//检查是否存在有多个皇后在同一行/列/对角线的情况
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
*function：N皇后问题-求解算法2
* 试探-回溯算法，递归实现
*
**************************************************/
//#include "iostream"  
//using namespace std;
//#include "time.h"  
//
//// sum用来记录皇后放置成功的不同布局数；upperlim用来标记所有列都已经放置好了皇后。  
//long sum = 0, upperlim = 1;
//
//// 试探算法从最右边的列开始。  
//void test(long row, long ld, long rd)
//{
//	if (row != upperlim)
//	{
//		// row，ld，rd进行“或”运算，求得所有可以放置皇后的列,对应位为0，  
//		// 然后再取反后“与”上全1的数，来求得当前所有可以放置皇后的位置，对应列改为1  
//		// 也就是求取当前哪些列可以放置皇后  
//		long pos = upperlim & ~(row | ld | rd);
//		while (pos)    // 0 -- 皇后没有地方可放，回溯  
//		{
//			// 拷贝pos最右边为1的bit，其余bit置0  
//			// 也就是取得可以放皇后的最右边的列  
//			long p = pos & -pos;
//
//			// 将pos最右边为1的bit清零  
//			// 也就是为获取下一次的最右可用列使用做准备，  
//			// 程序将来会回溯到这个位置继续试探  
//			pos -= p;
//
//			// row + p，将当前列置1，表示记录这次皇后放置的列。  
//			// (ld + p) << 1，标记当前皇后左边相邻的列不允许下一个皇后放置。  
//			// (ld + p) >> 1，标记当前皇后右边相邻的列不允许下一个皇后放置。  
//			// 此处的移位操作实际上是记录对角线上的限制，只是因为问题都化归  
//			// 到一行网格上来解决，所以表示为列的限制就可以了。显然，随着移位  
//			// 在每次选择列之前进行，原来N×N网格中某个已放置的皇后针对其对角线  
//			// 上产生的限制都被记录下来了  
//			test(row + p, (ld + p) << 1, (rd + p) >> 1);
//		}
//	}
//	else
//	{
//		// row的所有位都为1，即找到了一个成功的布局，回溯  
//		sum++;
//	}
//}
//
//int main(int argc, char *argv[])
//{
//	time_t tm;
//	int n = 16;
//	cout << "输入皇后数：";
//	cin >> n;
//	cout << endl;
//
//	if (argc != 1)
//		n = atoi(argv[1]);
//	tm = time(0);
//
//	// 因为整型数的限制，最大只能32位，  
//	// 如果想处理N大于32的皇后问题，需要  
//	// 用bitset数据结构进行存储  
//	if ((n < 1) || (n > 32))
//	{
//		printf(" 只能计算1-32之间\n");
//		exit(-1);
//	}
//	printf("%d 皇后\n", n);
//
//	// N个皇后只需N位存储，N列中某列有皇后则对应bit置1。  
//	upperlim = (upperlim << n) - 1;
//
//	test(0, 0, 0);
//	printf("共有%ld种排列, 计算时间%d秒 \n", sum, (int)(time(0) - tm));
//	system("pause");
//	return 0;
//}


/**************************************************
*
*function：数独解法
*
**************************************************/
//bool sign = false;
//
// //创建数独矩阵 
//int num[9][9];
//
// //函数声明 
//void Input();
//void Output();
//bool Check(int n, int key);
//void DFS(int n);
//
// //主函数 
//int main()
//{
//	cout << "请输入一个9*9的数独矩阵，空位以0表示:" << endl;
//	Input();
//	DFS(0);
//	Output();
//	system("pause");
//}
//
// //读入数独矩阵 
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
// //输出数独矩阵 
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
// //判断key填入n时是否满足条件 
//bool Check(int n, int key)
//{
//	 //判断n所在横列是否合法 
//	for (int i = 0; i < 9; i++)
//	{
//		 //j为n竖坐标 
//		int j = n / 9;
//		if (num[j][i] == key) return false;
//	}
//
//	 //判断n所在竖列是否合法 
//	for (int i = 0; i < 9; i++)
//	{
//		 //j为n横坐标 
//		int j = n % 9;
//		if (num[i][j] == key) return false;
//	}
//
//	 //x为n所在的小九宫格左顶点竖坐标 
//	int x = n / 9 / 3 * 3;
//
//	 //y为n所在的小九宫格左顶点横坐标 
//	int y = n % 9 / 3 * 3;
//
//	 //判断n所在的小九宫格是否合法 
//	for (int i = x; i < x + 3; i++)
//	{
//		for (int j = y; j < y + 3; j++)
//		{
//			if (num[i][j] == key) return false;
//		}
//	}
//
//	 //全部合法，返回正确 
//	return true;
//}
//
// //深搜构造数独 
//void DFS(int n)
//{
//	 //所有的都符合，退出递归 
//	if (n > 80)
//	{
//		sign = true;
//		return;
//	}
//	 //当前位不为空时跳过 
//	if (num[n / 9][n % 9] != 0)
//	{
//		DFS(n + 1);
//	}
//	else
//	{
//		 //否则对当前位进行枚举测试 
//		for (int i = 1; i <= 9; i++)
//		{
//			 //满足条件时填入数字 
//			if (Check(n, i) == true)
//			{
//				num[n / 9][n % 9] = i;
//				 //继续搜索 
//				DFS(n + 1);
//				 //返回时如果构造成功，则直接退出 
//				if (sign == true) return;
//				 //如果构造不成功，还原当前位 
//				num[n / 9][n % 9] = 0;
//			}
//		}
//	}
//}


/**************************************************
*
*function：完成端口实例
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
////接收数据存放类
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
////完成端口类型判断
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
////资源链表(保证关闭程序时，所有资源正常释放)
//IOCPLIST IOCPList;
////互斥量
//HANDLE hmutex;
//
//void main()
//{
//	//创建互斥量
//	hmutex = CreateMutex(NULL, TRUE, L"ICOPLIST");
//	ReleaseMutex(hmutex);
//
//	//完成端口句柄
//	HANDLE completionPort;
//
//	//初始化套接字
//	WSADATA wsaData;
//	DWORD ret;
//	if (ret = WSAStartup(0x0202, &wsaData) != 0)
//	{
//		std::cout << "WSAStartup failed. Error:" << ret << std::endl;
//		return;
//	}
//	//创建完成端口
//	completionPort = CreateIoCompletionPort(INVALID_HANDLE_VALUE, NULL, 0, 0);
//	if (completionPort == NULL)
//	{
//		std::cout << "CreateIoCompletionPort failed. Error:" << GetLastError() << std::endl;
//		return;
//	}
//	//获得系统信息
//	SYSTEM_INFO mySysInfo;
//	GetSystemInfo(&mySysInfo);
//	// 创建 2 * CPU核数 + 1 个线程  
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
//	//创建监听线程
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
//	//关闭完成端口，即退出所有线程
//	if (completionPort)
//		CloseHandle(completionPort);
//	//确保所有资源释放
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
//	// 启动一个监听socket  
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
//	// 绑定监听端口  
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
//	// 开始死循环，处理数据  
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
//		//客户端IP地址（可能指向中间路由器地址）
//		cout << inet_ntoa(cliaddr.sin_addr) << endl;
//		//客户端IP端口
//		cout << cliaddr.sin_port << endl;
//
//		pHandleData = (LPPER_HANDLE_DATA)GlobalAlloc(GPTR, sizeof(PER_HANDLE_DATA));
//		if (pHandleData == NULL)
//		{
//			std::cout << "GlobalAlloc( HandleData ) failed. Error:" << GetLastError() << std::endl;
//			continue;
//		}
//		//将套接字和完成端口绑定
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
//		// 检查数据是否已经传输完了  
//		if (bytesTransferred == 0)
//		{
//			//退出线程
//			if (!pHandleData && !pIoData)
//				break;
//
//			//客户端断开	//清除列表队列中对象
//			IOCPDATA Info;
//			Info.pHandleData = pHandleData;
//			Info.pIoData = pIoData;
//			WaitForSingleObject(hmutex, INFINITE);
//			IOCPLIST::iterator iter = find(IOCPList.begin(), IOCPList.end(), Info);
//			if (iter != IOCPList.end())
//				IOCPList.erase(iter);
//			ReleaseMutex(hmutex);
//
//			//释放资源
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
//				//数据未接收完、处理已接收数据
//				if (pIoData->pSocketData)
//				{
//					pIoData->pSocketData->add(pIoData->buffer);
//				}
//
//				//清除重叠结构
//				ZeroMemory(&(pIoData->overlapped), sizeof(OVERLAPPED));
//				//继续投放接收I/O请求
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
//				//清除重叠结构			
//				ZeroMemory(&(pIoData->overlapped), sizeof(OVERLAPPED));
//				//输出接收到的数据
//				string str = "Send : \t" + pIoData->pSocketData->m_str + "\r\n";
//				std::cout << str;
//				//投放发送I/O请求
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
//				//数据未发送完，继续投放发送I/O请求
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
//				//清除重叠结构			
//				ZeroMemory(&(pIoData->overlapped), sizeof(OVERLAPPED));
//				//投放接收I/O请求
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
*function：十六进制格式输出
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
*function：wchar_t和string转换
*
**************************************************/
////不要忘记在使用完wchar_t*后delete[]释放内存
//wchar_t *multiByteToWideChar(string& pKey)
//{
//	//第一次调用返回转换后的字符串长度，用于确认为wchar_t*开辟多大的内存空间
//	int pSize = MultiByteToWideChar(CP_OEMCP, 0, pKey.c_str(), strlen(pKey.c_str()) + 1, NULL, 0);
//	wchar_t *pWCStrKey = new wchar_t[pSize];
//	//第二次调用将单字节字符串转换成双字节字符串
//	MultiByteToWideChar(CP_OEMCP, 0, pKey.c_str(), strlen(pKey.c_str()) + 1, pWCStrKey, pSize);
//	return pWCStrKey;
//}
////不要忘记使用完char*后delete[]释放内存
//char* wideCharToMultiByte(wchar_t* pWCStrKey)
//{
//	//第一次调用确认转换后单字节字符串的长度，用于开辟空间
//	int pSize = WideCharToMultiByte(CP_OEMCP, 0, pWCStrKey, wcslen(pWCStrKey), NULL, 0, NULL, NULL);
//	char* pCStrKey = new char[pSize + 1];
//	//第二次调用将双字节字符串转换成单字节字符串
//	WideCharToMultiByte(CP_OEMCP, 0, pWCStrKey, wcslen(pWCStrKey), pCStrKey, pSize, NULL, NULL);
//	pCStrKey[pSize] = '\0';
//	return pCStrKey;
//
//	//如果想要转换成string，直接赋值即可
//	//string pKey = pCStrKey;
//}