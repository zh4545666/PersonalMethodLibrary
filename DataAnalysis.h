#pragma once

#include <map>
#include <queue>

#ifdef AFX_CLASS  
#define AFX_EX_CLASS _declspec(dllexport)  
#else  
#define AFX_EX_CLASS _declspec(dllimport)  
#endif 

class AFX_EX_CLASS CDataAnalysis
{
public:
	CDataAnalysis();
	~CDataAnalysis();

public:
	static void GetSTD(double* pData, size_t N, double pResult[3]);
	
	
	//
	// 循环冗余校验CRC32
	//

public:
	unsigned long crc(unsigned char *buf, int len);
private:
	unsigned long crc_table[256];
	int crc_table_computed = 0;
	void make_crc_table(void);
	unsigned long update_crc(unsigned long crc, unsigned char *buf, int len);

	//
	// Huffman树和Huffman编码
	//

public:
	enum { LEN = 512 };
	struct huffman_node
	{
		unsigned char c;
		int weight;
		char hufman_code[LEN];
		bool bValue;
		huffman_node* left;
		huffman_node* right;

		huffman_node()
		{
			c = '\0';
			weight = -1;
			left = NULL;
			right = NULL;
			bValue = false;
			memset(hufman_code, 0, sizeof(char)*LEN);
		}
	};
	huffman_node* huffman_tree_create(map<unsigned char, int> &word, map<unsigned char, string> &res);
	static void huffman_tree_del(huffman_node* pnode);
	static bool huffman_node_compare(const huffman_node* other, const huffman_node* other2);
private:
	int get_huffman_code(huffman_node* pnode);
	void print_huffman_pre(huffman_node* pnode, map<unsigned char, string>& res);

	//
	// 十六进制输出和转换
	//

public:
	void printHex(std::ifstream& ifs, std::ostream& ostream);
	static int UnsignedCharToInt(unsigned char c1, unsigned char c2);
	static void IntToUnsignedChar(int val, unsigned char res[2]);
	static string UnsignedCharToString(unsigned char c);


	//
	// KMP字符串匹配
	//

public:
	static int KmpMatch(char* s, char* p);
	static void GetNextval(char* p, int next[]);
};


