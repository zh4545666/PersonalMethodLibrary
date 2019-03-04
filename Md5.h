#ifndef CMD5_H  
#define CMD5_H  

/* Type define */
typedef unsigned char byte;
typedef unsigned int uint32;

#ifdef AFX_CLASS  
#define AFX_EX_CLASS _declspec(dllexport)  
#else  
#define AFX_EX_CLASS _declspec(dllimport)  
#endif

/* MD5 declaration. */
class AFX_EX_CLASS CMD5 {
public:
	//
	// 4个构造函数
	//

	CMD5();                                           //默认构造函数
	CMD5(const void *input, size_t length);           //输入内存地址与长度信息的构造函数
	CMD5(const string &str);                          //输入字符串的构造函数
	CMD5(ifstream &in);                               //输入流的构造函数

	//
	// 3个Update函数
	//

	void update(const void *input, size_t length);    //往CMD5对象内添加内存块
	void update(const string &str);                   //添加字符串
	void update(ifstream &in);                        //添加流

	//
	// 计算输出
	//

	const byte* digest();                             //计算MD5码,并返回指向它的指针
	string toString();                                //计算MD5码,并返回其对应的字符串
	void reset();                                     //重置
private:
	void update(const byte *input, size_t length);
	void final();
	void transform(const byte block[64]);
	void encode(const uint32 *input, byte *output, size_t length);
	void decode(const byte *input, uint32 *output, size_t length);
	string bytesToHexString(const byte *input, size_t length);

	/* class uncopyable */
	CMD5(const CMD5&);
	CMD5& operator=(const CMD5&);
private:
	uint32 _state[4];   /* state (ABCD) */
	uint32 _count[2];   /* number of bits, modulo 2^64 (low-order word first) */
	byte _buffer[64];   /* input buffer */
	byte _digest[16];   /* message digest */
	bool _finished;     /* calculate finished ? */

	static const byte PADDING[64];  /* padding for calculate */
	static const char HEX[16];
	static const size_t BUFFER_SIZE = 1024;
};

#endif/*CMD5_H*/ 
