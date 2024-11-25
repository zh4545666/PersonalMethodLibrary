﻿#include "stdafx.h"
#include "AES.h"

using namespace PersonalMethod;

//permutebox  
static unsigned char permutebox[] =  
{ /*  0    1    2    3    4    5    6    7    8    9    a    b    c    d    e    f */  
    0x63,0x7c,0x77,0x7b,0xf2,0x6b,0x6f,0xc5,0x30,0x01,0x67,0x2b,0xfe,0xd7,0xab,0x76, /*0*/  
    0xca,0x82,0xc9,0x7d,0xfa,0x59,0x47,0xf0,0xad,0xd4,0xa2,0xaf,0x9c,0xa4,0x72,0xc0, /*1*/  
    0xb7,0xfd,0x93,0x26,0x36,0x3f,0xf7,0xcc,0x34,0xa5,0xe5,0xf1,0x71,0xd8,0x31,0x15, /*2*/  
    0x04,0xc7,0x23,0xc3,0x18,0x96,0x05,0x9a,0x07,0x12,0x80,0xe2,0xeb,0x27,0xb2,0x75, /*3*/  
    0x09,0x83,0x2c,0x1a,0x1b,0x6e,0x5a,0xa0,0x52,0x3b,0xd6,0xb3,0x29,0xe3,0x2f,0x84, /*4*/  
    0x53,0xd1,0x00,0xed,0x20,0xfc,0xb1,0x5b,0x6a,0xcb,0xbe,0x39,0x4a,0x4c,0x58,0xcf, /*5*/  
    0xd0,0xef,0xaa,0xfb,0x43,0x4d,0x33,0x85,0x45,0xf9,0x02,0x7f,0x50,0x3c,0x9f,0xa8, /*6*/  
    0x51,0xa3,0x40,0x8f,0x92,0x9d,0x38,0xf5,0xbc,0xb6,0xda,0x21,0x10,0xff,0xf3,0xd2, /*7*/  
    0xcd,0x0c,0x13,0xec,0x5f,0x97,0x44,0x17,0xc4,0xa7,0x7e,0x3d,0x64,0x5d,0x19,0x73, /*8*/  
    0x60,0x81,0x4f,0xdc,0x22,0x2a,0x90,0x88,0x46,0xee,0xb8,0x14,0xde,0x5e,0x0b,0xdb, /*9*/  
    0xe0,0x32,0x3a,0x0a,0x49,0x06,0x24,0x5c,0xc2,0xd3,0xac,0x62,0x91,0x95,0xe4,0x79, /*a*/  
    0xe7,0xc8,0x37,0x6d,0x8d,0xd5,0x4e,0xa9,0x6c,0x56,0xf4,0xea,0x65,0x7a,0xae,0x08, /*b*/  
    0xba,0x78,0x25,0x2e,0x1c,0xa6,0xb4,0xc6,0xe8,0xdd,0x74,0x1f,0x4b,0xbd,0x8b,0x8a, /*c*/  
    0x70,0x3e,0xb5,0x66,0x48,0x03,0xf6,0x0e,0x61,0x35,0x57,0xb9,0x86,0xc1,0x1d,0x9e, /*d*/  
    0xe1,0xf8,0x98,0x11,0x69,0xd9,0x8e,0x94,0x9b,0x1e,0x87,0xe9,0xce,0x55,0x28,0xdf, /*e*/  
    0x8c,0xa1,0x89,0x0d,0xbf,0xe6,0x42,0x68,0x41,0x99,0x2d,0x0f,0xb0,0x54,0xbb,0x16  /*f*/  
};  
//inversepermutationbox  
static  unsigned char inversepermutationbox[]=   
{ /*  0    1    2    3    4    5    6    7    8    9    a    b    c    d    e    f  */  
    0x52,0x09,0x6a,0xd5,0x30,0x36,0xa5,0x38,0xbf,0x40,0xa3,0x9e,0x81,0xf3,0xd7,0xfb, /*0*/  
    0x7c,0xe3,0x39,0x82,0x9b,0x2f,0xff,0x87,0x34,0x8e,0x43,0x44,0xc4,0xde,0xe9,0xcb, /*1*/  
    0x54,0x7b,0x94,0x32,0xa6,0xc2,0x23,0x3d,0xee,0x4c,0x95,0x0b,0x42,0xfa,0xc3,0x4e, /*2*/  
    0x08,0x2e,0xa1,0x66,0x28,0xd9,0x24,0xb2,0x76,0x5b,0xa2,0x49,0x6d,0x8b,0xd1,0x25, /*3*/  
    0x72,0xf8,0xf6,0x64,0x86,0x68,0x98,0x16,0xd4,0xa4,0x5c,0xcc,0x5d,0x65,0xb6,0x92, /*4*/  
    0x6c,0x70,0x48,0x50,0xfd,0xed,0xb9,0xda,0x5e,0x15,0x46,0x57,0xa7,0x8d,0x9d,0x84, /*5*/  
    0x90,0xd8,0xab,0x00,0x8c,0xbc,0xd3,0x0a,0xf7,0xe4,0x58,0x05,0xb8,0xb3,0x45,0x06, /*6*/  
    0xd0,0x2c,0x1e,0x8f,0xca,0x3f,0x0f,0x02,0xc1,0xaf,0xbd,0x03,0x01,0x13,0x8a,0x6b, /*7*/  
    0x3a,0x91,0x11,0x41,0x4f,0x67,0xdc,0xea,0x97,0xf2,0xcf,0xce,0xf0,0xb4,0xe6,0x73, /*8*/  
    0x96,0xac,0x74,0x22,0xe7,0xad,0x35,0x85,0xe2,0xf9,0x37,0xe8,0x1c,0x75,0xdf,0x6e, /*9*/  
    0x47,0xf1,0x1a,0x71,0x1d,0x29,0xc5,0x89,0x6f,0xb7,0x62,0x0e,0xaa,0x18,0xbe,0x1b, /*a*/  
    0xfc,0x56,0x3e,0x4b,0xc6,0xd2,0x79,0x20,0x9a,0xdb,0xc0,0xfe,0x78,0xcd,0x5a,0xf4, /*b*/  
    0x1f,0xdd,0xa8,0x33,0x88,0x07,0xc7,0x31,0xb1,0x12,0x10,0x59,0x27,0x80,0xec,0x5f, /*c*/  
    0x60,0x51,0x7f,0xa9,0x19,0xb5,0x4a,0x0d,0x2d,0xe5,0x7a,0x9f,0x93,0xc9,0x9c,0xef, /*d*/  
    0xa0,0xe0,0x3b,0x4d,0xae,0x2a,0xf5,0xb0,0xc8,0xeb,0xbb,0x3c,0x83,0x53,0x99,0x61, /*e*/  
    0x17,0x2b,0x04,0x7e,0xba,0x77,0xd6,0x26,0xe1,0x69,0x14,0x63,0x55,0x21,0x0c,0x7d  /*f*/  
};   
  
  
  
//CAES::CAES(void)  
//{  
//}  
CAES::CAES(unsigned char* key , size_t N )   
{ 
	if (N < 16)
	{
		unsigned char keyData[16] = {0};
		memcpy(keyData,key,sizeof(unsigned char)*N);
		KeyExpansion(keyData, swapbox);
	}
	else
		KeyExpansion(key, swapbox);  
}  
  
CAES::~CAES(void)  
{  
}  

/************************************************************************/
/* create a Encrypt method                                              */
/* param : pSrcPath     file path                                       */
/* param : pDestPath    file path                                       */
/* return : bool                                                        */
/************************************************************************/
bool CAES::Encrypt(char* pSrcFilePath, char* pDestFilePath)
{
	ifstream fs;
	fs.open(pSrcFilePath, ios::in | ios::binary);
	if (!fs.is_open())
		return false;
	ofstream fp;
	fp.open(pDestFilePath, ios::out | ios::trunc | ios::binary);
	if (!fp.is_open())
	{
		fs.close();
		return false;
	}

	unsigned char input[16] = { 0 };
	while (!fs.eof())
	{
		fs.read((char *)input, sizeof(input));
		size_t DataCount = (size_t)fs.gcount();
		Encrypt(input, DataCount);
		fp.write((char *)input, DataCount);
		memset(input, 0, sizeof(input));
	}

	fs.close();
	fp.close();
	return true;
}

/************************************************************************/
/* create a Encrypt method                                              */
/* param : pSrcPath     file path                                       */
/* param : pDestPath    file path                                       */
/* return : bool                                                        */
/************************************************************************/
bool CAES::Decrypt(char* pSrcFilePath, char* pDestFilePath)
{
	ifstream fs;
	fs.open(pSrcFilePath, ios::in | ios::binary);
	if (!fs.is_open())
		return false;
	ofstream fp;
	fp.open(pDestFilePath, ios::out | ios::trunc | ios::binary);
	if (!fp.is_open())
	{
		fs.close();
		return false;
	}

	unsigned char input[16] = { 0 };
	while (!fs.eof())
	{
		fs.read((char *)input, sizeof(input));
		size_t DataCount = (size_t)fs.gcount();
		Decrypt(input, DataCount);
		fp.write((char *)input, DataCount);
		memset(input, 0, sizeof(input));
	}

	fs.close();
	fp.close();
	return true;
}
  
/************************************************************************/  
/* create a Encrypt method                                              */  
/* param : data     encrypt data                                        */  
/* param :encryptArrayencryptArray                                   */  
/* param :len encrypt data length                                      */  
/* return : void                                                        */  
/************************************************************************/  
void CAES::Encrypt(unsigned char* data ,unsigned char * encryptArray,size_t len)  
{  
    memcpy(encryptArray,data,len);  
    Cipher((void *)encryptArray,len);  
  
}
void CAES::Encrypt(unsigned char * encryptArray, size_t len)
{
	Cipher((void *)encryptArray, len);
}
  
/************************************************************************/  
/* create a Decrypt method                                              */  
/* param : data     decrypt data                                        */  
/* param :decryptArraydecryptArray                                   */  
/* param :len decrypt data length                                      */  
/* return : void                                                        */  
/************************************************************************/  
void CAES::Decrypt(unsigned char * data,unsigned char * decryptArray,size_t len)  
{     
    memcpy(decryptArray,data,len);  
    InvCipher((void *)decryptArray,len);  
}
void CAES::Decrypt(unsigned char * decryptArray, size_t len)
{
	InvCipher((void *)decryptArray, len);
}
  
/************************************************************************/  
/* create a Cipher  method  only one time  Encrypt                      */  
/* param : input     input encrypt data                                 */  
/* return : unsigned char *                                             */  
/************************************************************************/    
unsigned char* CAES::Cipher(unsigned char* input)  
{  
  
    unsigned char state[AES_KEY_ROW_NUMBER][AES_KEY_COLUMN_NUMBER];  
    int i,r,c;  
  
    for(r=0; r<AES_KEY_ROW_NUMBER; r++)  
    {  
        for(c=0; c<AES_KEY_COLUMN_NUMBER ;c++)  
        {  
            state[r][c] = input[c*AES_KEY_COLUMN_NUMBER+r];  
        }  
    }  
  
    AddRoundKey(state,swapbox[0]);  
  
    for(i=1; i<=AES_ROUND_COUNT; i++)  
    {  
        SubBytes(state);  
        ShiftRows(state);  
        if(i!=AES_ROUND_COUNT)MixColumns(state);  
        AddRoundKey(state,swapbox[i]);  
    }  
  
    for(r=0; r<AES_KEY_ROW_NUMBER; r++)  
    {  
        for(c=0; c<AES_KEY_COLUMN_NUMBER ;c++)  
        {  
            input[c*AES_KEY_COLUMN_NUMBER+r] = state[r][c];  
        }  
    }  
    return input;  
}  

/************************************************************************/  
/* create a InvCipher  method  only one time  decrypt                   */  
/* param : input     input decrypt data                                 */  
/* return : unsigned char *                                             */  
/************************************************************************/  
unsigned char* CAES::InvCipher(unsigned char* input)  
{  
    unsigned char state[AES_KEY_ROW_NUMBER][AES_KEY_COLUMN_NUMBER];  
    int i,r,c;  
  
    for(r=0; r<AES_KEY_ROW_NUMBER; r++)  
    {  
        for(c=0; c<AES_KEY_COLUMN_NUMBER ;c++)  
        {  
            state[r][c] = input[c*AES_KEY_COLUMN_NUMBER+r];  
        }  
    }  
  
    AddRoundKey(state, swapbox[10]);  
    for(i=9; i>=0; i--)  
    {  
        InvShiftRows(state);  
        InvSubBytes(state);  
        AddRoundKey(state, swapbox[i]);  
        if(i)  
        {  
            InvMixColumns(state);  
        }  
    }     
    for(r=0; r<AES_KEY_ROW_NUMBER; r++)  
    {  
        for(c=0; c<AES_KEY_COLUMN_NUMBER ;c++)  
        {  
            input[c*AES_KEY_COLUMN_NUMBER+r] = state[r][c];  
        }  
    }  
    return input;  
}  

/************************************************************************/  
/* Create a specified length of data encryption method                  */  
/* param : input     input data encryption                              */  
/* param : length    Input the length of the encrypted data             */  
/* return : unsigned char *                                             */  
/************************************************************************/  
unsigned char* CAES::Cipher(void * input, size_t length)  
{  
    unsigned char* in = (unsigned char*) input;  
    size_t i,j;  
    if(!length)  
    {  
        while(*(in+length++));  
        in = (unsigned char*) input;  
    }  
    for(i=0; i<length; i+=16)  
    {  
		if ((length - i) < 16)
		{
			for (j = 0; j<(length - i); j++)
			{
				in[i + j] = 255 - in[i + j];
			}
		}
		else
			Cipher(in + i);
    }  
    return (unsigned char*)input;  
}  

/************************************************************************/  
/* Create a specified length of InvCipher method                        */  
/* param : input     input data InvCipher                               */  
/* param : length    Input the length of the InvCipher data             */  
/* return : unsigned char *                                             */  
/************************************************************************/  
unsigned char* CAES::InvCipher(void * input, size_t length)  
{  
    unsigned char* in = (unsigned char*) input;  
    size_t i,j;  
    for(i=0; i<length; i+=16)  
    {  
		if ((length - i) < 16)
		{
			for (j = 0; j<(length - i); j++)
			{
				in[i + j] = 255 - in[i + j];
			}
		}
		else
			InvCipher(in + i);
    }  
    return (unsigned char*)input;  
}  

/************************************************************************/  
/*Create key method                                                     */  
/* param : key      input data encryption key                           */  
/* param :swapbox  Conversion of key array                             */  
/* return : void                                                        */  
/************************************************************************/  
void CAES::KeyExpansion(unsigned char* key, unsigned char swapbox[][4][4])  
{  
    int i,j,r,c;  
    unsigned char rc[] = {0x01, 0x02, 0x04, 0x08, 0x10, 0x20, 0x40, 0x80, 0x1b, 0x36};  
    for(r=0; r<AES_KEY_ROW_NUMBER; r++)  
    {  
        for(c=0; c<AES_KEY_COLUMN_NUMBER; c++)  
        {  
            swapbox[0][r][c] = key[r+c*AES_KEY_COLUMN_NUMBER];  
        }  
    }  
    for(i=1; i<=10; i++)  
    {  
        for(j=0; j<AES_KEY_COLUMN_NUMBER; j++)  
        {  
            unsigned char t[AES_KEY_ROW_NUMBER];  
            for(r=0; r<AES_KEY_ROW_NUMBER; r++)  
            {  
                t[r] = j ? swapbox[i][r][j-1] : swapbox[i-1][r][3];  
            }  
            if(j == 0)  
            {  
                unsigned char temp = t[0];  
                for(r=0; r<AES_KEY_ROW_NUMBER-1; r++)  
                {  
                    t[r] = permutebox[t[(r+1)%AES_KEY_ROW_NUMBER]];  
                }  
                t[3] = permutebox[temp];  
                t[0] ^= rc[i-1];  
            }  
            for(r=0; r<AES_KEY_ROW_NUMBER; r++)  
            {  
                swapbox[i][r][j] = swapbox[i-1][r][j] ^ t[r];  
            }  
        }  
    }  
}  
  
/************************************************************************/  
/*Create mixed operation method ranks                                   */  
/* param : a      row  char                                             */  
/* param : b      column  char                                          */  
/* return : unsigned char                                               */  
/************************************************************************/    
unsigned char CAES::FFmul(unsigned char a, unsigned char b)  
{  
    unsigned char bw[AES_KEY_ROW_NUMBER];  
    unsigned char res=0;  
    int i;  
    bw[0] = b;  
    for(i=1; i<AES_KEY_ROW_NUMBER; i++)  
    {  
        bw[i] = bw[i-1]<<1;  
        if(bw[i-1]&0x80)  
        {  
            bw[i]^=0x1b;  
        }  
    }  
    for(i=0; i<AES_KEY_ROW_NUMBER; i++)  
    {  
        if((a>>i)&0x01)  
        {  
            res ^= bw[i];  
        }  
    }  
    return res;  
}  

/************************************************************************/  
/* Create bytes alternatives                                            */  
/* param :  state[][]  Byte array alternative                           */  
/* return : void                                                        */  
/************************************************************************/  
void CAES::SubBytes(unsigned char state[][AES_KEY_COLUMN_NUMBER])  
{  
    int r,c;  
    for(r=0; r<AES_KEY_ROW_NUMBER; r++)  
    {  
        for(c=0; c<AES_KEY_COLUMN_NUMBER; c++)  
        {  
            state[r][c] = permutebox[state[r][c]];  
        }  
    }  
}  

/************************************************************************/  
/* Create rows transform method                                         */  
/* param :  state[][]  line array alternative                           */  
/* return : void                                                        */  
/************************************************************************/  
void CAES::ShiftRows(unsigned char state[][AES_KEY_COLUMN_NUMBER])  
{  
    unsigned char t[AES_KEY_COLUMN_NUMBER];  
    int r,c;  
    for(r=1; r<AES_KEY_ROW_NUMBER; r++)  
    {  
        for(c=0; c<AES_KEY_COLUMN_NUMBER; c++)  
        {  
            t[c] = state[r][(c+r)%AES_KEY_COLUMN_NUMBER];  
        }  
        for(c=0; c<AES_KEY_COLUMN_NUMBER; c++)  
        {  
            state[r][c] = t[c];  
        }  
    }  
}  

/************************************************************************/  
/* Create columns transform method                                      */  
/* param :  state[][]  columns array alternative                        */  
/* return : void                                                        */  
/************************************************************************/  
void CAES::MixColumns(unsigned char state[][AES_KEY_COLUMN_NUMBER])  
{  
    unsigned char t[AES_KEY_ROW_NUMBER];  
    int r,c;  
    for(c=0; c< AES_KEY_COLUMN_NUMBER; c++)  
    {  
        for(r=0; r<AES_KEY_ROW_NUMBER; r++)  
        {  
            t[r] = state[r][c];  
        }  
        for(r=0; r<AES_KEY_ROW_NUMBER; r++)  
        {  
            state[r][c] = FFmul(0x02, t[r])  
                ^ FFmul(0x03, t[(r+1)%AES_KEY_COLUMN_NUMBER])  
                ^ FFmul(0x01, t[(r+2)%AES_KEY_COLUMN_NUMBER])  
                ^ FFmul(0x01, t[(r+3)%AES_KEY_COLUMN_NUMBER]);  
        }  
    }  
}  
  
/************************************************************************/  
/*Create round keys plus transform method                               */  
/* param :  state[][]  keys plus array alternative                      */  
/* param :  k[][]  temp array alternative                               */  
/* return : void                                                        */  
/************************************************************************/  
void CAES::AddRoundKey(unsigned char state[][AES_KEY_COLUMN_NUMBER], unsigned char k[][AES_KEY_COLUMN_NUMBER])  
{  
    int r,c;  
    for(c=0; c<AES_KEY_COLUMN_NUMBER; c++)  
    {  
        for(r=0; r<AES_KEY_ROW_NUMBER; r++)  
        {  
            state[r][c] ^= k[r][c];  
        }  
    }  
}  

/************************************************************************/  
/* CreateInvSubBytes alternatives                                      */  
/* param :  state[][]  InvSubBytes array alternative                    */  
/* return : void                                                        */  
/************************************************************************/  
void CAES::InvSubBytes(unsigned char state[][AES_KEY_COLUMN_NUMBER])  
{  
    int r,c;  
    for(r=0; r<AES_KEY_ROW_NUMBER; r++)  
    {  
        for(c=0; c<AES_KEY_COLUMN_NUMBER; c++)  
        {  
            state[r][c] = inversepermutationbox[state[r][c]];  
        }  
    }  
}  

/************************************************************************/  
/* CreateInvShiftRows transform method                                 */  
/* param :  state[][]  InvShiftRows array alternative                   */  
/* return : void                                                        */  
/************************************************************************/  
void CAES::InvShiftRows(unsigned char state[][AES_KEY_COLUMN_NUMBER])  
{  
    unsigned char t[AES_KEY_COLUMN_NUMBER];  
    int r,c;  
    for(r=1; r<AES_KEY_ROW_NUMBER; r++)  
    {  
        for(c=0; c<AES_KEY_COLUMN_NUMBER; c++)  
        {  
            t[c] = state[r][(c-r+AES_KEY_COLUMN_NUMBER)%AES_KEY_COLUMN_NUMBER];  
        }  
        for(c=0; c<AES_KEY_COLUMN_NUMBER; c++)  
        {  
            state[r][c] = t[c];  
        }  
    }  
}  

/************************************************************************/  
/* CreateInvMixColumns transform method                                */  
/* param :  state[][]  InvMixColumns array alternative                  */  
/* return : void                                                        */  
/************************************************************************/  
void CAES::InvMixColumns(unsigned char state[][AES_KEY_COLUMN_NUMBER])  
{  
    unsigned char t[AES_KEY_ROW_NUMBER];  
    int r,c;  
    for(c=0; c< AES_KEY_COLUMN_NUMBER; c++)  
    {  
        for(r=0; r<AES_KEY_ROW_NUMBER; r++)  
        {  
            t[r] = state[r][c];  
        }  
        for(r=0; r<AES_KEY_ROW_NUMBER; r++)  
        {  
            state[r][c] = FFmul(0x0e, t[r])  
                ^ FFmul(0x0b, t[(r+1)%AES_KEY_ROW_NUMBER])  
                ^ FFmul(0x0d, t[(r+2)%AES_KEY_ROW_NUMBER])  
                ^ FFmul(0x09, t[(r+3)%AES_KEY_ROW_NUMBER]);  
        }  
    }  
}  