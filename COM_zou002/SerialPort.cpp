//////////////////////////////////////////////////////////////////////
// COPYRIGHT NOTICE  
// Copyright (c) 2016, 光巨力信息技术有限公司Algorithm Group  （版权声明）  
// All rights reserved.  
//
// @file    SerialPort.cpp    
// @brief   串口通信函数定义文件  
//
// 打开、配置、读写和关闭等操作
//
//
// @version 1.0     
// @author  LingWei zou    
// @E-mail：zlw2008ok@126.com  
// @date    2016/11/17  
//  
//
// 修订说明： 
////////////////////////////////////////////////////////////////////


#include "StdAfx.h"  
#include "SerialPort.h"  
#include <process.h>  
#include <iostream>  

#define buff_n 1   //读取串口buff大小
/** 线程退出标志 */   
bool SerialPort::s_bExit = false;  

SerialPort::SerialPort(void)
{
	//m_hComm = INVALID_HANDLE_VALUE;  
    m_hListenThread = INVALID_HANDLE_VALUE;  
}
SerialPort::~SerialPort(void)
{
	CloseListenThread();  
    ClosePort(); 
}


bool SerialPort::OpenPort(UINT PortNo )
{
    /** 把串口的编号转换为设备名 */   
    char szPort[50];  
    sprintf_s(szPort, "COM%d", PortNo);  
 
   //  Open a handle to the specified com port.
   hcom = CreateFileA( szPort,
                      GENERIC_READ | GENERIC_WRITE,
                      0,      //  must be opened with exclusive-access
                      NULL,   //  default security attributes
                      OPEN_EXISTING, //  must use OPEN_EXISTING
                      0,//FILE_ATTRIBUTE_NORMAL|FILE_FLAG_OVERLAPPED      //  not overlapped I/O
                      NULL ); //  hTemplate must be NULL for comm devices

   if (hcom == INVALID_HANDLE_VALUE) 
   {
       //  Handle the error.
       //printf ("CreateFile failed with error %d.\n", GetLastError());
       return false;
   } 
   return true;
}

bool SerialPort::InitPort(UINT BaudRate)
{
	DCB Dcb;
	bool fSuccess;
	//DCB Dcbb;
	fSuccess = GetCommState(hcom, &Dcb);	

	Dcb.BaudRate = BaudRate;     //  baud rate
	Dcb.ByteSize = 8;             //  data size, xmit and rcv
	Dcb.Parity   = NOPARITY;      //  parity bit
	Dcb.StopBits = ONESTOPBIT;    //  stop bit

	fSuccess = SetCommState(hcom, &Dcb);
	if (!fSuccess) 
	{
		//  Handle the error.
		printf ("SetCommState failed with error %d.\n", GetLastError());
		return false;
	}
	return true;
}

bool SerialPort::OpenListenThread()	//开启监听线程
{

    /** 检测线程是否已经开启了 */   
    if (m_hListenThread != INVALID_HANDLE_VALUE)  
    {  
        /** 线程已经开启 */   
        return false;  
    }  
 
    s_bExit = false;  
    /** 线程ID */   
    UINT threadId;  
    /** 开启串口数据监听线程 */   
    m_hListenThread = (HANDLE)_beginthreadex(NULL, 0, ListenThread, this, 0, &threadId);  
    if (!m_hListenThread)  
    {  
        return false;  
    }  

    /** 设置线程的优先级,高于普通线程 */   
    if (!SetThreadPriority(m_hListenThread, THREAD_PRIORITY_ABOVE_NORMAL))  
    {  
        return false;  
    }  

	return true;
}
UINT SerialPort::ListenThread(void* pParam)
{
	 /** 得到本类的指针 */   
	SerialPort *pSerialPort = reinterpret_cast<SerialPort*>(pParam);
    char str;
	while(!pSerialPort->s_bExit)
	{
		pSerialPort->ReadPort(&str);
 
		printf("%c",str);
	}

	return 0;
}
bool SerialPort::CloseListenThread()	//关闭监听线程
{
	if (m_hListenThread != INVALID_HANDLE_VALUE)  
    {  
        /** 通知线程退出 */   
        s_bExit = true;  
 
        /** 等待线程退出 */   
        Sleep(10);  

        /** 置线程句柄无效 */   
        CloseHandle( m_hListenThread );  
        m_hListenThread = INVALID_HANDLE_VALUE;  
    }
	return true;
}
bool SerialPort::ReadPort(char *str)//
{
	//char str[100];
	DWORD wCount;
	bool bReadStat;
	bReadStat=ReadFile(hcom,str,buff_n,&wCount,NULL);//第三个参数 为每次读取的字节数，读取满了就返回函数，没有读取满一直等待。。。，读取满后将数据存入str中
	
	if(!bReadStat)
	{
		printf("读取串口失败！");
		return false;
	}

	return true;
}

bool SerialPort::WritePort(unsigned char* WData, unsigned int length)
{
	return true;
}
bool SerialPort::ClearPort()
{
	return PurgeComm(hcom, PURGE_RXCLEAR | PURGE_RXABORT);  //清除串口缓冲区
}
bool SerialPort::ClosePort()
{
	return true;
}
