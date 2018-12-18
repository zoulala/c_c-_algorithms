//////////////////////////////////////////////////////////////////////
// COPYRIGHT NOTICE  
// Copyright (c) 2016, 光巨力信息技术有限公司Algorithm Group  （版权声明）  
// All rights reserved.  
//
// @file    SerialPort.h    
// @brief   串口通信头文件  
//
// 创建串口类
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

#ifndef SerialPort_h
#define SerialPort_h


#include <Windows.h>  

class SerialPort
{
public:
    SerialPort(void);  
    ~SerialPort(void); 

public:
	bool OpenPort(UINT PortNo );
	bool InitPort(UINT BaudRate);
	bool OpenListenThread();	//开启监听线程
    static UINT WINAPI ListenThread(void* pParam);//线程执行操作
	bool CloseListenThread();	//关闭监听线程
	bool ReadPort(char *str);
	bool WritePort(unsigned char* WData, unsigned int length);
	bool ClearPort();
	bool ClosePort();

private:
	HANDLE hcom;

	/** 线程退出标志变量 */   
   static bool s_bExit;  
 
    /** 线程句柄 */   
    volatile HANDLE    m_hListenThread;  

};



#endif