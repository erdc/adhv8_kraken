#pragma once

#ifdef _USRDLL
#define EXPORT __declspec(dllexport)
#else
#define EXPORT
#endif



#ifdef __cplusplus
#define EXTERN extern "C"
#else
#define EXTERN
#endif

//format:  EXTERN  EXPORT void print(...)







