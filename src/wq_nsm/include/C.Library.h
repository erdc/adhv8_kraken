#pragma once
#ifndef C_LIBRARY_H
#define C_LIBRARY_H

#include "LINKAGE.h"



 struct NSMListNode
{
	void* item;

	struct NSMListNode* next;
}; 

typedef struct NSMListNode NSMListNode;



// P R O T O T Y P E S

EXTERN  NSMListNode*  NSMList_PushBackItem(NSMListNode* list, void* ptr);

EXTERN  NSMListNode*  NSMList_DeleteItem(NSMListNode* list, void* ptr);

EXTERN  unsigned long  NSMList_Delete(NSMListNode* list);

EXTERN  unsigned long  NSMList_GetNodeCount(NSMListNode* list);

EXTERN  void  NSMList_EmptyList(NSMListNode* list);

EXTERN  void*  NSMList_GetNextItem(NSMListNode* list, void* lastItem);
#endif
