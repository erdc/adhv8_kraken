#!/usr/bin/env python
alphanum='ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789'

lookup = {}
for ii,key in enumerate(alphanum):
    lookup[key]=ii+1

def get_card_id(card):
     assert(len(card)<7)
     base=37
     id=0
     for ii,key in enumerate(card):
         id += lookup[key]*base**(ii)
     #
     return id
