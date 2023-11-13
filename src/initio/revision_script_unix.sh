#!/bin/sh
REV=`svnversion ./`
printf "#define ADHREV \"%s\"\n" $REV > $1
DATE=`date "+%Y.%m.%d"`
printf "#define ADHREVDATE \"%s\"\n" $DATE >> $1
TIME=`date "+%H:%M:%S"`
printf "#define ADHREVTIME \"%s\"\n" $TIME >> $1
