#!/bin/sh
##Checks if current directory matches expected current directory
dir_check () {
    correct_cur_Dir=$1
    cur_Dir=$(basename $(pwd))
    if [ $cur_Dir != $correct_cur_Dir ]
    then
        printf '%s\n' "In wrong directory should be in " >&1
        echo $correct_cur_Dir >&1
        printf '%s\n' "But instead is in  " >&1
        echo $cur_Dir >&1
        printf '%s\n' "In wrong directory should be in " >&2
        echo $correct_cur_Dir >&2
        printf '%s\n' "But instead is in  " >&2
        echo $cur_Dir >&2
        exit 1
    fi
    return 0
}

##Check if file exists and has size >0 in current directory
file_exist_check () {
    file_name_to_check=$1
    cur_Dir=$(basename $(pwd))
    flag=0
    if [ -s $file_name_to_check ]
    then
        flag=0
    else
        flag=1
        printf '%s\n' "File  " >&1
        echo $file_name_to_check >&1
        printf '%s\n' "was not found in directory  " >&1
        echo $cur_Dir >&1
        printf '%s\n' "File  " >&2
        echo $file_name_to_check >&2
        printf '%s\n' "was not found in directory  " >&2
        echo $cur_Dir >&2
        exit 1
    fi
    return 0
}

##Return truncated job id
get_trunc_hold_jid() {
    hold_jid=$1
    hold_jid_trunc=$2
    OLD_IFS=$IFS
    IFS="."
    newArray_1=($hold_jid)
    hold_jid_trunc=${newArray_1[0]}
    IFS=$OLD_IFS
}