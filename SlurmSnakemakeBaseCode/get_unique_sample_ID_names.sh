#!/bin/sh
#COUNTER=$SGE_TASK_ID
projectID=$1
left_read_file=$2
right_read_file=$3
## Function that reads in the  associated projected ID and left read file name
## Outputs a unique identifier for the left read (assuming that the left read file name is unique)
get_left_read_unique_ID () {
    projectID=$1
    left_read_file=$2
    left_read_file_base_name=$3
    echo "Left read file"
    echo $left_read_file
    rel_path_to_proj_dir="../"
    left_read_file_rel_path_from_script=${rel_path_to_proj_dir}$left_read_file
    echo 'Left read file rel path from script'
    echo $left_read_file_rel_path_from_script
    ## The command below does several things:
    ### i) Split based on ., and extract the portion to the left of the ., which will be the relative path for the read file without the 
    ### file type extension relative to the project directory. The path may still contain folders though, and for the unique identifier we 
    ### want to remove folder names and slashes. We reverse the order of the path string, then split based on '/'(corresponding to 
    ### directories in the path). We then take the first column of the reversed delimited string. This corresponds to the file name of the 
    ### input read file itself without the file type extension or any directory information, except that the file name is in reverse order.
    ### We then reverse it again, to return it to the correct order. 
    left_read_file_base_name=$(echo ${left_read_file} | cut -d'.' -f1| rev | cut -d'/' -f 1 | rev ) 
    echo "Left read file base name"
    echo $left_read_file_base_name
    echo "Done with left read base name print out"
    ## Check that input read files are either in .fastq or .fq.gz format
    left_file_format=$(echo ${left_read_file} | cut -d'.' -f2)
    if [ $left_file_format != "fastq" ]
    then
        if [ $left_file_format != "fq" ]
        then
            printf '%s\n' "left file format invalid " >&1
            printf '%s\n' "PathSeq only accepts .fq.gz and .fastq formats " >&1
            echo $left_file_format >&1
            echo $left_read_file_base_name >&1
            echo $left_read_file >&1
            printf '%s\n' "left file format invalid " >&2
            printf '%s\n' "PathSeq only accepts .fq.gz and .fastq formats " >&2
            echo $left_file_format >&2
            echo $left_read_file_base_name >&2
            echo $left_read_file >&2
            exit 1
        fi
        gz_left_format=$(echo ${left_read_file} | cut -d'.' -f3)
        if [ $gz_left_format != "gz" ]
        then
            printf '%s\n' "left file format invalid: has .fq but no .gz " >&1
            printf '%s\n' "PathSeq only accepts .fq.gz and .fastq formats " >&1
            echo $gz_left_format >&1
            echo $left_read_file_base_name >&1
            printf '%s\n' "left file format invalid: has .fq but no .gz " >&2
            printf '%s\n' "PathSeq only accepts .fq.gz and .fastq formats " >&2
            echo $gz_left_format >&2
            echo $left_read_file_base_name >&2
            exit 1
        fi
        echo "left file format is .fq.gz"
        echo $gz_left_format
    fi
    echo $left_read_file_base_name
    echo "Created sample unique ID for left read"
    #return $left_read_file_base_name
}


get_right_read_unique_ID () {
    projectID=$1
    right_read_file=$2
    right_read_file_base_name=$3


    echo "Line 28"
    echo $right_read_file

    rel_path_to_proj_dir="../"

    right_read_file_rel_path_from_script=${rel_path_to_proj_dir}$right_read_file

    ## The command below does several things:
    ### i) Split based on ., and extract the portion to the left of the ., which will be the relative path for the read file without the 
    ### file type extension relative to the project directory. The path may still contain folders though, and for the unique identifier we 
    ### want to remove folder names and slashes. We reverse the order of the path string, then split based on '/'(corresponding to 
    ### directories in the path). We then take the first column of the reversed delimited string. This corresponds to the file name of the 
    ### input read file itself without the file type extension or any directory information, except that the file name is in reverse order.
    ### We then reverse it again, to return it to the correct order. 
    right_read_file_base_name=$(echo ${right_read_file} | cut -d'.' -f1| rev | cut -d'/' -f 1 | rev ) 
    
    echo "Right read file base name"
    echo $right_read_file_base_name
    echo "Done with right read base name print out"

    ## Check that input read files are either in .fastq or .fq.gz format
    right_file_format=$(echo ${right_read_file} | cut -d'.' -f2)
    if [ $right_file_format != "fastq" ]
    then
        if [ $right_file_format != "fq" ]
        then
            printf '%s\n' "right file format invalid " >&1
            printf '%s\n' "PathSeq only accepts .fq.gz and .fastq formats " >&1
            echo $right_file_format >&1
            echo $right_read_file_base_name >&1
            echo $right_read_file >&1

            printf '%s\n' "right file format invalid " >&2
            printf '%s\n' "PathSeq only accepts .fq.gz and .fastq formats " >&2
            echo $right_file_format >&2
            echo $right_read_file_base_name >&2
            echo $right_read_file >&2

            exit 1
        fi
        gz_right_format=$(echo ${right_read_file} | cut -d'.' -f3)
        if [ $gz_right_format != "gz" ]
        then
        


            printf '%s\n' "right file format invalid: has .fq but no .gz " >&1
            printf '%s\n' "PathSeq only accepts .fq.gz and .fastq formats " >&1
            echo $gz_right_format >&1
            echo $right_read_file_base_name >&1

            printf '%s\n' "right file format invalid: has .fq but no .gz " >&2
            printf '%s\n' "PathSeq only accepts .fq.gz and .fastq formats " >&2
            echo $gz_right_format >&2
            echo $right_read_file_base_name >&2

            exit 1

        fi
        echo "right file format is .fq.gz"
        echo $gz_right_format
        
        
    fi


    echo $right_read_file_base_name
    echo "Created sample unique ID for right read"

    



}











