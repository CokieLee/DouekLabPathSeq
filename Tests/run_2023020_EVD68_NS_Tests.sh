#!/bin/bash

## read in tests to be performed

## TODO: check if expected output is present

## TODO: if expected output is not present, run necessary portion of code

## for now, pass in known ground truth file and tested file, and run given comparisons
python /hpcdata/vrc/vrc1_data/douek_lab/projects/PathSeq/cokie_pathseq_rewrite/PathSeq/Tests/runTestBowtieUnmasked.py -t /hpcdata/vrc/vrc1_data/douek_lab/projects/PathSeq/cokie_pathseq_rewrite/2023020_EVD68_NS/Sample_582760806_A5_S33_R1_001/ -i /hpcdata/vrc/vrc1_data/douek_lab/projects/PathSeq/cokie_pathseq_rewrite/2023020_EVD68_NS_TestRun/Sample_582760806_A5_S33_R1_001/
