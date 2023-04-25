# Analysis Plan
This repo is a templete for DCC Differential analysis(DEseq2).


## Target
No specific target


## Data source
Example data: /restricted/projectnb/ncrna/minty/samb_data4
Total 4 groups are there.

## Step 1 DCC step
Following the pipeline: ```https://github.com/90yearsoldcoder/DCC_protocol.git```

## Step 2 Prepare DCC result data
* Script: ```Script/concatenate.py```
* Manually reset the input parameter in the python file
* inputdata: DCCresult/CircRNACount
* outputdata: <group_name>_counts.tsv

## Step 3 Prepare the group information
* Example group information: ```Processing/group_info.tsv```
 
| sampleID                   | group    |
|----------------------------|----------|
| WB_LJ_PS_EX4_83_S83_L001  | Control1 |
| WB_LJ_PS_EX4_84_S84_L001  | Control1 |
| WB_LJ_PS_EX4_85_S85_L001  | Control1 |
| WB_LJ_PS_EX4_86_S86_L001  | otao1    |
| WB_LJ_PS_EX4_87_S87_L001  | otao1    |
| WB_LJ_PS_EX4_88_S88_L001  | otao1    |
| WB_LJ_PS_EX4_89_S89_L001  | Control2 |
| WB_LJ_PS_EX4_90_S90_L001  | Control2 |
| WB_LJ_PS_EX4_91_S91_L001  | Control2 |
| WB_LJ_PS_EX4_92_S92_L001  | otao2    |
| WB_LJ_PS_EX4_93_S93_L001  | otao2    |
| WB_LJ_PS_EX4_94_S94_L001  | otao2    |


## Step 4 DEseq2

....


## Expected Result

## Plot
