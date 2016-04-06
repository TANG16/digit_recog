#!/bin/bash
# --------------------------------------------------
#
##### JOB NAME
#$ -N script_nn
#
##### EMAIL AT INITIATION/COMPLETION/ABORTION/SUSPENSION
#$ -m beas
#
##### LLE-HSI DUKECLUSTER WORKING DIRECTORY
#$ -wd /home/home1/tl173/digit_recog/CNN
#
##### SPECIFY THAT THIS JOB IS SCRIPT-BASED
#$ -b y
#
##### RUN IN THE 'COMPSCI' QUEUE
#$ -q compsci
#
##### ARRAY JOB INDICES
#$ -t 1-1:1
#
##### MEMORY REQUIREMENTS
#$ -l mem_free=90G
#
##### THREAD REQUIREMENTS [INACTIVE]
# (-pe threaded 8-48)
#
##### OUTPUT/ERROR LOG PATHS
#$ -e logs-qsub/$JOB_NAME.e$JOB_ID.$TASK_ID
#$ -o logs-qsub/$JOB_NAME.o$JOB_ID.$TASK_ID
#
# --------------------------------------------------


##################################################
### PARAMETER ARRAYS


OUTPATH="/usr/xtmp/tl173/lle-hsi/workspace-logs/DIGITS_NN"


##################################################
### MATLAB SCRIPT INVOCATION

# bash arrays are 0-indexed
IDX=$(( SGE_TASK_ID -1 ))


# -------------------- start the MATLAB script --------------------
matlab -nodisplay \
    -r "script_nn $OUTPATH"
# -----------------------------------------------------------------


# bye-bye
exit 0



# ------------------------------------------------------------
#
# Alexandros-Stavros Iliopoulos		ailiop@cs.duke.edu
#
# April 16, 2015
#
# ------------------------------------------------------------
