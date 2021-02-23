#!/bin/bash

SCRIPT_HOME=$(dirname $(readlink -f $0))
shopt -s nullglob

# Number of sequences in each BLAST job
MAX_BLAST=10

source ${HPC_DIBIGTOOLS_DIR}/lib/sh/utils.sh

CMD=$1

# if [[ "$CMD" == "makedbOLD" ]];
# then
#   FASTA=$2
#   DBDIR=$3
#   shift 3
#   REMOVE=$*

#   mkdir -p $DBDIR
#   FNAME=$(basename $FASTA)
#   NEWFASTA=${DBDIR}/${FNAME%.*}.sub.fa
#   echo "Cleaning alignment $FASTA"
#   if [[ "$REMOVE" != "" ]]; then
#     echo "Removing sequences containing \"$REMOVE\""
#   fi
#   ${SCRIPT_HOME}/flacotools.py clean $FASTA $REMOVE > $NEWFASTA
#   echo "Clean alignment written to $NEWFASTA "
#   echo "Generating BLAST database"
#   submit -done makedb.@.done ${SCRIPT_HOME}/makedb.qsub $NEWFASTA
#   wait_for makedb 1
#   echo "BLAST database generated."
#   echo $NEWFASTA
#   exit 0
# fi

if [[ "$CMD" == "makedb" ]];
then
  FASTA=$2
  shift 2
  REMOVE=$*

  mkdir -p split-by-month
  FNAME=$(basename $FASTA)
  CLEAN=${FNAME%.*}.clean.fa
  EXCLUDED=${FNAME%.*}.excl.fa
  echo "Cleaning alignment $FASTA"
  if [[ "$REMOVE" != "" ]]; then
    echo "Removing sequences containing \"$REMOVE\""
  fi
  ${SCRIPT_HOME}/flacotools.py splitdb $FASTA $CLEAN split-by-month $EXCLUDED $REMOVE

  nbl=0
  echo "Generating BLAST database"
  for db in $(find split-by-month -name DB); do
      submit -done makedb.@.done ${SCRIPT_HOME}/makedb.qsub $db
      nbl=$((nbl+1))
  done
  wait_for makedb $nbl
  exit 0
fi

if [[ "$CMD" == "split" ]];
then
  FASTA=$2
  mkdir -p split-by-month
  echo "Splitting file $FASTA by month..."
  ${SCRIPT_HOME}/flacotools.py split $FASTA split-by-month
  exit 0
fi

if [[ "$CMD" == "blast" ]];
then
  INDEX=DB
  echo "Starting BLAST jobs..."
  nj=0
  rm -f blast.*.done
  for spdir in split-by-month/2*/;
  do
      blastcmd=$spdir/blast.cmd
      find $spdir -name \*.fa | xargs -n $MAX_BLAST echo submit -done blast.@.done -o --account=salemi,--qos=salemi-b ${SCRIPT_HOME}/blast.qsub $spdir/${INDEX} > $blastcmd
      nfastas=$(grep -c ^ $blastcmd)
      source $blastcmd
      # submit -done blast.@.done ${SCRIPT_HOME}/blast.qsub $INDEX $spdir
      nj=$((nj+nfastas))
  done
  wait_for blast $nj
  exit 0
fi

if [[ "$CMD" == "parse" ]];
then
  OUTFILE=$2
  MATCHES=$3
  ${SCRIPT_HOME}/flacotools.py parse split-by-month $OUTFILE $MATCHES
  exit 0
fi

if [[ "$CMD" == "extract" ]];
then
  FASTA=$2
  MATCHES=$3
  OUTFILE=$4
  ${SCRIPT_HOME}/flacotools.py extract $FASTA $MATCHES $OUTFILE
  exit 0
fi

if [[ "$CMD" == "run" ]];
then
  R_TARGETFASTA=$2
  R_DBFASTA=$3
  R_OUTFILE=${R_TARGETFASTA%.*}.out.txt
  R_MATCHES=${R_TARGETFASTA%.*}.matches.txt
  R_OUTFASTA=${R_TARGETFASTA%.*}.out.fa
  $0 split $R_TARGETFASTA
  $0 blast
  $0 parse $R_OUTFILE $R_MATCHES
  $0 extract $R_DBFASTA $R_MATCHES $R_OUTFASTA
  exit 0
fi

echo "flaco_blast.sh - BLAST pipeline for GISAID sequences."
echo
echo "Usage: flaco_blast.sh command [options...]"
echo
echo "where command is one of: run, makedb, split, blast, parse, extract."
echo
