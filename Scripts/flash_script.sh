#!/bin/bash
#for f in *R1_trimmed.fastq; do flash $f ${f%R1_trimmed.fastq}R2_trimmed.fastq -x 0.3 -M 250 -o ${f%R1.*} 2>&1 | tee ${f%R1.*}.flash.log;done
#for f in *R1.fastq; do flash $f ${f%R1.fastq}R2.fastq -x 0.5 -M 250 -o ${f%R1.*} 2>&1 | tee ${f%R1.*}.flash.log; done
#for f in *R1_merged.fastq; do flash $f ${f%R1_merged.fastq}R2_merged.fastq -x 0.3 -M 250 -o ${f%R1.*} 2>&1 | tee ${f%R1.*}.flash.log;done
#for f in *_R1_corrected_merged.fastq; do flash $f ${f%_R1_corrected_merged.fastq}_R2_corrected_merged.fastq -x 0.3 -M 250 -O -o ${f%_R1*} 2>&1 | tee ${f%_R1*}.flash.log;done
#for f in *_R1_corrected_merged.fastq; do flash $f ${f%_R1_corrected_merged.fastq}_R2_corrected_merged.fastq -x 0.3 -M 250 -m 40 -O -o ${f%_R1*} 2>&1 | tee ${f%_R1*}.flash.log;done
for f in Output/Flash_merging/*R1.corrected.merged.fastq; do flash $f ${f%R1.corrected.merged.fastq}R2.corrected.merged.fastq -x 0.3 -M 350 -m 40 -O -o ${f%R1*} 2>&1 | tee ${f%R1*}.flash.log;done 