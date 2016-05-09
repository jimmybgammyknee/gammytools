#!/bin/bash -l

vcf=$1

bgzip $vcf; tabix -p vcf $vcf.gz
