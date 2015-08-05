#!/bin/bash -l

# bgzip vcf and create tabix index

bgzip $1; tabix -p vcf ${1}.gz
