#!/bin/bash -l

perl -MCPAN -Mlocal::lib -e 'CPAN::install('$1')'
