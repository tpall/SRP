[![Build Status](https://travis-ci.org/tpall/SRP.svg?branch=master)](https://travis-ci.org/tpall/SRP)

# SRP

Calculate retrospective power of your omics experiment.

'''
## Import pvalues
library(qvalue)
data("hedenfalk")
pvalues <- hedenfalk$p

## calculate SRP
library(SRP)
pw <- srp(pvalues, FDR = 0.05)
pw
'''

