#!/usr/bin/env python
# -*- coding: utf8 -*-
import FASMA

if __name__ == '__main__':
    import sys

    if len(sys.argv) > 1:
        cfgfile = sys.argv[1]
    else:
        cfgfile = 'StarMe_synth.cfg'
    driver = FASMA.synthMethod(cfgfile=cfgfile, overwrite=None)
    driver.synthdriver()
