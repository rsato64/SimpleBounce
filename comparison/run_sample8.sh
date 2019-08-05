#!/bin/bash
date +"%T.%N"
../sample8.x 2>/dev/null | grep "S_E"
date +"%T.%N"

