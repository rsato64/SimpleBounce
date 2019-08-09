#!/bin/bash

echo "# Eq. 40, 41, 42, 43 in 1906.10829"
echo "# model runtime action"
for name in ../sample1b.x ../sample1c.x ../sample2a.x ../sample2b.x; do
	before=`date +"%s.%N"`
	result=`$name 2>/dev/null | grep "S_E"`
	after=`date +"%s.%N"`
	time=`echo "$after - $before" | bc`
	echo -n "$name  "
	echo -n "$time  "
	echo $result
done
echo

echo "# Table 1 in 1901.03714"
echo "# model runtime action"
for name in ../sample1d.x ../sample2c.x ../sample3.x ../sample4.x ../sample5.x ../sample6.x ../sample7.x ../sample8.x; do
	before=`date +"%s.%N"`
	result=`$name 2>/dev/null | grep "S_E"`
	after=`date +"%s.%N"`
	time=`echo "$after - $before" | bc`
	echo -n "$name  "
	echo -n "$time  "
	echo $result
done

