#!/bin/bash
for i in $(find quickff | egrep "\.pyc$|\.py~$|\.pyc~$|\.bak$|\.so$") ; do rm -v ${i}; done
for i in $(find quickff/tests | egrep "\.pyc$|\.py~$|\.pyc~$|\.bak$|\.so$") ; do rm -v ${i}; done
for i in $(find doc | egrep "\.pyc$|\.py~$|\.pyc~$|\.bak$|\.so$") ; do rm -v ${i}; done
for i in $(find share/systems | egrep "\.pdf$|\.xyz$|\.txt$|\.chk$|\.pps$|\.zip$|\.log$") ; do rm -v ${i}; done

(cd doc; make clean)
rm -vr dist
rm -vr build
