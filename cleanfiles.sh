#!/bin/bash
for i in $(find quickff | egrep "\.pyc$|\.py~$|\.pyc~$|\.bak$|\.so$") ; do rm -v ${i}; done

(cd doc; make clean)
rm -v MANIFEST
rm -vr dist
rm -vr build
