#!/bin/bash
echo Cleaning python code in \'`pwd`\' and subdirectories
for file in `find quickff scripts *.py | egrep '(\.py$)|(\.c$)|(\.h$)|(\.pxd$)|(\.pyx$)|(\.rst$)|(^scripts/)'`; do
  echo "Cleaning $file"
  sed -i -e $'s/\t/    /' ${file}
  sed -i -e $'s/[ \t]\+$//' ${file}
  #sed -i -e $'s/^# --$/#--/' ${file}
  #sed -i -e $'s/^\/\/ --$/\/\/--/' ${file}
  sed -i -e :a -e '/^\n*$/{$d;N;ba' -e '}' ${file}
  #./updateheaders.py ${file}
done
exit 0
