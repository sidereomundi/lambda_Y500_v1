#!/bin/bash 
bfile=tmp.pro
echo .r "call_likelihood_library.pro" >> $bfile
echo debugt >> $bfile
echo exit >> $bfile
idl << EOF > out.txt
@$bfile
EOF
rm $bfile
cat out.txt