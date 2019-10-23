#!/bin/bash
ln -s CRYSTAL-calculation/TiTe2_sm200.hessfreq .
i=1
n=12
while ((i <= ${n}))
do
j=1
while ((j <= ${n}))
do
if ((j <= 9))
then
grep -A6 "          ${i}           ${j}" TiTe2_sm200.hessfreq |grep -v "          ${i}           ${j}"  > hessian_${i}-${j}.tmp
sed -e '/^$/d' -e 'N;s/\(.*\)\n\(.*\)/\1\2/' hessian_${i}-${j}.tmp >  hessian_${i}-${j}.dat
rm hessian_${i}-${j}.tmp
else
#echo ${i}    ${j}
grep -A6 "          ${i}          ${j}" TiTe2_sm200.hessfreq | grep -v  "          ${i}          ${j}"   > hessian_${i}-${j}.tmp

sed -e '/^$/d' -e 'N;s/\(.*\)\n\(.*\)/\1\2/' hessian_${i}-${j}.tmp >  hessian_${i}-${j}.dat
rm hessian_${i}-${j}.tmp
fi
((j++))
done
((i++))
done

