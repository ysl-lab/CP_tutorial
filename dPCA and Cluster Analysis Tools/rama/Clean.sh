for f in */*.xvg; do 
   out=${f/.xvg/.txt}
   tail -n +18 $f > $out
done
