for f in *.xvg; do 
   out=${f/.xvg/.txt}
   tail -n +13 $f > $out
done
