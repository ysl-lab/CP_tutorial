#!/bin/env bash 


function Clean () {
  for inp in *.xvg; do
    out=${inp/xvg/txt}
    echo "Processing $inp --> $out ..."
    tail -n +13 $inp > $out
  done
}

Clean

