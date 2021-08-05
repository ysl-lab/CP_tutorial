#!/bin/env bash

ff=rsff2
sol=tip3p
prot=s2c
NX=50
NY=50
NZ=50
input=${prot}.txt
TITLE="RSFF2 + TIP3P, BE-2D-META, S2 50-100ns, ${NX}x${NY}x${NZ} grids"

XMIN=
XMAX=
YMIN=
YMAX=
ZMIN=
ZMAX=
DEN_MAX=5
FES_MAX=2.25

clean=0.05
##########################################################################################
function calc_den () {
  den=${input/.txt/.den}
  echo "Calcuting 2D density profile using input data: $input..."
  #python Py_CalcDensit2D_v3.py $input $XMIN $XMAX $YMIN $YMAX $ZMIN $ZMAX $NX $NY $NZ > $den
  python Py_CalcDensit3D.py $input $XMIN $XMAX $YMIN $YMAX $ZMIN $ZMAX $NX $NY $NZ > $den

  den_png=${input/.txt/_den.png}
  echo "Ploting the 2D Density using input data: $den ..."
  title="$TITLE, 2D density"
  #gnuplot -e "TITLE='$title'; INPUT='$den'; XMIN='$XMIN'; XMAX='$XMAX'; YMIN='$YMIN'; YMAX='$YMAX'; ZMIN='$ZMIN'; ZMAX='$ZMAX'; DEN_MAX='$DEN_MAX'" Gp_PlotDensit2D_v1.gplt
  gnuplot -e "TITLE='$title'; INPUT='$den'; XMIN='$XMIN'; XMAX='$XMAX'; YMIN='$YMIN'; YMAX='$YMAX'; ZMIN='$ZMIN'; ZMAX='$ZMAX'; DEN_MAX='$DEN_MAX'" Gp_PlotDensit3D.gplt
  convert -density 300 tmp.eps $den_png
}

##########################################################################################
function clean() {
  # Remove grids with density lower than pre-defined cutoff 
  den=${input/.txt/.den}
  den2=${input/.txt/_kept.den}
  echo "Removing grids with 0 density from data file: $den .."
  python Py_RemoveLowDensitGrids.py $den $clean > $den2

  den_png=${input/.txt/_keptden.png}
  echo "Ploting the 2D Density using input data: $den ..."
  title="$TITLE, 2D density"
  gnuplot -e "TITLE='$title'; INPUT='$den2'; XMIN='$XMIN'; XMAX='$XMAX'; YMIN='$YMIN'; YMAX='$YMAX'; ZMIN='$ZMIN'; ZMAX='$ZMAX'; DEN_MAX='$DEN_MAX'" Gp_PlotKeptDensit3D.gplt
  convert -density 300 tmp.eps $den_png


  # Construct the distance matrix with kept grids 
  echo "Constructing distance matrix using input data: $den2 ..."
  dmtx=${input/.txt/_kept.dmtx}
  python Py_CalcDistMatrixWithDensit.py $den2 > $dmtx 
}

##########################################################################################
function calc_pop () {
  # Combine the Assignment file with density grids
  den2=${input/.txt/_kept.den}
  python Py_CombineDensitAndAssign.py $den2 CLUSTER_ASSIGNATION > GRID_ASSIGNATION
  # Plot the grid assignment without halo
  gcl_png=${input/.txt/_grid_cluster_v1.png}
  title="$TITLE, grid cluster"
  #gnuplot -e "TITLE='$title'; INPUT='GRID_ASSIGNATION'; XMIN='$XMIN'; XMAX='$XMAX'; YMIN='$YMIN'; YMAX='$YMAX'" Gp_PlotGridCluster2D_v1.gplt
  gnuplot -e "TITLE='$title'; INPUT='GRID_ASSIGNATION'; XMIN='$XMIN'; XMAX='$XMAX'; YMIN='$YMIN'; YMAX='$YMAX'; ZMIN='$ZMIN'; ZMAX='$ZMAX';" Gp_PlotGridCluster3D.gplt
  convert -density 300 tmp.eps $gcl_png
  
  # Calculate the population
  python Py_CalcPopFromGrid_v1.py $XMIN $XMAX $YMIN $YMAX $ZMIN $ZMAX $NX $NY $NZ
  
  # Plot the grid assignment according to the population
  gcl_png2=${input/.txt/_grid_cluster_v2.png}
  title="$TITLE, grid cluster"
  gnuplot -e "TITLE='$title'; INPUT='GRID_ASSIGNATION2'; XMIN='$XMIN'; XMAX='$XMAX'; YMIN='$YMIN'; YMAX='$YMAX'; ZMIN='$ZMIN'; ZMAX='$ZMAX';" Gp_PlotGridCluster3D.gplt
  convert -density 300 tmp.eps $gcl_png2
}
##########################################################################################
calc_den
clean
# Do cluster analysis using the kept grids
#calc_pop


