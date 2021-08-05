
IC=8
NC=12
IT=12
NT=16

function CalcEnergySum () {
  for ((it=$IT; it<=$NT; it++ )); do
    for ((ic=$IC; ic<=$NC; ic++)); do
      echo "Calculating energy sum for it=$it, ic=$ic ..."
      python Py_CalcSumEnergy.py $it $ic
    done
  done
}


function CalcMeanEnergy () {
  for it in LJ-14 LJ-SR Coul-14 Coul-SR; do
    for ((ic=$IC; ic<=$NC; ic++)); do
      echo "Calculate mean energies for cluster $ic, energy component = $it ..."
      python Py_CalcMean.py $it $ic
    done
  done
}


function CalcMeanEnergyAll () {
  for ((ic=$IC; ic<=$NC; ic++)); do
    echo "Calculate mean energies for cluster $ic ..."
    python Py_CalcMeanEnergyAll.py $ic
  done
}

CalcEnergySum
CalcMeanEnergy
CalcMeanEnergyAll
