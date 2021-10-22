
def get_psi(num_res, resids, residues, atoms, indexes):
    psi = []
    for i in range(1, num_res + 1):
        if i == num_res:
            next_res = 1
        else: 
            next_res = i + 1
        temp = []
        for n in range(len(residues)):
            if resids[n] == i and atoms[n] == 'N':
                temp.append(indexes[n])
        for n in range(len(residues)):
            if resids[n] == i and atoms[n] == 'CA':
                temp.append(indexes[n])
        for n in range(len(residues)):
            if resids[n] == i and atoms[n] == 'C': 
                temp.append(indexes[n])
        for n in range(len(residues)):
            if resids[n] == next_res and atoms[n] == 'N':
                temp.append(indexes[n])
        psi.append(temp)
    return psi


def get_phi(num_res, resids, residues, atoms, indexes):
    phi = []
    for i in range(1, num_res + 1):
        if i == 1:
            prev_res = num_res
        else:
            prev_res = i - 1
        temp = []
        for n in range(len(residues)):
            if resids[n] == prev_res and atoms[n] == 'C':
                temp.append(indexes[n])
        for n in range(len(residues)):
            if resids[n] == i and atoms[n] == 'N':
                temp.append(indexes[n])
        for n in range(len(residues)):
            if resids[n] == i and atoms[n] == 'CA':
                temp.append(indexes[n])
        for n in range(len(residues)):
            if resids[n] == i and atoms[n] == 'C':
                temp.append(indexes[n])
        phi.append(temp)
    return phi

def get_phiP(num_res,resids,residues,atoms,indexes):
   phiP=[]
   for i in range(1,num_res+1):
      temp=[]
      for n in range(len(residues)):
         if resids[n] == i and residues[n] == 'ALA' and atoms[n] == 'H': temp.append(indexes[n])
         elif resids[n] == i and residues[n] == 'NMA' and atoms[n] == 'CN': temp.append(indexes[n])
      for n in range(len(residues)):
         if resids[n] == i and atoms[n] == 'N': temp.append(indexes[n])
      for n in range(len(residues)):
         if resids[n] == i and atoms[n] == 'CA': temp.append(indexes[n])
      for n in range(len(residues)):
         if resids[n] == i and atoms[n] == 'C': temp.append(indexes[n])
      phiP.append(temp)
   return phiP

def get_psiP(num_res,resids,residues,atoms,indexes):
   psiP=[]
   for i in range(1,num_res+1):
      temp=[]
      for n in range(len(residues)):
         if resids[n] == i and atoms[n] == 'N': temp.append(indexes[n])
      for n in range(len(residues)):
         if resids[n] == i and atoms[n] == 'CA': temp.append(indexes[n])
      for n in range(len(residues)):
         if resids[n] == i and atoms[n] == 'C': temp.append(indexes[n])
      for n in range(len(residues)):
         if resids[n] == i and atoms[n] == 'O': temp.append(indexes[n])
      psiP.append(temp)
   return psiP

def get_omega(num_res,resids,residues,atoms,indexes):
   NMA=[]
   for i in range(len(residues)):
      if residues[i] == 'NMA' and resids[i] not in NMA: NMA.append(resids[i])

   omega=[]
   #for i in range(1,num_res+1):
   for i in NMA:
      if i == 1: prev_res=num_res
      else: prev_res=i-1
      temp=[]
      for n in range(len(residues)):
         if resids[n] == prev_res and atoms[n] == 'CA': temp.append(indexes[n])
      for n in range(len(residues)):
         if resids[n] == prev_res and atoms[n] == 'C': temp.append(indexes[n])
      for n in range(len(residues)):
         if resids[n] == i and atoms[n] == 'N': temp.append(indexes[n])
      for n in range(len(residues)):
         if resids[n] == i and atoms[n] == 'CA': temp.append(indexes[n])
      omega.append(temp)
   #print omega
   return omega,NMA