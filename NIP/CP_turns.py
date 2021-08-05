import numpy as np

# common constant lists of turn combinations
TURN_DIHEDRALS={'I': (-60, -30, -90, 0), 'II': (-60, 120, 80, 0), \
                'Ip': (60, 30, 90, 0), 'IIp': (60, -120, -80, 0)}

GAMMA_DIHEDRALS={'G': (75, -65), 'iG': (-75, 65)}

# (turn1, turn2) for those with 2 turns
TURN_TYPES=TURN_DIHEDRALS.keys()
TURN_COMBS10=[]
for t1 in TURN_TYPES:
   for t2 in TURN_TYPES:
      if (t1, t2) not in TURN_COMBS10 and (t2, t1) not in TURN_COMBS10: TURN_COMBS10.append((t1, t2))

TURN_TRANS=[]
for t1 in TURN_COMBS10:
   for t2 in TURN_COMBS10:
      if t2 != t1:
         TURN_TRANS.append((t1, t2))

# (turn1, turn2) includes (turn, '-') and ('-', '-')
TURN_TYPES5=TURN_DIHEDRALS.keys()
TURN_TYPES5.append('-')
TURN_COMBS14=[]
for t1 in TURN_TYPES5:
   for t2 in TURN_TYPES5:
      if (t1, t2) not in TURN_COMBS14 and (t2, t1) not in TURN_COMBS14: TURN_COMBS14.append((t1, t2))

XPM_DIHEDRALS={'A': [-60, -30, -90, 0], 'B': [-60, 120, 80, 0], \
                'C': [60, 30, 90, 0], 'D': [60, -120, -80, 0]}

XPM_TURNS={'I': 'A', 'II': 'B', 'Ip': 'C', 'IIp': 'D', '-': '-'}

PSUEDO_XPM_TURNS={'I': 0, 'II': 1, 'Ip': 2, 'IIp': 3, '-': 4, 'G': 5, 'iG': 6}
#PSUEDO_XPM_TURNS={'I': 0, 'II': 1, 'Ip': 2, 'IIp': 3, '-': 4}

# ---------------------------------------------------------------------- #

def get_phi_psi(path):
   """
   reads dihedrals from xvg file
   """
   num=''
   with open(path, 'r') as f:
      line1=f.readline()
      num=len(line1.split())
   frames=np.loadtxt(path, usecols=(0,)) #, comments=['#','@'])
   temp=np.loadtxt(path, usecols=range(2, num))

   frame_phi_psi=dict((frames[i], temp[i]) for i in range(len(frames)))

   return frame_phi_psi

def get_only_phi_psi(path):
   num=''
   with open(path, 'r') as f:
      line1=f.readline()
      num=len(line1.split())
   return np.loadtxt(path, usecols=range(2,num))

# MODIFIED 12-13-15 to make sure once a residue is assigned it doesn't change #
def get_turn_lst(dihed, TOL):
   """
   takes list of all dihedrals as phi/psi in residue order
   returns list of len(residues) with turn classification for each residue
   '-' if residue not classified as a turn
   """
   from math import fabs

   temp=['-' for i in range(len(dihed)/2)]
   for i in range(0, len(dihed)-1, 2):
      res=i/2
      nextres=res+1
      if nextres > len(temp)-1: nextres=0 # makes circular
      if temp[res] == '-' and temp[nextres] == '-':
         # test as i+1 position
         if res == len(dihed)/2-1:
            test=[dihed[i], dihed[i+1], dihed[0], dihed[1]]
         else:
            test=dihed[i:i+4]
         for turn in TURN_DIHEDRALS.keys(): # check all turn types
            match=True
            for k in range(4):
               if fabs(TURN_DIHEDRALS[turn][k]-test[k]) > TOL:
                  match=False
                  break
            if match:
               temp[res]=turn
               if res == len(dihed)/2-1: temp[0]=turn
               else: temp[res+1]=turn
               break
            # if doesn't match any turns as i+1, try as i+2
         prevres=res-1
         if prevres < 0: prevres=len(temp)-1 # makes circular
         if temp[res] == '-' and temp[prevres] == '-':
            if res == 0:
               test=[dihed[-2], dihed[-1], dihed[0], dihed[1]]
            else: test=dihed[i-2:i+2]
            for turn in TURN_DIHEDRALS.keys(): # check all turn types
               match=True
               for k in range(4):
                  if fabs(TURN_DIHEDRALS[turn][k]-test[k]) > TOL:
                     match=False
                     break
               if match:
                  temp[res]=turn
                  if res == 0: temp[-1]=turn
                  else: temp[res-1]=turn
                  break
   return temp

def check_gamma(dihed, TOL, turns):
   """
   checks for gamma or inverse gamma
   should only be used for those already checked for two turns
   that only have 1 turn or lack turns
   returns updated turn list
   """
   from math import fabs
   newTurns=['-' for i in range(len(turns))]
   for i in range(len(turns)):
      if turns[i] != '-': newTurns[i]=turns[i]
      else:
         test=[dihed[i*2], dihed[i*2+1]]
         for gamma in GAMMA_DIHEDRALS.keys():
            match=True
            for j in range(len(test)):
               if fabs(GAMMA_DIHEDRALS[gamma][j]-test[j]) > TOL:
                  match=False
                  break
            if match:
               newTurns[i]=gamma
               break
   return newTurns
         

def get_turns(cp):
   """
   takes list with turn classification for each residue
   NOTE: must have two turns/use check_linkers() first
   returns string of turn1_turn2
   """
   for i in range(len(cp)):
      if cp[i] == '-':
         if i > 0: it1=i-1
         else: it1=len(cp)-1
         if i < len(cp)-1: it2=i+1
         else: it2=0
         return cp[it1]+'_'+cp[it2]

def get_intermed(cp):
   """
   takes list with turn classification for each residue
   NOTE: expects that only 1 turn will be found in either
   """

def check_linkers(turns):
   """
   takes list with turn classifications for each residue
   checks for only 2 residues without turn classifications
   checks if two non-turn residues are next to each other
   NOTE: only really for 6mers
   """
   if turns.count('-') != 2:
      return False
   else:
      for i in range(len(turns)-1):
         if turns[i] == '-' and turns[i+1] == '-':
            return False
      if turns[-1] == '-' and turns[0] == '-': return False
   return True




