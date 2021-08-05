import sys

# would be beter with module psutils... couldn't easily install properly
def memory_usage_resource():
   """
   returns current usage in MB
   note: doesn't seem to see that orphaned arrays are liberated by the python interpreter
   source: http://fa.bianp.net/blog/2013/different-ways-to-get-memory-consumption-or-lessons-learned-from-memory_profiler/
   """
   import resource
   rusage_denom = 1024.
   if sys.platform == 'darwin':
      rusage_denom=rusage_denom*rusage_denom
   mem=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/rusage_denom
   return mem

# ------------ file management ----------- #

def make_dir(name):
   """
   makes dir
   if fails asks if want to exit
   """
   import commands
   cmd='mkdir '+name
   status, output=commands.getstatusoutput(cmd)
   if status != 0:
      print cmd+' failed'
      cont=raw_input('continue anyways? (yes/no) ')
      if cont == 'no': sys.exit()
      
def mv_file(directory, name):
   """
   moves file name to given dir
   * NOTE: does not exit if fails *
   """
   import commands
   cmd='mv '+name+' '+directory
   commands.getoutput(cmd)

def get_basename(pwd_name):
   """
   returns base file or dir of given path
   file name returned with file extension
   given path can have a / at the end; will be ignored
   """
   import os
   if pwd_name[-1] == '/': return os.path.basename(pwd_name[:-1])
   else: return os.path.basename(pwd_name)

def get_files_from_dir(dirname):
   """
   equivalent to 'ls'
   returns list of only files in given dir
   """
   import commands, os.path
   cmd='ls '+dirname
   output=commands.getoutput(cmd)

   fils=[f for f in output.split('\n') if os.path.isfile(dirname+'/'+f)]
   return fils

def get_dirs_from_dir(dirname):
   """
   returns list of subdirs in given dirname
   """
   import os
   sub=os.listdir(dirname)
   #for d in sub:
    #  print d, os.path.isdir(dirname+'/'+d)
   dirs=[d for d in sub if os.path.isdir(dirname+'/'+d)]
   return dirs

def rename_file(curr, new):
   """
   renames a file; uses 'mv'
   """
   import commands
   cmd='mv '+curr+' '+new
   commands.getoutput(cmd)

def get_lines(f_name):
   """
   opens, reads, returns lines
   """
   f=open(f_name)
   lines=f.readlines()
   f.close()
   return lines

def determine_version(name):
   """
   funciton for mk_outputdir to determine verison number
   modified 10-8-15 to create versions > 10
   """
   import commands
   cmd='ls'
   output=commands.getoutput(cmd)
   progs=[p for p in output.split('\n') if p[:len(name)] == name and p[len(name):len(name)+2] == '_v']
   versions=[]
   for p in progs:
      curr=p[-1]
      for i in range(len(p)-2, 0, -1):
         if is_number(p[i]): curr+=p[i]
         else:
            versions.append(int(curr[::-1]))
            break
   #versions=[int(p[-1]) for p in progs]
   if len(versions) != 0: return max(versions)+1
   else: return 1

def mk_outputdir(temp_name):
   """
   assumes program name is base of temp_name given
   makes dir intended for output files of given program name
   determines version num based on dirs in current dir
   calls determine_version, make_dir
   """
   import re
   match=re.search(r'\w+', get_basename(temp_name))
   prog_name=match.group()
   v=determine_version('out_'+prog_name)
   name='out_'+prog_name+'_v'+str(v)
   make_dir(name)
   return name

# ------------ for info files ------------ #

def determine_file_version(name):
   """
   function for mk_saved_file to determine version
   * for .txt files *
   modified 1-4-16 to create versions > 10
   """
   #print name
   import commands
   cmd='ls'
   output=commands.getoutput(cmd)
   #files=[f for f in output.split('\n') if f[:len(name)] == name]
   files=[f for f in output.split('\n') if f[:len(name)] == name and f[len(name):len(name)+2] == '_v']
   #print files
   if len(files) == 1: return 1
   versions=[]
   for f in files:
      curr=f[-5]
      for i in range(len(f)-6, 0, -1):
         if is_number(f[i]): curr+=f[i]
         else:
            if is_number(curr): versions.append(int(curr[::-1]))
            break
   #print versions
   if len(versions) != 0: return max(versions)+1
   else: return 1
   #else:
    #  versions=[f[-5] for f in files if len(f) > len(name)+4]
     # versions=[int(f[-5]) for f in files if len(f) > len(name)+4]
   #return max(versions)+1

def mk_saved_file(temp_name):
   """
   moves given file to file_name_v#.txt
   oldest files are lowest versions
   """
   import re
   match=re.search(r'\w+', get_basename(temp_name))
   fil_name=match.group()
   #print temp_name, fil_name
   v=determine_file_version(fil_name)
   name=fil_name+'_v'+str(v)+'.txt'
   rename_file(temp_name, name)
   return

def split_by_pdb(lines):
   """
   for info files, splits entries into dict with key as pdbid
   """
   pdb_lines={}
   pdb=''
   for l in lines:
      if l[:2] == 'ID':
         pdb=l[10:14]
         if pdb in pdb_lines.keys():
            pdb_lines[pdb].append(l)
         else:
            pdb_lines[pdb]=[l]
      else: pdb_lines[pdb].append(l)
   return pdb_lines

def split_by_site(lines):
   """
   for info files, splits by site/entry as [[lines/site], ...]
   """
   sites=[]
   temp=[]
   for l in lines:
      if l[:2] == '\\\\':
         temp.append(l)
         sites.append(temp)
         temp=[]
      else:
         temp.append(l)
   return sites

# ------------ generally useful ---------- #

def make_combinations(key_lst):
   """
   from values in list, makes all possible combinations
   of all possible lengths; order does not matter
   returns list of created strings
   """
   import itertools
   lst_combs=[]
   for i in xrange(1, len(key_lst)+1):
      temp=[list(x) for x in itertools.combinations(key_lst, i)]
      lst_combs.extend(temp)

   combs=[]
   for c in lst_combs:
      temp=''
      for i in range(len(c)):
         temp+=c[i]
      combs.append(temp)
   return combs

def is_number(s):
   """
   determines if string is number
   """
   try:
      float(s)
      return True
   except ValueError:
      return False

# ------------ for lists ----------------- #

def extend_lst(lst, num):
   """
   extends given lst by num
   """
   for i in range(len(lst), len(lst)+num+1):
      lst.append(0)

   return

def trim_lst(lst):
   """
   trims given list to remove indices on end that contain 0's
   """
   for i in range(len(lst)-1, -1, -1):
      if i == len(lst)-1 and lst[i] != 0: return
      elif lst[i] != 0:
         del lst[i+1:len(lst)]
         return

def shuffle_lst(lst):
   """
   shuffles the elements of the given list
   also works for sizes over 2080 which random.shuffle cannot handle
   * not sure how random?? efficient?? *
   """
   import random
   length=len(lst)
   copy=[i for i in lst]

   new=[]
   for i in range(length):
      new.append(random.choice(copy))
      copy.remove(new[-1])
   return new

# ------------------ for typical 2D array like output ----------- #

def get_data_2D(file_name):
   """
   takes file name (will open and read)
   to read in 2D array like data ouput
   returns 2D list of floats: [[values corresponding to index2], ...] with len(lst) = len(index1)
   returns list for index1 (for each output row), index2 (for each column)
   """
   full=get_lines(file_name)
   temp=full[1:]

   index2=[i.strip() for i in full[0].split('\t') if i != '']
   index1=[l.split('\t')[0] for l in temp]

   data=[[] for i in range(len(temp))]
   for i in range(len(temp)):
      for j in temp[i].split('\t')[1:]:
         data[i].append(float(j.strip()))

   return data, index1, index2

def convert_lst_dict2(lst, index1, index2):
   """
   takes 2D list, returns 2D dictionary
   index1 is correctly ordered list of key1; index2 is correctly ordered list of key2
   """
   new=dict((i, dict((j, lst[index1.index(i)][index2.index(j)]) for j in index2)) for i in index1)
   return new      




