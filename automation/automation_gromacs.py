import os
import shutil

os.system('clear')

copy1 = shutil.copy('../00_itp/simu_sysm_step.mdp', 'simu_sysm_step.mdp')
copy2 = shutil.copy('../00_itp/simu_sysm_nvtt.mdp', 'simu_sysm_nvtt.mdp')
copy3 = shutil.copy('../00_itp/simu_sysm_nptt.mdp', 'simu_sysm_nptt.mdp')
copy4 = shutil.copy('../00_itp/simu_sysm_stat.mdp', 'simu_sysm_stat.mdp')

######################################################## pbs file

os.system('qhost')

print('')
inode = input('---------- Which node to use? ')
print('you input is ' + inode)

ipath = os.path.dirname(os.path.realpath(__file__))

fin = open("../00_itp/01_mmk_grogpu_gromacs.pbs", "rt")
data = fin.read()
data = data.replace('xnodex', inode)
data = data.replace('xpathx', ipath)
fin.close()

fout = open("01_mmk_grogpu_gromacs.pbs", "wt")
fout.write(data)
fout.close()


######################################################## packmol file

fout = open("packmol.inp", "wt")

head = '''tolerance 2.5
filetype pdb
output packmol.pdb
'''
fout.write(head)

print('')
imolecules = input('---------- How many ions/molecules? ')
idlxyz = input('---------- Size of modelling box: ')
print('There are ' + imolecules  + ' type molecules in (' + idlxyz + 'A)^3')
print('')

structure = '''
structure ../00_itp/{0}/{0}.pdb
  number {1}
  inside box 0.0 0.0 0.0 {2} {2} {2}
end structure
'''

iname = []
imoll = []

for inow in range(int(imolecules)):
    print('')
    print('Name and No. of molecules for Type ', inow+1)
    iaaaa = input('---------- Name of molecules: ')
    ibbbb = input('---------- No. of molecules: ')
    fout.write(structure.format(iaaaa, ibbbb, idlxyz))
    iname.append(iaaaa)
    imoll.append(ibbbb)
fout.close()

######################################################## top file

fout = open("simu_sysm_topp.top", "wt")

head = '''[ defaults ]
; nbfunc  comb-rule  gen-pairs  fudgeLJ  fudgeQQ
  1       2          yes        0.50     0.50

'''
fout.write(head)

inc_sigepi = '#include "../00_itp/{0}/{0}_sigepi.itp"'
inc_normal = '#include "../00_itp/{0}/{0}.itp"'

for inow in range(int(imolecules)):
    fout.write(inc_sigepi.format(iname[inow]) + '\n')

fout.write('\n')

for inow in range(int(imolecules)):
    fout.write(inc_normal.format(iname[inow]) + '\n')

simu_sysm = '''
[ system ]
simu_sysm

[ molecules ]
'''
fout.write(simu_sysm)

for inow in range(int(imolecules)):
    fout.write(iname[inow] + ' ' + imoll[inow] + '\n')
fout.close()

######################################################## run packmol

print('')
ipackmol = input('---------- Run packmol to generte pdb file? ')
if ipackmol[0] == 'y':
    os.system('packmol < packmol.inp')
    os.system('editconf -f packmol.pdb -o simu_sysm_pkml.gro -d 0')

######################################################## run packmol

print('')
print('----- The End =====')
print('')
