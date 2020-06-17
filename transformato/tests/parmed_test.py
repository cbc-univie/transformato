from parmed.charmm import CharmmParameterSet

print ('std')
params_std = CharmmParameterSet('/home/master/debug/std.prm', '/home/master/debug/std.rtf')
print (params_std.bond_types[('CG331' , 'CG331')])
print (params_std.bond_types[('CG331' , 'HGA3')])

print (params_std.angle_types[('CG331' , 'CG331' , 'HGA3')])
print (params_std.angle_types[('HGA3' , 'CG331' , 'HGA3')])

print (params_std.dihedral_types[('HGA3' , 'CG331' , 'CG331' ,'HGA3')])


print ('non std')
params_non_std = CharmmParameterSet('/home/master/debug/non_std.prm', '/home/master/debug/non_std.rtf')
print (params_non_std.bond_types[('CG331' , 'CG3311')])
print (params_non_std.bond_types[('CG331' , 'HGA3')])
print (params_non_std.bond_types[('CG3311' , 'HGA3')])

print (params_non_std.angle_types[('CG3311' , 'CG331' , 'HGA3')])
print (params_non_std.angle_types[('CG331' , 'CG3311' , 'HGA3')])
print (params_non_std.angle_types[('HGA3' , 'CG331' , 'HGA3')])
print (params_non_std.angle_types[('HGA3' , 'CG3311' , 'HGA3')])

print (params_non_std.dihedral_types[('HGA3' , 'CG331' , 'CG3311' ,'HGA3')])
print (params_non_std.dihedral_types[('HGA3' , 'CG3311' , 'CG331' ,'HGA3')])