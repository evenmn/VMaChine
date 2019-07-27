from pyscf import gto, scf
mol = gto.M(atom='he 0 0 0')
mol.basis = '4-31g'
mol.build()

mf = scf.RHF(mol)
mf.kernel()
print(mf.mo_coeff.shape)
print(mf.mo_coeff)
