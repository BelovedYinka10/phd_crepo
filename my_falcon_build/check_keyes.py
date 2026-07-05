
from keygen import keygen

sk, pk = keygen(512)
print('Secret key f (first 10 coeffs):', sk.f[:10])
print('Secret key g (first 10 coeffs):', sk.g[:10])
print('Public key h (first 10 coeffs):', pk.h[:10])
print('Public key size:', len(pk.h), 'coefficients')


