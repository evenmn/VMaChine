import mpi4py
import vmc


rng = vmc.MersenneTwister()

print(rng)
print(rng.next_int(10))
print(rng.next_double())
print(rng.next_gaussian(0, 1))
