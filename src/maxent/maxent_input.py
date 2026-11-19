import numpy as np


def read_inputs():
    with open('maxent_input.dat') as stream:
        values = [line.split('#', 1)[0].strip() for line in stream]
    values = [float(value.replace('D', 'e').replace('d', 'e'))
              for value in values if value]
    if len(values) != 10:
        raise ValueError('maxent_input.dat requires 10 parameters')
    ntau, idg, deltag, beta, nomega, domega, steps, alpha, nr, iflat = values
    if iflat not in (0, 1):
        raise ValueError("iflat must be 0 (flat) or 1 (Gaussian)")
    ntau, idg, nomega, steps, nr, iflat = map(int, (ntau, idg, nomega, steps, nr, iflat))
    if min(ntau, nomega, steps) < 1 or min(deltag, beta, domega, alpha) <= 0 or nr < 0:
        raise ValueError('Invalid MaxEnt grid, error, or iteration parameters')
    data = np.loadtxt('Gtau.dat', ndmin=2)
    if data.shape != (ntau, 2) or not np.isfinite(data).all():
        raise ValueError('Gtau.dat must contain Nt finite rows: tau, G(tau)')
    tau, g = data.T
    if np.any(tau < 0) or np.any(tau > beta + 1e-12 * beta):
        raise ValueError('tau must lie in [0, beta]')
    return tau, g, (ntau, idg, deltag, beta, nomega, domega, steps, alpha, nr, iflat)
