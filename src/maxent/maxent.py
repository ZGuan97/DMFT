#!/usr/bin/env python3

import numpy as np
from scipy.linalg import eigvalsh, inv
from scipy.optimize import minimize
from scipy.special import logsumexp

from maxent_input import read_inputs


try:
    from numba import njit
except ImportError:
    def njit(func):
        return func


@njit
def seed_rng(seed):
    np.random.seed(seed)


@njit
def rand():
    return np.random.random()


@njit
def entropy_term(value, reference, domega):
    if value == 0.0:
        return 0.0
    return -domega * value * np.log(value / reference)


@njit
def entropy(spectrum, model, domega):
    value = 0.0
    for j in range(len(spectrum)):
        value += entropy_term(spectrum[j], model[j], domega)
    return value


@njit
def maxent(spectrum, rfactor, alpha, temp, kernel, weight, target, model, domega, steps):
    am = np.max(spectrum)
    delta_a = (spectrum / 10.0 + 0.01) * rfactor
    n_acc = 0
    n_trial = 0
    for step in range(1, steps + 1):
        fit = kernel @ spectrum
        chi1 = np.sum(weight * (target - fit)**2)
        for trial in range(len(spectrum)):
            while True:
                j = np.random.randint(len(spectrum))
                k = np.random.randint(len(spectrum) - 1)
                if k >= j:
                    k += 1
                while True:
                    dj = delta_a[j] * (rand() - 0.5)
                    if spectrum[j] + dj >= 0:
                        break
                dk = -dj + 0.05 * delta_a[k] * (rand() - 0.5)
                if spectrum[k] + dk >= 0:
                    break
            if spectrum[j] > 0.1 * am:
                n_trial += 1
            newfit = fit + dj * kernel[:, j] + dk * kernel[:, k]
            chi2 = np.sum(weight * (target - newfit)**2)
            aj_new, ak_new = spectrum[j] + dj, spectrum[k] + dk
            ds = (entropy_term(aj_new, model[j], domega)
                  - entropy_term(spectrum[j], model[j], domega)
                  + entropy_term(ak_new, model[k], domega)
                  - entropy_term(spectrum[k], model[k], domega))
            delta_Q = ((chi1 - chi2) / 2 + alpha * ds) / temp

            if rand() < np.exp(delta_Q):
                if spectrum[j] > 0.1 * am:
                    n_acc += 1
                spectrum[j], spectrum[k] = aj_new, ak_new
                fit = newfit
                chi1 = chi2
        if step % 100 == 0:
            rate = n_acc / n_trial if n_trial else np.nan
            if rate > 0.1:
                if rfactor < 0.01:
                    rfactor *= 1.5
            else:
                if rfactor > 0.001:
                    rfactor /= 1.5
            am = np.max(spectrum)
            delta_a = (0.1 * spectrum + 0.01 * am) * rfactor
            if step % (steps // 10) == 0:
                print('ann-step=', step, 'chi2=', chi1, 'accepted=', n_acc)
            temp /= 1.5
            n_acc, n_trial = 0, 0
    return spectrum, rfactor, temp


def solve_anneal(omega, kernel, target, sigma, model, domega, steps, alpha, nr, seed=0):
    seed_rng(seed)
    spectrum = np.array([rand() for j in range(len(omega))])
    spectrum *= 1 / (spectrum.sum() * domega)
    weight = 1 / sigma**2 # C^-1
    H = (kernel.T * weight) @ kernel
    temp, rfactor = 10.0, 1.0
    while True:
        print('Restarting maxent with rfactor=', rfactor, 'alpha=', alpha, flush=True)
        spectrum, rfactor, temp = maxent(spectrum, rfactor, alpha, temp, kernel, weight, target, model, domega, steps)
        s = entropy(spectrum, model, domega)
        scale = np.sqrt(spectrum / domega)
        Lambda = scale[:, None] * H * scale[None, :]
        trace = np.sum(Lambda * inv(Lambda + alpha * np.eye(len(spectrum))).T)
        ls = -2 * s * alpha
        ratio = ls / trace
        print('Finished maxent with alpha=', alpha, '-2*alpha*S=', ls, 'Trace=', trace, 'ratio=', ratio, flush=True)
        temp, rfactor = 0.001, 0.05
        if abs(ratio - 1) < 0.001:
            break
        elif trace / ls < 0.05:
            alpha *= 0.05
        else:
            alpha *= (trace / ls) * (1 + 0.001 * (rand() - 0.5))
    for st in range(nr):
        print('Starting smoothing', st + 1, flush=True)
        spectrum = np.array([spectrum[max(0, j - 3):j + 4].mean() for j in range(len(spectrum))])
        spectrum *= 1 / (domega * spectrum.sum())
        spectrum, rfactor, temp = maxent(spectrum, 0.005, alpha, 0.005, kernel, weight, target, model, domega, steps)
    return spectrum


def solve_alpha(kernel, target, sigma, model, domega, alpha, maxiter=1000):
    """Return the fixed-alpha optimum and a primal-dual certificate.

    p=domega*A, q=domega*m, F=K/(sigma*domega), y=target/sigma.
    Primal: min .5*||F p-y||^2 + alpha*KL(p,q), sum(p)=1.
    Negative dual: .5*z.z + y.z + alpha*logsumexp(log(q)-F.T*z/alpha).
    """
    f = kernel / (sigma[:, None] * domega)
    y = target / sigma
    logq = np.log(model * domega)

    def terms(z):
        logits = logq - f.T @ z / alpha
        lognorm = logsumexp(logits)
        p = np.exp(logits - lognorm)
        mean = f @ p
        value = 0.5 * z @ z + y @ z + alpha * lognorm
        gradient = z + y - mean
        return value, gradient, p, mean

    def objective(z):
        value, gradient, _, _ = terms(z)
        return value, gradient

    def hessian(z):
        _, _, p, mean = terms(z)
        centered = f - mean[:, None]
        return np.eye(len(y)) + (centered * p) @ centered.T / alpha

    result = minimize(objective, np.zeros(len(y)), jac=True, hess=hessian,
                      method='trust-exact', options={'gtol': 1e-7, 'maxiter': maxiter})
    value, gradient, p, _ = terms(result.x)
    residual = f @ p - y
    logits = logq - f.T @ result.x / alpha
    logp = logits - logsumexp(logits)
    kl = float(p @ (logp - logq))
    primal = float(0.5 * residual @ residual + alpha * kl)
    # Exact duality identity: Q(p(z))-D(z) = .5*||gradient(z)||^2.
    gap = float(0.5 * gradient @ gradient)
    info = {'alpha': float(alpha), 'chi2': float(residual @ residual),
            'kl': kl, 'objective': primal, 'dual_objective': float(-value),
            'duality_gap': gap, 'dual_gradient_max': float(np.max(np.abs(gradient))),
            'iterations': int(result.nit), 'solver_message': str(result.message),
            'normalization': float(p.sum())}
    if gap > 1e-8 * max(1.0, abs(primal)):
        raise RuntimeError(f'Fixed-alpha solve failed certificate: {info}')
    return p / domega, info


def solve_conv(omega, kernel, target, sigma, model, domega, steps, alpha, nr):
    for iteration in range(60):
        spectrum, info = solve_alpha(kernel, target, sigma, model, domega, alpha, steps)
        b = kernel / sigma[:, None] * np.sqrt(spectrum / domega)
        gram = b @ b.T if b.shape[0] < b.shape[1] else b.T @ b
        eigenvalues = np.maximum(eigvalsh(gram), 0)
        trace = float(np.sum(eigenvalues / (eigenvalues + alpha)))
        ratio = 2 * alpha * info['kl'] / trace
        info['alpha_ratio'] = ratio
        print(', '.join(f'{key}={value}' for key, value in info.items()), flush=True)
        if abs(ratio - 1) < 0.001:
            return spectrum
        if ratio <= 0 or not np.isfinite(ratio):
            raise RuntimeError('Invalid outer alpha ratio')
        alpha *= max(1 / ratio, 0.05)
    raise RuntimeError('Outer alpha iteration did not converge in 60 iterations')


def main(method='anneal'):
    if method not in ('anneal', 'convex'):
        raise ValueError("method must be 'anneal' or 'convex'")
    tau, g, params = read_inputs()
    ntau, idg, deltag, beta, nomega, domega, steps, alpha, nr, iflat = params
    sigma = np.maximum(np.abs(g * deltag) if idg == 0
                       else np.full(ntau, deltag), 1e-5)
    omega = np.arange(-nomega, nomega + 1) * domega
    model = np.ones(len(omega)) if iflat == 0 else np.exp(-omega**2 / (nomega * domega))
    model = model / (model.sum() * domega)
    kernel = domega * np.exp(-tau[:, None] * omega - np.logaddexp(0, -beta * omega))
    inputs = (omega, kernel, -g, sigma, model, domega, steps, alpha, nr)
    if method == 'anneal':
        spectrum = solve_anneal(*inputs)
        output_format = '%15.8f'
    elif method == 'convex':
        spectrum = solve_conv(*inputs)
        output_format = '%.15e'
    np.savetxt('dos', np.column_stack((omega, spectrum)), fmt=output_format)


if __name__ == '__main__':
    main()
