
import json
import numpy as np
import lieroy.lieroy_core as core
import sys

def log(T):
    return core.se3_log(T)

def exp(xsi):
    return core.se3_exp(xsi)

def sample_normal(mean=np.identity(4), covariance=np.identity(6), n_samples=100):
    return core.se3_sample_normal_distribution(mean, covariance, n_samples)

def sample_uniform(mean=np.identity(4), range_translation=1.0, range_rotation=1.0, n_samples=100):
    return core.se3_sample_uniform_distribution(mean, range_translation, range_rotation, n_samples);

def adjoint(T):
    return core.se3_adjoint(T)

def log_cli():
    T = np.array(json.load(sys.stdin))

    xsi = log(T)

    json.dump(xsi.tolist(), sys.stdout)

def gaussian_from_sample(sample):
    return core.se3_gaussian_distribution_of_sample(sample)

def compound_poses(t1, sigma1, t2, sigma2, second_order=False):
    if second_order:
        print('Second order!!')
        adj_t1 = adjoint(t1)
        return t1 @ t2, sigma1 + adj_t1 @ (sigma2 @ adj_t1.T)
    else:
        return core.se3_compound_poses(t1, sigma1, t2, sigma2)


def compound_poses_py(t1, sigma1, t2, sigma2, second_order=False):
    """
    4th order compounding of se3 poses.
    see barfoot book p. 276, Eq. 7.204.
    """
    def french_quote(t):
        return -1.0 * np.trace(t) * np.identity(3) + t

    def double_french_quote(a, b):
        return french_quote(a) @ french_quote(b) + french_quote(b @ a)

    t1_adj = adjoint(t1)
    sigma2_prime = t1_adj @ (sigma2 @ t1_adj.T)
    second_order_cov = sigma1 + sigma2_prime

    if second_order:
        return t1 @ t2, second_order_cov

    sigma1_rho_rho = sigma1[0:3, 0:3]
    sigma1_rho_phi = sigma1[0:3, 3:6]
    sigma1_phi_phi = sigma1[3:6, 3:6]

    sigma2_rho_rho = sigma2[0:3, 0:3]
    sigma2_rho_phi = sigma2[0:3, 3:6]
    sigma2_phi_phi = sigma2[3:6, 3:6]

    a1 = np.zeros((6,6))
    a1[0:3,0:3] = french_quote(sigma1_phi_phi)
    a1[0:3,3:6] = french_quote(sigma1_rho_phi + sigma1_rho_phi.T)
    a1[3:6,3:6] = french_quote(sigma1_phi_phi)

    a2 = np.zeros((6,6))
    a2[0:3,0:3] = french_quote(sigma2_phi_phi)
    a2[0:3,3:6] = french_quote(sigma2_rho_phi + sigma2_rho_phi.T)
    a2[3:6,3:6] = french_quote(sigma2_phi_phi)

    b_rho_rho = (
        double_french_quote(sigma1_phi_phi, sigma2_rho_rho) +
        double_french_quote(sigma1_rho_phi.T, sigma2_rho_phi) +
        double_french_quote(sigma1_rho_phi, sigma2_rho_phi.T) +
        double_french_quote(sigma1_rho_rho, sigma2_phi_phi)
    )

    b_rho_phi = (
        double_french_quote(sigma1_phi_phi, sigma2_rho_phi.T) +
        double_french_quote(sigma1_rho_phi.T, sigma2_phi_phi)
    )

    b_phi_phi = double_french_quote(sigma1_phi_phi, sigma2_phi_phi)


    a_terms = a1 @ sigma2_prime + sigma2_prime @ a1.T + a2 @ sigma1 + sigma1 @ a1.transpose()
    b = np.zeros((6,6))
    b[0:3,0:3] = b_rho_rho
    b[0:3,3:6] = b_rho_phi
    b[3:6, 0:3] = b_rho_phi.T
    b[3:6,3:6] = b_phi_phi

    return t1 @ t2, second_order_cov + 0.25 * b + (1.0 / 12.0) * a_terms

