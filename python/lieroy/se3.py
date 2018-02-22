
import numpy as np
import lieroy.lieroy_core as core

def log(T):
    return core.se3_log(T)

def exp(xsi):
    return core.se3_exp(T)

def sample(mean=np.identity(4), covariance=np.identity(6), n_samples=100):
    return core.se3_sample_normal_distribution(mean, covariance, n_samples)
