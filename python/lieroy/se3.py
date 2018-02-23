
import numpy as np
import lieroy.lieroy_core as core

def log(T):
    return core.se3_log(T)

def exp(xsi):
    return core.se3_exp(xsi)

def sample_normal(mean=np.identity(4), covariance=np.identity(6), n_samples=100):
    return core.se3_sample_normal_distribution(mean, covariance, n_samples)

def sample_uniform(mean=np.identity(4), range_translation=1.0, range_rotation=1.0, n_samples=100):
    return core.se3_sample_uniform_distribution(mean, range_translation, range_rotation, n_samples);
