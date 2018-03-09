
import numpy as np

import lieroy.parallel
import lieroy.lieroy_core


if __name__ == '__main__':
    wrapper = lieroy.parallel.FunctionWrapper('se3_log', 'lieroy.lieroy_core')
    i = wrapper(np.identity(4))
    print(i)

    print(pk.dumps(wrapper))
