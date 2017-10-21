# From http://stackoverflow.com/a/3620972

import time
from functools import wraps

PROF_DATA = {}

def profile(fn):
    @wraps(fn)
    def with_profiling(*args, **kwargs):
        start_time = time.time()

        ret = fn(*args, **kwargs)

        elapsed_time = time.time() - start_time

        if fn.__name__ not in PROF_DATA:
            PROF_DATA[fn.__name__] = [0, []]
        PROF_DATA[fn.__name__][0] += 1
        PROF_DATA[fn.__name__][1].append(elapsed_time)

        return ret

    return with_profiling

def print_prof_data():
    for fname, data in PROF_DATA.items():
        max_time = max(data[1])
        avg_time = sum(data[1]) / len(data[1])
        tot_time = sum(data[1])
        print("Function %s called %d times. " % (fname, data[0]))
        print('Execution time total: %.3f, max: %.6f, average: %.6f' % (tot_time, max_time, avg_time))

def clear_prof_data():
    global PROF_DATA
    PROF_DATA = {}