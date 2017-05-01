from tco import *
from memory_profiler import *


def factorial_rec(n):
    return n * factorial_rec(n-1) if n > 1 else 1


@with_continuations()
def factorial_tail_call_rec(n, k, self=None):
    return self(n-1, k*n) if n > 1 else k


def factorial_tail_call(n, acc):
    while True:
        if n == 0:
            return acc
        n, acc = n - 1, acc * n

n = 20000
sys.setrecursionlimit(max(2 * n, 100))


@profile
def rec_profiler():
    return factorial_rec(n)


@profile
def tail_call_rec_profiler():
    return factorial_tail_call_rec(n, 1)


@profile
def tail_call_profiler():
    return factorial_tail_call(n, 1)


print("{0}!={1} (using recursion)".format(n, rec_profiler()))
print("{0}!={1} (using tail-call recursion)".format(n, tail_call_rec_profiler()))
print("{0}!={1} (using tail-call)".format(n, tail_call_profiler()))