def f(**kwargs):

    return g(**kwargs)

def g(n, m):

    print(n)
    print(m)

    return [n,m]

test = {"n":4, "m":5}
print(f(n=4, m=5))
print(f(**test))

