from random import Random

def factorization_success_estimate(n, b):
    r = 0
    for se in range(n+1):
        lp = n - se
        if se <= b:
            rho = 1
        else:
            rho = dickman_rho.approximate(se/b)
        if lp <= 1:
            rp = 1
        else:
            rp = float(li(2^lp) - li(2^(lp-1)))
        r += rho * rp / 2^lp
    return r

if __name__=="__main__":
    b = 20
    for n in [266, 322, 400, 500, 600]:
        print(n, factorization_success_estimate(n, b))

