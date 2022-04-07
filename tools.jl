include("primes.jl")

# floor(sqrt(n))
function IntegerSquareRoot(n::T) where T <: Integer
    n == 0 && return T(0)
    x = n + 1
    y = n
    while y < x
        x = y
        y = div(x + div(n, x), 2)
    end
    return x
end

# psuedo-prime
function is_prime(n::Integer)
    n==1 && return false

    for p in SmallPrimes
        n == p && return true
        n % p == 0 && return false
    end

    return powermod(2, n-1, n) == 1
end

# Return a, b such that a^2 + b^2 = n or nothing if not found.
function SumOfTwoSqures(n::Integer)
    n <= 0 && return nothing
    n == 1 && return 1, 0
    a, b = 1, 0
    for l in SmallPrimes
        e = 0
        while n % l == 0
            n = div(n, l)
            e += 1
        end
        s = l^(div(e, 2))
        a *= s
        b *= s
        if e % 2 == 1
            l % 4 == 3 && return nothing
            s, t = CornacchiaSmith(l)
            a, b = a*s - b*t, a*t + b*s
        end
    end
    if n % 4 == 1 && is_prime(n)
        s, t = CornacchiaSmith(n)
        a, b = a*s - b*t, a*t + b*s
    elseif n > 1
        return nothing
    end
    return a, b
end

# Return a, b such that a^2 + b^2 = q, where q is 2 or a prime satisfing q = 1 mod 4.
function CornacchiaSmith(q::Integer)
    q == 2 && return [1, 1]
    x = SquareRootMod(-1, q)
    a = q
    b = x
    c = IntegerSquareRoot(q)
    while b > c
        a, b = b, a % b
    end
    return b, IntegerSquareRoot(q - b^2)
end

# LegendreSymbol (a|p)
function LegendreSymbol(a::Integer, p::T) where T <: Integer
    r = powermod(a, T((p-1)//2), p)
    (r - 1) % p == 0 ? (return 1) : return -1
end

# square root of a mod p. p is a prime.
function SquareRootMod(a::Integer, p::T) where T <: Integer
    p == 2 && return 1
    p % 4 == 3 && return powermod(a, T((p+1)//4), p)
    if p % 8 == 5
        x = powermod(a, T((p+3)//8), p)
        (a - x^2) % p == 0 ? (return x) : return powermod(2, T((p-1)//4), p)*x % p
    elseif p % 8 == 1
        d = 2
        while (LegendreSymbol(d, p) + 1) % p != 0
            d += 1
        end
        e = 0
        t = p - 1
        while t % 2 == 0
            e += 1
            t >>= 1
        end
        A = powermod(a, t, p)
        D = powermod(d, t, p)
        m = 0
        for i in 0:e-1
            (powermod(A * powermod(D, m, p), 2^(e-1-i), p) + 1) % p == 0 && (m += 2^i)
        end
        x = powermod(a, T((t+1)//2), p) * powermod(D, T(m//2), p)
        return x % p
    end
    error("p is not prime")
end

# Gauss elimination modulo a prime N by row transformations
function GaussEliminationMod(M::Matrix{T}, N::Integer) where T <: Integer
    M = mod.(M, N)
    m, n = size(M)
    r = 1
    for i in 1:n
        k = 0
        for j in r:m
            M[j,i] != 0 && (k = j)
        end
        if k != 0
            M[k, :], M[r, :] = M[r, :], M[k, :]
            M[r, :] = mod.(M[r, :] * invmod(M[r, i], N), N)
            for j in 1:m
                j != i && (M[j, :] = mod.(M[j, :] - M[j, i]*M[r, :], N))
            end
            r += 1
        end
    end
    return M
end

# Chinease remainder theorem
function CRT(ms::Vector{T1}, xs::Vector{T2}) where T1 <: Integer where T2 <: Integer
    m = ms[1]
    x = xs[1]
    for i in 2:length(ms)
        _, u, v = gcdx(m, ms[i])
        x = u*m*xs[i] + v*ms[i]*x
        m *= ms[i]
        x %= m
    end
    return x
end

# return x s.t. x^2 = a mod 2^e
function SquareRootModPowerOfTwo(a::T, e::Int) where T <: Integer
    e <= 0 && error("exponent must be positive")
    a % T(2)^e == 0 && return 0
    d = 0
    while a % 2 == 0
        a = div(a, 2)
        d += 1
        e -= 1
    end
    d % 2 == 1 && return nothing
    a % T(2)^e == 1 && return T(2)^div(d, 2)
    e == 2 && (a % 4 == 3) && return nothing

    k = 3
    t = a % 8
    if t == 1
        x = T(1)
    elseif t == 4
        x = T(2)
    else
        return nothing
    end
    while k < e
        x += (x^2 - a) % T(2)^(k + 1) == 0 ? T(0) : T(2)^(k-1)
        k += 1
    end
    return x * T(2)^div(d, 2)
end

function IsPowerOfTwo(n::Integer)
    while n % 2 == 0
        n = div(n, 2)
    end
    return n == 1
end
