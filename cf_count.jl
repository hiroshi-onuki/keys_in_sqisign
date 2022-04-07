using Random
include("tools.jl")

function SmoothTimesPrimes(N::BigInt)
    for l in SmallPrimes
        while N % l == 0
            N = div(N, l)
        end
    end
    return is_prime(N) || N == 1
end

function FactorizationSuccess(N::BigInt, count::Int)
    ret = 0
    for _ in 1:count
        n = Random.rand(1:N)
        SmoothTimesPrimes(n) && (ret += 1)
    end
    return ret
end

function FactorizationSuccessSumOfSquare(N::BigInt, count::Int)
    ret = 0
    for _ in 1:count
        a = Random.rand(1:BigInt(floor(sqrt(N))))
        b = Random.rand(1:N - a^2)
        SmoothTimesPrimes(a^2 + b^2) && (ret += 1)
    end
    return ret
end

count = 100000
for b in [266, 322, 400, 500, 600]
    println("$b: $(FactorizationSuccess(BigInt(2)^b, count)), $(FactorizationSuccessSumOfSquare(BigInt(2)^b, count))")
end