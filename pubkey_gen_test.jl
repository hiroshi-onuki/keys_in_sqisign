using Nemo
using Statistics
using Random
using Printf

prime_name = "p256"

include("tools.jl")
include("quaternion.jl")
include("lattice.jl")

p = BigInt(2)^33*5^21*7^2*11*31*83*107*137*751*827*3691*4019*6983*517434778561*26602537156291 - 1

B = Quaternion(p)
basis = [
        B([1,0,0,0]),
        B([0,1,0,0]),
        B([1//2,0,1//2,0]),
        B([0,1//2,0,1//2])
    ]
O, id = Order(basis)

# [alpha]*I
function pushforward(alpha::order_elem, I::LeftIdeal)
    Oalpha = LeftIdeal(alpha, Nrd(alpha))
    J = Oalpha*ToOrderElement([Nrd(I),0,0,0], id)
    K = I*ToOrderElement([Nrd(alpha),0,0,0], id)
    bs = GetBasis([b.v for b in vcat(J.basis, K.basis)])
    L = LeftIdeal(ToOrderElement.(bs, id))
    return L * involution(alpha) / Nrd(alpha)
end

# prime ideal of norm N twisted by alpha = 1 + ni
function PrimeIdeal(I::LeftIdeal, N::BigInt, n::BigInt)
    alpha = ToOrderElement([0,1,0,0], id)
    if n != N
        alpha = ToOrderElement([1,0,0,0], id) + n*alpha
    end
    return pushforward(alpha, I)
end

function secret_key(I::LeftIdeal)
    N = 1
    Nmax = BigInt(floor(p^(1//4)))
    Nmin = BigInt(1)
    while !(N % 8 == 3) || TDataT1.degree % N == 0 || !is_prime(N) 
        N = Random.rand(Nmin:Nmax)
    end
    n = Random.rand(0:N)
    return PrimeIdeal(I, N, n)
end

function RepresentIntegerRandom(p::BigInt, N::BigInt)
    B = BigInt(floor(sqrt(N//p)))
    while true
        c, d = Random.rand(-B:B, 2)
        ret = SumOfTwoSqures(N - p*(c^2 + d^2))
        ret != nothing && return quat_elem([ret[1], ret[2], c, d], p)
    end
end

function FullRepresentIntegerRandom(p::BigInt, N::BigInt)
    B = BigInt(floor(sqrt(4N//p)))
    while true
        d = Random.rand(-B:B)
        Bd = BigInt(floor(sqrt(4N//p - d^2)))
        c = Random.rand(-Bd:Bd)
        ret = SumOfTwoSqures(4N - p*(c^2 + d^2))
        ret != nothing && (ret[1] - c) % 2 == 0 && (ret[2] - d) % 2 == 0 && return quat_elem([ret[1], ret[2], c, d], p)/BigInt(2)
    end
end

function seckey_sqisign(p::BigInt, N::BigInt, is_full::Bool, ex::Int)
    n = Int(ceil(log(2, p))) + ex
    if is_full
        gamma = ToOrderElement(FullRepresentIntegerRandom(p, BigInt(2)^n*N), id)
    else
        gamma = ToOrderElement(RepresentIntegerRandom(p, BigInt(2)^n*N), id)
    end
    return LeftIdeal(gamma, N)
end

function save_to_file(data::Vector{Int}, filename::String)
    fp = open(filename, "w")
    for d in data
        println(fp, d)
    end
    close(fp)
end

function IdealList(N::BigInt)
    n = Int(ceil(log(2, p)))
    alpha = ToOrderElement(RepresentIntegerRandom(p, BigInt(2)^n*N), id)
    I = LeftIdeal(alpha, N)

    Is = LeftIdeal[]
    for s in 0:N
        push!(Is, PrimeIdeal(I, N, s))
    end
    return Is
end

function make_data(N::BigInt, Is::Vector{LeftIdeal}, n::Int, is_full::Bool, ex::Int)
    data = Int[]
    for _ in 1:n
        I = seckey_sqisign(p, N, is_full, ex)
        idx = findfirst(x->x==I, Is)
        push!(data, idx)
    end
    return data
end

function save_to_file(data::Vector{Int}, filename::String)
    fp = open(filename, "w")
    for d in data
        println(fp, d)
    end
    close(fp)
end


# main()
n = 100
for N in BigInt[211, 223, 227, 239, 251]
    println("N=",N)

    Is = IdealList(N) 
    println("all ideals are listed.")

    ps = []
    for ex in [0, 5, 11]
        for is_full in [false, true]
            data, runtime = @timed make_data(N, Is, Int(N*n), is_full, ex)
            println("calc. time: ", runtime)

            data_n = Int[]
            for s in 1:N+1
                push!(data_n, length(filter(x->x==s, data)))
            end
            println("num of generated keys: $(length(filter(x->x>0, data_n))) ($(@sprintf("%.2f", round(length(filter(x->x>0, data_n))/(N+1)*100,digits=2)))%), ex=$(ex), $(is_full ? "FullRep." : "Rep")")
            save_to_file(data, "data/secret_$(N)_$(ex)_$(is_full ? "full" : "").txt")
        end
    end
end
