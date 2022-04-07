import Base.+, Base.-, Base.*, Base./, Base.==
import LinearAlgebra

# Quaternion algebra over Q ramified at infinity and p
function Quaternion(p::T) where T <: Integer
    return x -> quat_elem(x, p)
end

struct quat_elem
    v::Vector{Rational{BigInt}}
    p::BigInt
    quat_elem(v::Vector{Rational{BigInt}}, p::BigInt) = length(v) == 4 ? new(v, p) : error("length of vector is not 4")
    quat_elem(v::Vector{T1}, p::T2) where T1 <: Union{Integer, Rational} where T2 <: Integer = quat_elem(Vector{Rational{BigInt}}(v), BigInt(p))
end

+(x::quat_elem, y::quat_elem) = x.p == y.p ? quat_elem(x.v + y.v, x.p) : error("different ramified primes")
-(x::quat_elem, y::quat_elem) = x.p == y.p ? quat_elem(x.v - y.v, x.p) : error("different ramified primes")
-(x::quat_elem) = quat_elem(-x.v, x.p)
/(x::quat_elem, n::BigInt) = quat_elem(x.v ./ n, x.p)
==(x::quat_elem, y::quat_elem) = x.p == y.p && x.v == y.v

function *(x::quat_elem, y::quat_elem)
    if x.p != y.p
        error("different ramified primes")
    else
        c1 = x.v[1]*y.v[1] - x.v[2]*y.v[2] - x.p*x.v[3]*y.v[3] - x.p*x.v[4]*y.v[4]
        ci = x.v[1]*y.v[2] + x.v[2]*y.v[1] + x.p*x.v[3]*y.v[4] - x.p*x.v[4]*y.v[3]
        cj = x.v[1]*y.v[3] - x.v[2]*y.v[4] + x.v[3]*y.v[1] + x.v[4]*y.v[2]
        ck = x.v[1]*y.v[4] + x.v[2]*y.v[3] - x.v[3]*y.v[2] + x.v[4]*y.v[1]
        return quat_elem([c1, ci, cj, ck], x.p)
    end
end

involution(x::quat_elem) = quat_elem([x.v[1], -x.v[2], -x.v[3], -x.v[4]], x.p)
Nrd(x::quat_elem) = x.v[1]^2 + x.v[2]^2 + x.p*x.v[3]^2 + x.p*x.v[4]^2
Trd(x::quat_elem) = 2*x.v[1]
bilinear(x::quat_elem, y::quat_elem) = Trd(involution(x) * y)

function LeftProductMatrix(x::quat_elem)
    return hcat([
        (x*quat_elem([1,0,0,0], x.p)).v,
        (x*quat_elem([0,1,0,0], x.p)).v,
        (x*quat_elem([0,0,1,0], x.p)).v,
        (x*quat_elem([0,0,0,1], x.p)).v
        ]...)
end

function RightProductMatrix(x::quat_elem)
    return hcat([
        (quat_elem([1,0,0,0], x.p)*x).v,
        (quat_elem([0,1,0,0], x.p)*x).v,
        (quat_elem([0,0,1,0], x.p)*x).v,
        (quat_elem([0,0,0,1], x.p)*x).v
        ]...)
end

function Base.show(io::IO, x::quat_elem)
    print(io, "$(x.v[1]) + $(x.v[2])i + $(x.v[3])j + $(x.v[4])k")
end

# Order
struct OrderData
    p::BigInt
    basis::Vector{quat_elem}
    OrdtoQuat::Matrix{Rational{BigInt}}
    QuattoOrd::Matrix{Rational{BigInt}}
    prod_matrices::Vector{Matrix{BigInt}}
    involution_martix::Matrix{BigInt}
    trace_vector::Vector{BigInt}
    bilinear_matrix::Matrix{BigInt}
end
global_order_data = OrderData[]

function Order(b::Vector{quat_elem})
    if length(b) != 4
        error("length of vector is not 4")
    elseif !prod([b[1].p == b[i].p for i in 2:4])
        error("different ramified primes")
    end

    # transformation matrix between order and quaternion
    OrdtoQuat = hcat([b[i].v for i in 1:4]...)
    if LinearAlgebra.det(OrdtoQuat) == 0
        error("input basis is not linearly independent")
    end
    QuattoOrd = LinearAlgebra.inv(OrdtoQuat)

    # xy = sum(x[i]*prod_matrices[i])*y
    prod_matrices = [QuattoOrd*hcat([(b[i]*b[j]).v for j in 1:4]...) for i in 1:4]
    if !prod([denominator(x) == 1 for x in hcat(prod_matrices...)])
        error("not muliplicatively closed")
    end

    # bar(x) = involution_martix*x
    involution_martix = QuattoOrd * hcat([involution(x).v for x in b]...)

    # Trd(x) = (trace_vector)^t * x
    trace_vector = [BigInt(Trd(x)) for x in b]

    # Nrd(x + y) - Nrd(x) - Nrd(y) = x*bilinear_matrix*y 
    bilinear_matrix = hcat([[BigInt(bilinear(x, y)) for y in b] for x in b]...)

    O = OrderData(b[1].p, b, OrdtoQuat, QuattoOrd, prod_matrices, involution_martix, trace_vector, bilinear_matrix)
    push!(global_order_data, O)
    id = length(global_order_data)
    return x->order_elem(x, id), id
end


struct order_elem
    v::Vector{BigInt}
    id::Int
    order_elem(v::Vector{T}, id::Int) where T <: Integer = length(v) == 4 ? new(v, id) : error("length of vector is not 4")
end

+(x::order_elem, y::order_elem) = x.id == y.id ? order_elem(x.v + y.v, x.id) : error("elements in different orders")
-(x::order_elem, y::order_elem) = x.id == y.id ? order_elem(x.v - y.v, x.id) : error("elements in different orders")
-(x::order_elem) = order_elem(-x.v, x.id)
*(x::order_elem, y::order_elem) = x.id == y.id ? order_elem(sum(x.v .* global_order_data[x.id].prod_matrices)*y.v, x.id) : error("elements in different orders")
*(x::BigInt, y::order_elem) = order_elem(x*y.v, y.id)
/(x::order_elem, n::Integer) = order_elem(BigInt.(x.v//n), x.id)
==(x::order_elem, y::order_elem) = x.id == y.id && x.v == y.v

involution(x::order_elem) = order_elem(global_order_data[x.id].involution_martix * x.v, x.id)
Trd(x::order_elem) = transpose(global_order_data[x.id].trace_vector) * x.v
Nrd(x::order_elem) = div(transpose(x.v) * global_order_data[x.id].bilinear_matrix * x.v,  2)
bilinear(x::order_elem, y::order_elem) = x.id == y.id ? transpose(x.v) * global_order_data[x.id].bilinear_matrix * y.v : error("elements in different orders")
ToQuaternion(x::order_elem) = quat_elem(global_order_data[x.id].OrdtoQuat * x.v, global_order_data[x.id].basis[1].p)
ToOrderElement(x::quat_elem, id::Int) = id <= length(global_order_data) ? order_elem(BigInt.(global_order_data[id].QuattoOrd * x.v), id) : error("order id out of range")
ToOrderElement(x::Vector{T}, id::Int) where T <: Integer = id <= length(global_order_data) ? order_elem(x, id) : error("order id out of range")
bilinear_matrix(id::Int) = id <= length(global_order_data) ? global_order_data[id].bilinear_matrix : error("order id out of range")
ramified_prime(id::Int) = id <= length(global_order_data) ? global_order_data[id].basis[1].p : error("order id out of range")

# Ideal
struct LeftIdeal
    basis::Vector{order_elem}
    order_id::Int
    p::BigInt
    function LeftIdeal(b::Vector{order_elem})
        id = b[1].id
        if length(b) != 4
            error("length of vector is not 4")
        elseif !(prod([id == b[i].id for i in 2:4]))
            error("elements in different orders")
        end
        M = hcat([x.v for x in b]...)
        if LinearAlgebra.det(Rational{BigInt}.(M)) == 0
            error("input basis is not linearly independent")
        end
        trans_mat = LinearAlgebra.inv(Rational{BigInt}.(M))
        for A in global_order_data[id].prod_matrices
            if abs(prod(denominator.(trans_mat * A * M))) != 1
                error("not left ideal")
            end
        end

        return new(b, id, global_order_data[id].p)
    end
end

*(I::LeftIdeal, x::order_elem) = LeftIdeal([b*x for b in I.basis])
/(I::LeftIdeal, n::Integer) = LeftIdeal([b/n for b in I.basis])
Nrd(I::LeftIdeal) = IntegerSquareRoot(BigInt(abs(LinearAlgebra.det(Rational{BigInt}.(hcat([x.v for x in I.basis]...))))))

function isin(x::order_elem, I::LeftIdeal)
    x.id == I.basis[1].id || return false
    trans_mat = LinearAlgebra.inv(Rational{BigInt}.(hcat([b.v for b in I.basis]...)))
    return prod(denominator.(trans_mat * x.v)) == 1
end

isin(I::LeftIdeal, J::LeftIdeal) = prod([isin(b, J) for b in I.basis])
==(I::LeftIdeal, J::LeftIdeal) = isin(I, J) && isin(J, I)

# left O-ideal Ox + ON
function LeftIdeal(x::order_elem, N::BigInt)
    units = [[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1]]
    basis = GetBasis(vcat([(ToOrderElement(u, x.id)*x).v for u in units], units*N))
    return LeftIdeal(ToOrderElement.(basis, x.id))
end

# Return a quaternion a such that I = (a, Nrd(I)).
function IdealGenerator(I::LeftIdeal)
    N = Nrd(I)
    a = ToOrderElement([0,0,0,0], I.basis[1].id)
    i = 1
    while gcd(a.v) != 1 || gcd(Nrd(a), N^2) != N
        a += I.basis[i]
        i = (i % 4) + 1
    end
    return a
end

# Return a primitive element in I
function PrimitiveElement(I::LeftIdeal)
    a = ToOrderElement([0,0,0,0], I.basis[1].id)
    i = 1
    while gcd(a.v) != 1
        a += I.basis[i]
        i = (i % 4) + 1
    end
    return a
end

# left O-ideal I + ON
function GeneratedBy(I::LeftIdeal, N::Integer)
    bs = vcat([b.v for b in I.basis], Vector{BigInt}[[N, 0, 0, 0],[0,N,0,0],[0,0,N,0],[0,0,0,N]])
    bs = GetBasis(bs)
    return LeftIdeal([ToOrderElement(b, I.order_id) for b in bs])
end

# smallest element in I
function SmallestElement(I::LeftIdeal)
    # LLL reduction
    H = IntegralLLL([b.v for b in I.basis], bilinear_matrix(id))
    LLLmat = hcat([b.v for b in I.basis]...) * H
    red_basis = [LLLmat[:, i] for i in 1:4]

    Nbound = 2*Nrd(ToOrderElement(red_basis[1], id))
    vecs = ShortVectors(red_basis, bilinear_matrix(id), Nbound)
    sort!(vecs, lt=(x,y)->x[2]<y[2])
    return ToOrderElement(vecs[1][1], I.order_id)
end