using HSL
using LinearOperators

K = diagm(1 ./ [1.0; 2.0; 3.0; 4.0; 5.0])
rhs = ones(5)
x1=K\rhs
spK = sparse(K)

K[3,3]=eps()

LBL = Ma97(K)
ma97_factorize(LBL)
x2 = ma97_solve(LBL, rhs)  # or x = LBL \ b
ret = ma97_inquire(LBL)
d=ret[2]
piv=ret[1]
show(piv);println()
show(d);println()
d[1,3] = 3.0
show(d);println()

ma97_alter(LBL, d)

x3 = ma97_solve(LBL, rhs)

show(x3)

#for i=1:5
    K = sprand(15,15,0.4)
    K = K + K'
    LBL = Ma97(K)
    ma97_factorize(LBL)
    ret = ma97_inquire(LBL)
    d=ret[2]
    piv=ret[1]
    println(" pivots : $piv")
    apiv = abs(piv)
    D1=spdiagm((d[1,apiv],d[2,apiv][1:end-1],d[2,apiv][1:end-1]),[0 -1 1])
        D2=spdiagm((d[1,:],d[2,[1:end-1]],d[2,[1:end-1]]),[0 -1 1])
D=D2[apiv,apiv]
    b=ones(15)
    x2 = ma97_solve(LBL, b)
    p1=copy(b)
    ma97_solve!(LBL, p1, job=:PL)
    p2 = D*p1
    p3 = copy(p2)
    ma97_solve!(LBL, p3, job=:LPS)
                                                                      
show([x2-p3])
println("\n")

    b=ones(15)
    x2 = ma97_solve(LBL, b)
    p1=copy(b)
    ma97_solve!(LBL, p1, job=:PL)
    p2 = (D*p1)[apiv]
    p3 = copy(p2)
    ma97_solve!(LBL, p3, job=:LPS)
                                                                      
show([x2-p3])
println("\n")

    b=ones(15)
    x2 = ma97_solve(LBL, b)
    p1=copy(b)
    ma97_solve!(LBL, p1, job=:PL)
    p2[apiv] = D*p1[apiv]
    p3 = copy(p2)
    ma97_solve!(LBL, p3, job=:LPS)
                                                                      
show([x2-p3])
println("\n")
    b=ones(15)
    x2 = ma97_solve(LBL, b)
    p1=copy(b)
    ma97_solve!(LBL, p1, job=:PL)
    p2[apiv] = D*p1
    p3 = copy(p2)
    ma97_solve!(LBL, p3, job=:LPS)
                                                                      
[x2-p3]


#end
