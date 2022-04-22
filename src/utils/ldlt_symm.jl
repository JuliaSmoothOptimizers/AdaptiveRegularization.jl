export ldlt_symm

"""
Higams' ldlt_symm translated from matlab. Performs  a so called
 BKK  Bounded Bunch Kaufman factorization of A0, that means
 ||L|| is bounded and bounded away from zero.
"""
function ldlt_symm(A0::Array{T,2}, piv::Char = 'r') where {T}

    #LDLT_SYMM  Block LDL^T factorization for a symmetric indefinite matrix.
    #     Given a Hermitian matrix A,
    #     [L, D, P, RHO, NCOMP] = LDLT_SYMM(A, PIV) computes a permutation P,
    #     a unit lower triangular L, and a real block diagonal D
    #     with 1x1 and 2x2 diagonal blocks, such that  P*A*P' = L*D*L'.
    #     PIV controls the pivoting strategy:
    #       PIV = 'p': partial pivoting (Bunch and Kaufman),
    #       PIV = 'r': rook pivoting (Ashcraft, Grimes and Lewis).
    #     The default is partial pivoting.
    #     RHO is the growth factor.
    #     For PIV = 'r' only, NCOMP is the total number of comparisons.

    #     References:
    #     J. R. Bunch and L. Kaufman, Some stable methods for calculating
    #        inertia and solving symmetric linear systems, Math. Comp.,
    #        31(137):163-179, 1977.
    #     C. Ashcraft, R. G. Grimes and J. G. Lewis, Accurate symmetric
    #        indefinite linear equation solvers. SIAM J. Matrix Anal. Appl.,
    #        20(2):513-561, 1998.
    #     N. J. Higham, Accuracy and Stability of Numerical Algorithms,
    #        Second edition, Society for Industrial and Applied Mathematics,
    #        Philadelphia, PA, 2002; chap. 11.

    #    This routine does not exploit symmetry and is not designed to be
    #     efficient.
    #
    #    Adapted from N. Higham Matlab code.
    #    JPD april 23 2015


    # Since array contents are mutable and modified,
    # copy the input matrix to keep it unchanged outside the
    # scope of the function

    A = copy(A0)

    # minimal checks for conforming inputs
    isequal(triu(A)', tril(A)) || error("Must supply Hermitian matrix.")
    piv in ['p', 'r'] || error("Pivoting must be \"'p'\" or \" 'r' \".")

    n, = size(A)

    k = 1
    # D = eye(n,n);
    # L = eye(n,n);
    global D = Matrix{T}(1.0I, n, n)
    global L = Matrix{T}(1.0I, n, n)
    if n == 1
        D = A
    end

    global pp = collect(1:n)

    maxA = norm(A, Inf)
    global ρ = maxA

    global ncomp = 0
    s = 1

    α = T((1 + sqrt(17)) / 8)
    while k < n
        (λ, vr) = findmax(abs.(A[k+1:n, k]))
        r = vr[1] + k
        if λ > 0
            swap = false
            if abs.(A[k, k]) >= α * λ
                s = 1
            else
                if piv == 'p'
                    temp = A[k:n, r]
                    temp[r-k+1] = 0
                    σ = norm(temp, Inf)
                    if α * λ^2 <= abs.(A[k, k]) * σ
                        s = 1
                    elseif abs.(A[r, r]) >= α * σ
                        swap = true
                        m1 = k
                        m2 = r
                        s = 1
                    else
                        swap = true
                        m1 = k + 1
                        m2 = r
                        s = 2
                    end
                    if swap
                        A[[m1, m2], :] = A[[m2, m1], :]
                        L[[m1, m2], :] = L[[m2, m1], :]
                        A[:, [m1, m2]] = A[:, [m2, m1]]
                        L[:, [m1, m2]] = L[:, [m2, m1]]

                        pp[[m1, m2]] = pp[[m2, m1]]
                    end
                elseif piv == 'r'
                    j = k
                    pivot = false
                    λ_j = λ
                    while ~pivot
                        (temp1, vr) = findmax(abs.(A[k:n, j]))
                        global ncomp = ncomp + n - k
                        r = vr[1] + k - 1
                        temp = A[k:n, r]
                        temp[r-k+1] = 0.0
                        λᵣ = max(maximum(temp), -minimum(temp))
                        global ncomp = ncomp + n - k
                        if α * λᵣ <= abs.(A[r, r])
                            pivot = true
                            s = 1
                            A[k, :], A[r, :] = A[r, :], A[k, :]
                            L[k, :], L[r, :] = L[r, :], L[k, :]
                            A[:, k], A[:, r] = A[:, r], A[:, k]
                            L[:, k], L[:, r] = L[:, r], L[:, k]
                            pp[k], pp[r] = pp[r], pp[k]
                        elseif λ_j == λᵣ
                            pivot = true
                            s = 2
                            A[k, :], A[j, :] = A[j, :], A[k, :]
                            L[k, :], L[j, :] = L[j, :], L[k, :]
                            A[:, k], A[:, j] = A[:, j], A[:, k]
                            L[:, k], L[:, j] = L[:, j], L[:, k]
                            pp[k], pp[j] = pp[j], pp[k]
                            k1 = k + 1
                            A[k1, :], A[r, :] = A[r, :], A[k1, :]
                            L[k1, :], L[r, :] = L[r, :], L[k1, :]
                            A[:, k1], A[:, r] = A[:, r], A[:, k1]
                            L[:, k1], L[:, r] = L[:, r], L[:, k1]
                            pp[k1], pp[r] = pp[r], pp[k1]
                        else
                            j = r
                            s = 1
                            λ_j = λᵣ
                        end
                    end
                end
            end

            if s == 1

                D[k, k] = A[k, k]
                A[k+1:n, k] = A[k+1:n, k] / A[k, k]
                L[k+1:n, k] = A[k+1:n, k]
                i = k+1:n
                A[i, i] = A[i, i] - A[i, k:k] * A[k:k, i]
                A[i, i] = 0.5 * (A[i, i] + A[i, i]')

            elseif s == 2

                E = A[k:k+1, k:k+1]
                D[k:k+1, k:k+1] = E
                i = k+2:n
                C = A[i, k:k+1]
                temp = C / E
                L[i, k:k+1] = temp
                A[i, k+2:n] = A[i, k+2:n] - temp * C'
                A[i, i] = 0.5 * (A[i, i] + A[i, i]')
            end

            # Ensure diagonal real (see LINPACK User's Guide, p. 5.17).
            for i = k+s:n
                A[i, i] = real(A[i, i])
            end

            if k + s <= n
                val, = findmax(abs.(A[k+s:n, k+s:n]))
                global ρ = max(ρ, val)
            end

        else  # Nothing to do.
            s = 1
            D[k, k] = A[k, k]
        end

        k = k + s

        if k == n
            D[n, n] = A[n, n]
            break
        end

    end
    #P=eye(n)
    #P = P[pp,:]
    return L, D, pp, ρ, ncomp
end
