function get_Blackwell_TR_eqs()
    vars = @variables CaM0_R, CaM2_R, CaM4_R, CaM0_T, CaM2_T, CaM4_T

    Ca, CaM_T = @variables Ca, CaM_T
    kon_1_T, koff_1_T, kon_2_T, koff_2_T = @variables kon_1_T, koff_1_T, kon_2_T, koff_2_T
    kon_1_R, koff_1_R, kon_2_R, koff_2_R = @variables kon_1_R, koff_1_R, kon_2_R, koff_2_R
    k_TR_0, k_RT_0, k_TR_2, k_RT_2, k_TR_4, k_RT_4 =
        @variables k_TR_0, k_RT_0, k_TR_2, k_RT_2, k_TR_4, k_RT_4

    eqs = [
        0 ~
            -kon_1_R * CaM0_R * Ca^2 +
            koff_1_R * CaM2_R +
            -k_RT_0 * CaM0_R +
            k_TR_0 * CaM0_T,
        0 ~
            kon_1_R * CaM0_R * Ca^2 - koff_1_R * CaM2_R +
            -kon_2_R * CaM2_R * Ca^2 +
            koff_2_R * CaM4_R +
            -k_RT_2 * CaM2_R +
            k_TR_2 * CaM2_T,
        0 ~
            kon_2_R * CaM2_R * Ca^2 - koff_2_R * CaM4_R +
            -k_RT_4 * CaM4_R +
            k_TR_4 * CaM4_T,
        0 ~
            -kon_1_T * CaM0_T * Ca^2 + koff_1_T * CaM2_T + k_RT_0 * CaM0_R -
            k_TR_0 * CaM0_T,
        0 ~
            kon_1_T * CaM0_T * Ca^2 - koff_1_T * CaM2_T +
            -kon_2_T * CaM2_T * Ca^2 +
            koff_2_T * CaM4_T +
            k_RT_2 * CaM2_R - k_TR_2 * CaM2_T,
        0 ~
            kon_2_T * CaM2_T * Ca^2 - koff_2_T * CaM4_T + k_RT_4 * CaM4_R - k_TR_4 * CaM4_T,
        CaM_T ~ CaM0_R + CaM2_R + CaM4_R + CaM0_T + CaM2_T + CaM4_T,
    ]

    return eqs, vars
end


function solve_Blackwell_TR_eqs(
    kon_1_T,
    koff_1_T,
    kon_2_T,
    koff_2_T,
    kon_1_R,
    koff_1_R,
    kon_2_R,
    koff_2_R,
    k_TR_0,
    k_RT_0,
    k_TR_2,
    k_RT_2,
    k_TR_4,
    k_RT_4,
    Ca,
    CaM_T,
    CaM_equil,
)

    if CaM_equil
        A = [
            -k_RT_0-kon_1_R*(Ca^2) koff_1_R 0 k_TR_0 0 0
            kon_1_R*(Ca^2) -k_RT_2-koff_1_R-kon_2_R*(Ca^2) koff_2_R 0 k_TR_2 0
            0 kon_2_R*(Ca^2) -k_RT_4-koff_2_R 0 0 k_TR_4
            k_RT_0 0 0 -k_TR_0-kon_1_T*(Ca^2) koff_1_T 0
            0 k_RT_2 0 kon_1_T*(Ca^2) -k_TR_2-koff_1_T-kon_2_T*(Ca^2) koff_2_T
            ###   0 0 k_RT_4 0 kon_2_T*(Ca^2) -k_TR_4 - koff_2_T; 
            1 1 1 1 1 1
        ]
        b = [zeros(5); CaM_T]

        sol = A \ b

        return sol

    else
        return [CaM_T; zeros(5)]
    end
end


function get_Blackwell_CN_eqs()
    vars = @variables CaM0, CaM2C, CaM2N, CaM4

    Ca, CaM_T, kon_2C, koff_2C, kon_2N, koff_2N =
        @variables Ca, CaM_T, kon_2C, koff_2C, kon_2N, koff_2N

    eqs = [
        0 ~
            -kon_2C * CaM0 * Ca^2 +
            koff_2C * CaM2C +
            -kon_2N * CaM0 * Ca^2 +
            koff_2N * CaM2N,
        0 ~
            kon_2C * CaM0 * Ca^2 - koff_2C * CaM2C +
            -kon_2N * CaM2C * Ca^2 +
            koff_2N * CaM4,
        0 ~
            kon_2N * CaM0 * Ca^2 - koff_2N * CaM2N +
            -kon_2C * CaM2N * Ca^2 +
            koff_2C * CaM4,
        0 ~
            kon_2N * CaM2C * Ca^2 - koff_2N * CaM4 + kon_2C * CaM2N * Ca^2 - koff_2C * CaM4,
        CaM_T ~ CaM0 + CaM2C + CaM2N + CaM4,
    ]

    return eqs, vars
end


function solve_Blackwell_CN_eqs(kon_2C, koff_2C, kon_2N, koff_2N, Ca, CaM_t, CaM_equil)

    if CaM_equil

        CaM0 =
            (-CaM_t * koff_2C * koff_2N * kon_2C - CaM_t * koff_2C * koff_2N * kon_2N) / (
                (kon_2C + kon_2N) * (
                    -koff_2C * koff_2N - koff_2C * kon_2N * (Ca^2) -
                    koff_2N * kon_2C * (Ca^2) - kon_2C * kon_2N * (Ca^4)
                )
            )
        CaM2C =
            (-CaM_t * koff_2N * kon_2C * (Ca^2)) / (
                -koff_2C * koff_2N - koff_2C * kon_2N * (Ca^2) - koff_2N * kon_2C * (Ca^2) -
                kon_2C * kon_2N * (Ca^4)
            )
        CaM2N =
            (-CaM_t * koff_2C * kon_2N * (Ca^2)) / (
                -koff_2C * koff_2N - koff_2C * kon_2N * (Ca^2) - koff_2N * kon_2C * (Ca^2) -
                kon_2C * kon_2N * (Ca^4)
            )
        CaM4 =
            (-CaM_t * kon_2C * kon_2N * (Ca^4)) / (
                -koff_2C * koff_2N - koff_2C * kon_2N * (Ca^2) - koff_2N * kon_2C * (Ca^2) -
                kon_2C * kon_2N * (Ca^4)
            )

        return [CaM0; CaM2C; CaM2N; CaM4]
    else
        return [CaM_T; zeros(3)]
    end
end


function get_Pepke_eqs()

    x = @variables CaM0, CaM1C, CaM1N, CaM2C, CaM1N1C, CaM2N, CaM1N2C, CaM2N1C, CaM4
    CaM_T = (@variables CaM_T)[1]

    Ca, k_off_C1, k_on_C1, k_off_N1, k_on_N1, k_off_C2, k_on_C2, k_off_N2, k_on_N2 =
        @variables Ca,
        k_off_C1,
        k_on_C1,
        k_off_N1,
        k_on_N1,
        k_off_C2,
        k_on_C2,
        k_off_N2,
        k_on_N2

    A = Symbolics.variables(:a, 1:10, 1:9)
    A[1, 1] = -k_on_C1 * Ca - k_on_N1 * Ca
    A[1, 2] = k_off_C1
    A[1, 3] = k_off_N1
    A[1, 4] = 0
    A[1, 5] = 0
    A[1, 6] = 0
    A[1, 7] = 0
    A[1, 8] = 0
    A[1, 9] = 0


    A[2, 1] = k_on_C1 * Ca
    A[2, 2] = -k_off_C1 - k_on_N1 * Ca - k_on_C2 * Ca
    A[2, 3] = 0
    A[2, 4] = k_off_C2
    A[2, 5] = k_off_N1
    A[2, 6] = 0
    A[2, 7] = 0
    A[2, 8] = 0
    A[2, 9] = 0


    A[3, 1] = k_on_N1 * Ca
    A[3, 2] = 0
    A[3, 3] = -k_off_N1 - k_on_C1 * Ca - k_on_N2 * Ca
    A[3, 4] = 0
    A[3, 5] = k_off_C1
    A[3, 6] = k_off_N2
    A[3, 7] = 0
    A[3, 8] = 0
    A[3, 9] = 0


    A[4, 1] = 0
    A[4, 2] = k_on_C2 * Ca
    A[4, 3] = 0
    A[4, 4] = -k_off_C2 - k_on_N1 * Ca
    A[4, 5] = 0
    A[4, 6] = 0
    A[4, 7] = k_off_N1
    A[4, 8] = 0
    A[4, 9] = 0


    A[5, 1] = 0
    A[5, 2] = k_on_N1 * Ca
    A[5, 3] = k_on_C1 * Ca
    A[5, 4] = 0
    A[5, 5] = -k_off_N1 - k_off_C1 - k_on_C2 * Ca - k_on_N2 * Ca
    A[5, 6] = 0
    A[5, 7] = k_off_C2
    A[5, 8] = k_off_N2
    A[5, 9] = 0


    A[6, 1] = 0
    A[6, 2] = 0
    A[6, 3] = k_on_N2 * Ca
    A[6, 4] = 0
    A[6, 5] = 0
    A[6, 6] = -k_off_N2 - k_on_C1 * Ca
    A[6, 7] = 0
    A[6, 8] = k_off_C1
    A[6, 9] = 0


    A[7, 1] = 0
    A[7, 2] = 0
    A[7, 3] = 0
    A[7, 4] = k_on_N1 * Ca
    A[7, 5] = k_on_C2 * Ca
    A[7, 6] = 0
    A[7, 7] = -k_off_N1 - k_off_C2 - k_on_N2 * Ca
    A[7, 8] = 0
    A[7, 9] = k_off_N2


    A[8, 1] = 0
    A[8, 2] = 0
    A[8, 3] = 0
    A[8, 4] = 0
    A[8, 5] = k_on_N2 * Ca
    A[8, 6] = k_on_C1 * Ca
    A[8, 7] = 0
    A[8, 8] = -k_off_C1 - k_off_N2 - k_on_C2 * Ca
    A[8, 9] = k_off_C2


    A[9, 1] = 0
    A[9, 2] = 0
    A[9, 3] = 0
    A[9, 4] = 0
    A[9, 5] = 0
    A[9, 6] = 0
    A[9, 7] = k_on_N2 * Ca
    A[9, 8] = k_on_C2 * Ca
    A[9, 9] = -k_off_N2 - k_off_C2


    A[10, 1] = 1
    A[10, 2] = 1
    A[10, 3] = 1
    A[10, 4] = 1
    A[10, 5] = 1
    A[10, 6] = 1
    A[10, 7] = 1
    A[10, 8] = 1
    A[10, 9] = 1

    b = Symbolics.variables(:b, 1:10)
    b[1] = 0
    b[2] = 0
    b[3] = 0
    b[4] = 0
    b[5] = 0
    b[6] = 0
    b[7] = 0
    b[8] = 0
    b[9] = 0
    b[10] = CaM_T

    eqs = [
        (A*x)[1] ~ b[1],
        (A*x)[2] ~ b[2],
        (A*x)[3] ~ b[3],
        (A*x)[4] ~ b[4],
        (A*x)[5] ~ b[5],
        (A*x)[6] ~ b[6],
        (A*x)[7] ~ b[7],
        (A*x)[8] ~ b[8],
        (A*x)[9] ~ b[9],
        (A*x)[10] ~ b[10],
    ]

    return A, x, b, eqs
end


function get_Byrne_lobe_eqs()

    vars = @variables N0, N1, N2, N3
    CaM_T = (@variables CaM_t)[1]

    Ca, k_01, k_10, k_02, k_20, k_13, k_31, k_23, k_32 =
        @variables Ca, k01, k10, k02, k20, k13, k31, k23, k32

    eqs = [
        0 ~ -k01 * Ca * N0 + k10 * N1 + -k02 * Ca * N0 + k20 * N2
        0 ~ k01 * Ca * N0 - k10 * N1 + -k13 * Ca * N1 + k31 * N3
        0 ~ k02 * Ca * N0 - k20 * N2 + -k23 * Ca * N2 + k32 * N3
        0 ~ k13 * Ca * N1 - k31 * N3 + k23 * Ca * N2 - k32 * N3
        CaM_T ~ N0 + N1 + N2 + N3
    ]

    return eqs, vars
end


function get_Pepke_m1_lobe_eqs(Ca_free, K1, K2, CaM_t, CaM_equil)

    if CaM_equil

        X1 = CaM_t / (K1 / Ca_free + 1 + Ca_free / K2)
        X0 = K1 / Ca_free * X1
        X2 = Ca_free / K2 * X1
        return [X0; X1; X2]
    else
        return [CaM_t; 0.0; 0.0]
    end
end


function get_Faas_lobe_eqs(Ca_free, Kd_T, Kd_R, CaM_t, CaM_equil)

    if CaM_equil
        c1 = (2 * Ca_free) / Kd_T
        c2 = (2 * Kd_R) / Ca_free
        tt₀ = -(c2 * CaM_t) / (1 - (c2 + 1) * (c1 + 1))
        tr₀ = c1 * tt₀
        rr₀ = tr₀ / c2
        return [tt₀; tr₀; rr₀]
    else
        return [CaM_t; 0.0; 0.0]
    end
end


function get_Pepke_TR_eqs()

    x = @variables T0,
    T1C,
    T1N,
    T2C,
    T1N1C,
    T2N,
    T1N2C,
    T2N1C,
    T4,
    R0,
    R1C,
    R1N,
    R2C,
    R1N1C,
    R2N,
    R1N2C,
    R2N1C,
    R4
    Ca, CaM_T = @variables Ca, CaM_T

    k_off_C1_T,
    k_on_C1_T,
    k_off_N1_T,
    k_on_N1_T,
    k_off_C2_T,
    k_on_C2_T,
    k_off_N2_T,
    k_on_N2_T = @variables k_off_C1_T,
    k_on_C1_T,
    k_off_N1_T,
    k_on_N1_T,
    k_off_C2_T,
    k_on_C2_T,
    k_off_N2_T,
    k_on_N2_T

    k_off_C1_R,
    k_on_C1_R,
    k_off_N1_R,
    k_on_N1_R,
    k_off_C2_R,
    k_on_C2_R,
    k_off_N2_R,
    k_on_N2_R = @variables k_off_C1_R,
    k_on_C1_R,
    k_off_N1_R,
    k_on_N1_R,
    k_off_C2_R,
    k_on_C2_R,
    k_off_N2_R,
    k_on_N2_R

    k0_TR, k0_RT, k1_TR, k1_RT, k2_TR, k2_RT, k3_TR, k3_RT, k4_TR, k4_RT =
        @variables k0_TR, k0_RT, k1_TR, k1_RT, k2_TR, k2_RT, k3_TR, k3_RT, k4_TR, k4_RT

    eqs = [
        0 ~
            -k_on_C1_T * T0 * Ca +
            k_off_C1_T * T1C +    ### start of T states
            -k_on_N1_T * T0 * Ca +
            k_off_N1_T * T1N +
            -k0_TR * T0 +
            k0_RT * R0,
        0 ~
            k_on_C1_T * T0 * Ca - k_off_C1_T * T1C +
            -k_on_C2_T * T1C * Ca +
            k_off_C2_T * T2C +
            -k_on_N1_T * T1C * Ca +
            k_off_N1_T * T1N1C +
            -k1_TR * T1C +
            k1_RT * R1C,
        0 ~
            k_on_N1_T * T0 * Ca - k_off_N1_T * T1N +
            -k_on_N2_T * T1N * Ca +
            k_off_N2_T * T2N +
            -k_on_C1_T * T1N * Ca +
            k_off_C1_T * T1N1C +
            -k1_TR * T1N +
            k1_RT * R1N,
        0 ~
            k_on_C2_T * T1C * Ca - k_off_C2_T * T2C +
            -k_on_N1_T * T2C * Ca +
            k_off_N1_T * T1N2C +
            -k2_TR * T2C +
            k2_RT * R2C,
        0 ~
            k_on_N1_T * T1C * Ca - k_off_N1_T * T1N1C + k_on_C1_T * T1N * Ca -
            k_off_C1_T * T1N1C +
            -k_on_C2_T * T1N1C * Ca +
            k_off_C2_T * T1N2C +
            -k_on_N2_T * T1N1C * Ca +
            k_off_N2_T * T2N1C +
            -k2_TR * T1N1C +
            k2_RT * R1N1C,
        0 ~
            k_on_N2_T * T1N * Ca - k_off_N2_T * T2N +
            -k_on_C1_T * T2N * Ca +
            k_off_C1_T * T2N1C +
            -k2_TR * T2N +
            k2_RT * R2N,
        0 ~
            k_on_N1_T * T2C * Ca - k_off_N1_T * T1N2C + k_on_C2_T * T1N1C * Ca -
            k_off_C2_T * T1N2C +
            -k_on_N2_T * T1N2C * Ca +
            k_off_N2_T * T4 +
            -k3_TR * T1N2C +
            k3_RT * R1N2C,
        0 ~
            k_on_C1_T * T2N * Ca - k_off_C1_T * T2N1C + k_on_N2_T * T1N1C * Ca -
            k_off_N2_T * T2N1C +
            -k_on_C2_T * T2N1C * Ca +
            k_off_C2_T * T4 +
            -k3_TR * T2N1C +
            k3_RT * R2N1C,
        0 ~
            k_on_N2_T * T1N2C * Ca - k_off_N2_T * T4 + k_on_C2_T * T2N1C * Ca -
            k_off_C2_T * T4 +
            -k4_TR * T4 +
            k4_RT * R4,
        0 ~
            -k_on_C1_R * R0 * Ca +
            k_off_C1_R * R1C +            ### start of R states
            -k_on_N1_R * R0 * Ca +
            k_off_N1_R * R1N +
            k0_TR * T0 - k0_RT * R0,
        0 ~
            k_on_C1_R * R0 * Ca - k_off_C1_R * R1C +
            -k_on_C2_R * R1C * Ca +
            k_off_C2_R * R2C +
            -k_on_N1_R * R1C * Ca +
            k_off_N1_R * R1N1C +
            k1_TR * T1C - k1_RT * R1C,
        0 ~
            k_on_N1_R * R0 * Ca - k_off_N1_R * R1N +
            -k_on_N2_R * R1N * Ca +
            k_off_N2_R * R2N +
            -k_on_C1_R * R1N * Ca +
            k_off_C1_R * R1N1C +
            k1_TR * T1N - k1_RT * R1N,
        0 ~
            k_on_C2_R * R1C * Ca - k_off_C2_R * R2C +
            -k_on_N1_R * R2C * Ca +
            k_off_N1_R * R1N2C +
            k2_TR * T2C - k2_RT * R2C,
        0 ~
            k_on_N1_R * R1C * Ca - k_off_N1_R * R1N1C + k_on_C1_R * R1N * Ca -
            k_off_C1_R * R1N1C +
            -k_on_C2_R * R1N1C * Ca +
            k_off_C2_R * R1N2C +
            -k_on_N2_R * R1N1C * Ca +
            k_off_N2_R * R2N1C +
            k2_TR * T1N1C - k2_RT * R1N1C,
        0 ~
            k_on_N2_R * R1N * Ca - k_off_N2_R * R2N +
            -k_on_C1_R * R2N * Ca +
            k_off_C1_R * R2N1C +
            k2_TR * T2N - k2_RT * R2N,
        0 ~
            k_on_N1_R * R2C * Ca - k_off_N1_R * R1N2C + k_on_C2_R * R1N1C * Ca -
            k_off_C2_R * R1N2C +
            -k_on_N2_R * R1N2C * Ca +
            k_off_N2_R * R4 +
            k3_TR * T1N2C - k3_RT * R1N2C,
        0 ~
            k_on_C1_R * R2N * Ca - k_off_C1_R * R2N1C + k_on_N2_R * R1N1C * Ca -
            k_off_N2_R * R2N1C +
            -k_on_C2_R * R2N1C * Ca +
            k_off_C2_R * R4 +
            k3_TR * T2N1C - k3_RT * R2N1C,
        0 ~
            k_on_N2_R * R1N2C * Ca - k_off_N2_R * R4 + k_on_C2_R * R2N1C * Ca -
            k_off_C2_R * R4 + k4_TR * T4 - k4_RT * R4,
        CaM_T ~
            T0 +
            T1C +
            T1N +
            T2C +
            T1N1C +
            T2N +
            T1N2C +
            T2N1C +
            T4 +
            R0 +
            R1C +
            R1N +
            R2C +
            R1N1C +
            R2N +
            R1N2C +
            R2N1C +
            R4,
    ]

    return eqs, x
end


function solve_Pepke_TR_eq(
    k_on_C1_T,
    k_off_C1_T,
    k_on_C2_T,
    k_off_C2_T,
    k_on_N1_T,
    k_off_N1_T,
    k_on_N2_T,
    k_off_N2_T,
    k_on_C1_R,
    k_off_C1_R,
    k_on_C2_R,
    k_off_C2_R,
    k_on_N1_R,
    k_off_N1_R,
    k_on_N2_R,
    k_off_N2_R,
    Kd_0_TR,
    k0_TR,
    k0_RT,
    Kd_1_TR,
    k1_TR,
    k1_RT,
    Kd_2_TR,
    k2_TR,
    k2_RT,
    Kd_3_TR,
    k3_TR,
    k3_RT,
    Kd_4_TR,
    k4_TR,
    k4_RT,
    Ca,
    CaM_t,
)
    A = [
        -k0_TR-Ca*k_on_C1_T-Ca*k_on_N1_T k_off_C1_T k_off_N1_T 0 0 0 0 0 0 k0_RT 0 0 0 0 0 0 0 0
        Ca*k_on_C1_T -k1_TR-k_off_C1_T-Ca*k_on_C2_T-Ca*k_on_N1_T 0 k_off_C2_T k_off_N1_T 0 0 0 0 0 k1_RT 0 0 0 0 0 0 0
        Ca*k_on_N1_T 0 -k1_TR-k_off_N1_T-Ca*k_on_C1_T-Ca*k_on_N2_T 0 k_off_C1_T k_off_N2_T 0 0 0 0 0 k1_RT 0 0 0 0 0 0
        0 Ca*k_on_C2_T 0 -k2_TR-k_off_C2_T-Ca*k_on_N1_T 0 0 k_off_N1_T 0 0 0 0 0 k2_RT 0 0 0 0 0
        0 Ca*k_on_N1_T Ca*k_on_C1_T 0 -k2_TR-k_off_C1_T-k_off_N1_T-Ca*k_on_C2_T-Ca*k_on_N2_T 0 k_off_C2_T k_off_N2_T 0 0 0 0 0 k2_RT 0 0 0 0
        0 0 Ca*k_on_N2_T 0 0 -k2_TR-k_off_N2_T-Ca*k_on_C1_T 0 k_off_C1_T 0 0 0 0 0 0 k2_RT 0 0 0
        0 0 0 Ca*k_on_N1_T Ca*k_on_C2_T 0 -k3_TR-k_off_C2_T-k_off_N1_T-Ca*k_on_N2_T 0 k_off_N2_T 0 0 0 0 0 0 k3_RT 0 0
        0 0 0 0 Ca*k_on_N2_T Ca*k_on_C1_T 0 -k3_TR-k_off_C1_T-k_off_N2_T-Ca*k_on_C2_T k_off_C2_T 0 0 0 0 0 0 0 k3_RT 0
        0 0 0 0 0 0 Ca*k_on_N2_T Ca*k_on_C2_T -k4_TR-k_off_C2_T-k_off_N2_T 0 0 0 0 0 0 0 0 k4_RT
        k0_TR 0 0 0 0 0 0 0 0 -k0_RT-Ca*k_on_C1_R-Ca*k_on_N1_R k_off_C1_R k_off_N1_R 0 0 0 0 0 0
        0 k1_TR 0 0 0 0 0 0 0 Ca*k_on_C1_R -k1_RT-k_off_C1_R-Ca*k_on_C2_R-Ca*k_on_N1_R 0 k_off_C2_R k_off_N1_R 0 0 0 0
        0 0 k1_TR 0 0 0 0 0 0 Ca*k_on_N1_R 0 -k1_RT-k_off_N1_R-Ca*k_on_C1_R-Ca*k_on_N2_R 0 k_off_C1_R k_off_N2_R 0 0 0
        0 0 0 k2_TR 0 0 0 0 0 0 Ca*k_on_C2_R 0 -k2_RT-k_off_C2_R-Ca*k_on_N1_R 0 0 k_off_N1_R 0 0
        0 0 0 0 k2_TR 0 0 0 0 0 Ca*k_on_N1_R Ca*k_on_C1_R 0 -k2_RT-k_off_C1_R-k_off_N1_R-Ca*k_on_C2_R-Ca*k_on_N2_R 0 k_off_C2_R k_off_N2_R 0
        0 0 0 0 0 k2_TR 0 0 0 0 0 Ca*k_on_N2_R 0 0 -k2_RT-k_off_N2_R-Ca*k_on_C1_R 0 k_off_C1_R 0
        0 0 0 0 0 0 k3_TR 0 0 0 0 0 Ca*k_on_N1_R Ca*k_on_C2_R 0 -k3_RT-k_off_C2_R-k_off_N1_R-Ca*k_on_N2_R 0 k_off_N2_R
        0 0 0 0 0 0 0 k3_TR 0 0 0 0 0 Ca*k_on_N2_R Ca*k_on_C1_R 0 -k3_RT-k_off_C1_R-k_off_N2_R-Ca*k_on_C2_R k_off_C2_R
        ###          0 0 0 0 0 0 0 0 k4_TR 0 0 0 0 0 0 Ca*k_on_N2_R Ca*k_on_C2_R -k4_RT - k_off_C2_R - k_off_N2_R;
        1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
    ]




    b = [zeros(17); CaM_t]
    sol = A \ b
    return sol
end



function solve_Byrne_eq(k01, k10, k02, k20, k13, k31, k23, k32, Ca, CaM_T, CaM_equil)

    if CaM_equil

        A = [
            -Ca*k01-Ca*k02 k10 k20 0
            Ca*k01 -k10-Ca*k13 0 k31
            Ca*k02 0 -k20-Ca*k23 k32
            ###0 Ca*k13 Ca*k23 -k31 - k32; 
            1 1 1 1
        ]

        b = [zeros(3); CaM_T]

        sol = A \ b

        return sol

    else

        return [CaM_T; zeros(3)]
    end
end


function solve_Bwell_eq(Kd_1, Kd_2, Ca_free, CaM_t, CaM_equil)
    c1 = 1 + Kd_2 / Ca_free^2
    c2 = 1 + Ca_free^2 / Kd_1
    CaM₀ = (c1 - 1) / (c2 * c1 - 1) * CaM_t
    CaM2Ca₀ = CaM₀ * Ca_free^2 / Kd_1
    CaM4Ca₀ = (c2 - 1) / (c2 * c1 - 1) * CaM_t
    return [CaM₀; CaM2Ca₀; CaM4Ca₀]
end


function get_bhalla_equilibrium_eqs()

    CaM0, CaM2Ca, CaM3Ca, CaM4Ca, CaM_T, Ca =
        @variables CaM0, CaM2Ca, CaM3Ca, CaM4Ca, CaM_T, Ca
    kon_1, koff_1, kon_2, koff_2, kon_3, koff_3 =
        @variables kon_1, koff_1, kon_2, koff_2, kon_3, koff_3

    eq1 = -CaM0 * Ca^2 * kon_1 + CaM2Ca * koff_1 ~ 0
    eq2 = CaM0 * Ca^2 * kon_1 - CaM2Ca * koff_1 - CaM2Ca * Ca * kon_2 + CaM3Ca * koff_2 ~ 0
    eq3 = CaM2Ca * Ca * kon_2 - CaM3Ca * koff_2 - CaM3Ca * Ca * kon_3 + CaM4Ca * koff_3 ~ 0
    eq4 = CaM3Ca * Ca * kon_3 - CaM4Ca * koff_3 ~ 0
    eq5 = CaM0 + CaM2Ca + CaM3Ca + CaM4Ca ~ CaM_T

    sol = Symbolics.solve_for([eq1; eq2; eq3; eq5], [CaM0; CaM2Ca; CaM3Ca; CaM4Ca])

    return sol
end


function solve_Bhalla_eq(kon_1, koff_1, kon_2, koff_2, kon_3, koff_3, Ca, CaM_T)

    CaM0 =
        (CaM_T * koff_1 * koff_2 * koff_3) / (
            koff_1 * koff_2 * koff_3 +
            koff_2 * koff_3 * kon_1 * (Ca^2) +
            koff_3 * kon_1 * kon_2 * (Ca^3) +
            kon_1 * kon_2 * kon_3 * (Ca^4)
        )
    CaM2Ca =
        (CaM_T * koff_2 * koff_3 * kon_1 * (Ca^2)) / (
            koff_1 * koff_2 * koff_3 +
            koff_2 * koff_3 * kon_1 * (Ca^2) +
            koff_3 * kon_1 * kon_2 * (Ca^3) +
            kon_1 * kon_2 * kon_3 * (Ca^4)
        )
    CaM3Ca =
        (CaM_T * koff_3 * kon_1 * kon_2 * (Ca^3)) / (
            koff_1 * koff_2 * koff_3 +
            koff_2 * koff_3 * kon_1 * (Ca^2) +
            koff_3 * kon_1 * kon_2 * (Ca^3) +
            kon_1 * kon_2 * kon_3 * (Ca^4)
        )
    CaM4Ca =
        (CaM_T * kon_1 * kon_2 * kon_3 * (Ca^4)) / (
            koff_1 * koff_2 * koff_3 +
            koff_2 * koff_3 * kon_1 * (Ca^2) +
            koff_3 * kon_1 * kon_2 * (Ca^3) +
            kon_1 * kon_2 * kon_3 * (Ca^4)
        )

    return [CaM0; CaM2Ca; CaM3Ca; CaM4Ca]
end


function get_shifman_equilibrium_eqs()

    CaM0, CaM1Ca, CaM2Ca, CaM3Ca, CaM4Ca, CaM_T, Ca =
        @variables CaM0, CaM1Ca, CaM2Ca, CaM3Ca, CaM4Ca, CaM_T, Ca
    kon_1, koff_1, kon_2, koff_2, kon_3, koff_3, kon_4, koff_4 =
        @variables kon_1, koff_1, kon_2, koff_2, kon_3, koff_3, kon_4, koff_4

    eq1 = -CaM0 * Ca * kon_1 + CaM1Ca * koff_1 ~ 0
    eq2 = CaM0 * Ca * kon_1 - CaM1Ca * koff_1 - CaM1Ca * Ca * kon_2 + CaM2Ca * koff_2 ~ 0
    eq3 = CaM1Ca * Ca * kon_2 - CaM2Ca * koff_2 - CaM2Ca * Ca * kon_3 + CaM3Ca * koff_3 ~ 0
    eq4 = CaM2Ca * Ca * kon_3 - CaM3Ca * koff_3 - CaM3Ca * Ca * kon_4 + CaM4Ca * koff_4 ~ 0
    eq5 = CaM3Ca * Ca * kon_4 - CaM4Ca * koff_4 ~ 0
    eq6 = CaM0 + CaM1Ca + CaM2Ca + CaM3Ca + CaM4Ca ~ CaM_T

    sol = Symbolics.solve_for(
        [eq1; eq2; eq3; eq4; eq6],
        [CaM0; CaM1Ca; CaM2Ca; CaM3Ca; CaM4Ca],
    )

    return sol
end

function solve_Shifman_eq(
    kon_1,
    koff_1,
    kon_2,
    koff_2,
    kon_3,
    koff_3,
    kon_4,
    koff_4,
    Ca,
    CaM_T,
)

    CaM0 =
        (-CaM_T * koff_1 * koff_2 * koff_3 * koff_4) / (
            -koff_1 * koff_2 * koff_3 * koff_4 - Ca * koff_2 * koff_3 * koff_4 * kon_1 -
            koff_3 * koff_4 * kon_1 * kon_2 * (Ca^2) -
            koff_4 * kon_1 * kon_2 * kon_3 * (Ca^3) -
            kon_1 * kon_2 * kon_3 * kon_4 * (Ca^4)
        )
    CaM1Ca =
        (-Ca * CaM_T * koff_2 * koff_3 * koff_4 * kon_1) / (
            -koff_1 * koff_2 * koff_3 * koff_4 - Ca * koff_2 * koff_3 * koff_4 * kon_1 -
            koff_3 * koff_4 * kon_1 * kon_2 * (Ca^2) -
            koff_4 * kon_1 * kon_2 * kon_3 * (Ca^3) -
            kon_1 * kon_2 * kon_3 * kon_4 * (Ca^4)
        )
    CaM2Ca =
        (-CaM_T * koff_3 * koff_4 * kon_1 * kon_2 * (Ca^2)) / (
            -koff_1 * koff_2 * koff_3 * koff_4 - Ca * koff_2 * koff_3 * koff_4 * kon_1 -
            koff_3 * koff_4 * kon_1 * kon_2 * (Ca^2) -
            koff_4 * kon_1 * kon_2 * kon_3 * (Ca^3) -
            kon_1 * kon_2 * kon_3 * kon_4 * (Ca^4)
        )
    CaM3Ca =
        (-CaM_T * koff_4 * kon_1 * kon_2 * kon_3 * (Ca^3)) / (
            -koff_1 * koff_2 * koff_3 * koff_4 - Ca * koff_2 * koff_3 * koff_4 * kon_1 -
            koff_3 * koff_4 * kon_1 * kon_2 * (Ca^2) -
            koff_4 * kon_1 * kon_2 * kon_3 * (Ca^3) -
            kon_1 * kon_2 * kon_3 * kon_4 * (Ca^4)
        )
    CaM4Ca =
        (-CaM_T * kon_1 * kon_2 * kon_3 * kon_4 * (Ca^4)) / (
            -koff_1 * koff_2 * koff_3 * koff_4 - Ca * koff_2 * koff_3 * koff_4 * kon_1 -
            koff_3 * koff_4 * kon_1 * kon_2 * (Ca^2) -
            koff_4 * kon_1 * kon_2 * kon_3 * (Ca^3) -
            kon_1 * kon_2 * kon_3 * kon_4 * (Ca^4)
        )

    return [CaM0; CaM1Ca; CaM2Ca; CaM3Ca; CaM4Ca]
end
