function get_Blackwell_TR_eqs()
    vars = @variables CaM0_R, CaM2_R, CaM4_R, CaM0_T, CaM2_T, CaM4_T

    Ca, CaM_T = @variables Ca, CaM_T
    kon_1_T, koff_1_T, kon_2_T, koff_2_T = @variables kon_1_T, koff_1_T, kon_2_T, koff_2_T
    kon_1_R, koff_1_R, kon_2_R, koff_2_R = @variables kon_1_R, koff_1_R, kon_2_R, koff_2_R
    k_TR_0, k_RT_0, k_TR_2, k_RT_2, k_TR_4, k_RT_4 = @variables k_TR_0, k_RT_0, k_TR_2, k_RT_2, k_TR_4, k_RT_4

    eqs = [
        0 ~ -kon_1_R * CaM0_R * Ca^2 + koff_1_R * CaM2_R + 
            -k_RT_0  * CaM0_R + k_TR_0 * CaM0_T,
        0 ~  kon_1_R * CaM0_R * Ca^2 - koff_1_R * CaM2_R +
            -kon_2_R * CaM2_R * Ca^2 + koff_2_R * CaM4_R +
            -k_RT_2  * CaM2_R + k_TR_2 * CaM2_T,
        0 ~  kon_2_R * CaM2_R * Ca^2 - koff_2_R * CaM4_R +
            -k_RT_4  * CaM4_R + k_TR_4 * CaM4_T,
        0 ~ -kon_1_T * CaM0_T * Ca^2 + koff_1_T * CaM2_T + 
             k_RT_0  * CaM0_R - k_TR_0 * CaM0_T,
        0 ~  kon_1_T * CaM0_T * Ca^2 - koff_1_T * CaM2_T +
            -kon_2_T * CaM2_T * Ca^2 + koff_2_T * CaM4_T +
             k_RT_2  * CaM2_R - k_TR_2 * CaM2_T,
        0 ~  kon_2_T * CaM2_T * Ca^2 - koff_2_T * CaM4_T +
             k_RT_4  * CaM4_R - k_TR_4 * CaM4_T,
    CaM_T ~ CaM0_R + CaM2_R + CaM4_R + CaM0_T + CaM2_T + CaM4_T
    ]

    return eqs, vars
end


function solve_Blackwell_TR_eqs(kon_1_T, koff_1_T, kon_2_T, koff_2_T, 
                                kon_1_R, koff_1_R, kon_2_R, koff_2_R, 
                                k_TR_0, k_RT_0, k_TR_2, k_RT_2, k_TR_4, k_RT_4,
                                Ca, CaM_T, CaM_equil
)

    if CaM_equil
        A = [-k_RT_0 - kon_1_R*(Ca^2) koff_1_R 0 k_TR_0 0 0; 
              kon_1_R*(Ca^2) -k_RT_2 - koff_1_R - kon_2_R*(Ca^2) koff_2_R 0 k_TR_2 0; 
              0 kon_2_R*(Ca^2) -k_RT_4 - koff_2_R 0 0 k_TR_4; 
              k_RT_0 0 0 -k_TR_0 - kon_1_T*(Ca^2) koff_1_T 0; 
              0 k_RT_2 0 kon_1_T*(Ca^2) -k_TR_2 - koff_1_T - kon_2_T*(Ca^2) koff_2_T; 
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

    Ca, CaM_T, kon_2C, koff_2C, kon_2N, koff_2N = @variables Ca, CaM_T, kon_2C, koff_2C, kon_2N, koff_2N

    eqs = [
        0 ~ -kon_2C * CaM0  * Ca^2 + koff_2C * CaM2C +
            -kon_2N * CaM0  * Ca^2 + koff_2N * CaM2N,
        0 ~  kon_2C * CaM0  * Ca^2 - koff_2C * CaM2C +
            -kon_2N * CaM2C * Ca^2 + koff_2N * CaM4,
        0 ~  kon_2N * CaM0  * Ca^2 - koff_2N * CaM2N +
            -kon_2C * CaM2N * Ca^2 + koff_2C * CaM4,
        0 ~  kon_2N * CaM2C * Ca^2 - koff_2N * CaM4 +
             kon_2C * CaM2N * Ca^2 - koff_2C * CaM4,
    CaM_T ~ CaM0 + CaM2C + CaM2N + CaM4
    ]

    return eqs, vars
end


function solve_Blackwell_CN_eqs(kon_2C, koff_2C, kon_2N, koff_2N, Ca, CaM_T, CaM_equil)

    if CaM_equil
        A = [
            -kon_2C*(Ca^2) - kon_2N*(Ca^2) koff_2C koff_2N 0; 
             kon_2C*(Ca^2) -koff_2C - kon_2N*(Ca^2) 0 koff_2N; 
            kon_2N*(Ca^2) 0 -koff_2N - kon_2C*(Ca^2) koff_2C; 
       ###  0 kon_2N*(Ca^2) kon_2C*(Ca^2) -koff_2C - koff_2N; 
            1 1 1 1
            ]

        b = [zeros(3); CaM_T]

        sol = A \ b

        return sol
    else
        return [CaM_T; zeros(3)]
    end
end


function get_Pepke_eqs()

    x = @variables CaM0, CaM1C, CaM1N, CaM2C, CaM1N1C, CaM2N, CaM1N2C, CaM2N1C, CaM4;
    CaM_T = (@variables CaM_T)[1];

    Ca, k_off_C1, k_on_C1, k_off_N1, k_on_N1, k_off_C2, k_on_C2, k_off_N2, k_on_N2 =
        @variables Ca, k_off_C1, k_on_C1, k_off_N1, k_on_N1, k_off_C2, k_on_C2, k_off_N2, k_on_N2;

    A = Symbolics.variables(:a, 1:10, 1:9);
    A[1, 1] = -k_on_C1*Ca - k_on_N1*Ca;
    A[1, 2] = k_off_C1;
    A[1, 3] = k_off_N1;
    A[1, 4] = 0;
    A[1, 5] = 0;
    A[1, 6] = 0;
    A[1, 7] = 0;
    A[1, 8] = 0;
    A[1, 9] = 0;
    
    
    A[2, 1] = k_on_C1 * Ca;
    A[2, 2] = -k_off_C1 - k_on_N1*Ca - k_on_C2*Ca;
    A[2, 3] = 0;
    A[2, 4] = k_off_C2;
    A[2, 5] = k_off_N1;
    A[2, 6] = 0;
    A[2, 7] = 0;
    A[2, 8] = 0;
    A[2, 9] = 0;
    
    
    A[3, 1] = k_on_N1 * Ca;
    A[3, 2] = 0 ;
    A[3, 3] = -k_off_N1 - k_on_C1*Ca - k_on_N2*Ca ;
    A[3, 4] = 0;
    A[3, 5] = k_off_C1;
    A[3, 6] = k_off_N2;
    A[3, 7] = 0;
    A[3, 8] = 0;
    A[3, 9] = 0;
    
    
    A[4, 1] = 0;
    A[4, 2] = k_on_C2 * Ca;
    A[4, 3] = 0;
    A[4, 4] = -k_off_C2 - k_on_N1*Ca;
    A[4, 5] = 0;
    A[4, 6] = 0;
    A[4, 7] = k_off_N1 ;
    A[4, 8] = 0;
    A[4, 9] = 0;
    
    
    A[5, 1] = 0;
    A[5, 2] = k_on_N1*Ca;
    A[5, 3] = k_on_C1*Ca;
    A[5, 4] = 0;
    A[5, 5] = -k_off_N1 -k_off_C1 - k_on_C2*Ca - k_on_N2*Ca;
    A[5, 6] = 0;
    A[5, 7] = k_off_C2;
    A[5, 8] = k_off_N2;
    A[5, 9] = 0;
    
    
    A[6, 1] = 0;
    A[6, 2] = 0;
    A[6, 3] = k_on_N2 * Ca;
    A[6, 4] = 0;
    A[6, 5] = 0;
    A[6, 6] = -k_off_N2 - k_on_C1 * Ca;
    A[6, 7] = 0;
    A[6, 8] = k_off_C1;
    A[6, 9] = 0;
    
    
    A[7, 1] = 0;
    A[7, 2] = 0;
    A[7, 3] = 0;
    A[7, 4] = k_on_N1*Ca;
    A[7, 5] = k_on_C2*Ca;
    A[7, 6] = 0;
    A[7, 7] = -k_off_N1 - k_off_C2 - k_on_N2 * Ca;
    A[7, 8] = 0;
    A[7, 9] = k_off_N2;
    
    
    A[8, 1] = 0;
    A[8, 2] = 0;
    A[8, 3] = 0;
    A[8, 4] = 0;
    A[8, 5] = k_on_N2*Ca;
    A[8, 6] = k_on_C1*Ca;
    A[8, 7] = 0;
    A[8, 8] = -k_off_C1 - k_off_N2 - k_on_C2*Ca;
    A[8, 9] = k_off_C2;
    
    
    A[9, 1] = 0;
    A[9, 2] = 0;
    A[9, 3] = 0;
    A[9, 4] = 0;
    A[9, 5] = 0;
    A[9, 6] = 0;
    A[9, 7] = k_on_N2*Ca;
    A[9, 8] = k_on_C2*Ca;
    A[9, 9] = -k_off_N2 - k_off_C2;
    
    
    A[10, 1] = 1;
    A[10, 2] = 1;
    A[10, 3] = 1;
    A[10, 4] = 1;
    A[10, 5] = 1;
    A[10, 6] = 1;
    A[10, 7] = 1;
    A[10, 8] = 1;
    A[10, 9] = 1;

    b = Symbolics.variables(:b, 1:10);
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
    
    eqs = [ (A*x)[1] ~ b[1],
            (A*x)[2] ~ b[2],
            (A*x)[3] ~ b[3],
            (A*x)[4] ~ b[4],
            (A*x)[5] ~ b[5],
            (A*x)[6] ~ b[6],
            (A*x)[7] ~ b[7],
            (A*x)[8] ~ b[8],
            (A*x)[9] ~ b[9],
            (A*x)[10] ~ b[10]];
    
    return A, x, b, eqs
end


function get_Byrne_lobe_eqs()

    vars = @variables N0, N1, N2, N3
    CaM_T = (@variables CaM_t)[1];

    Ca, k_01, k_10, k_02, k_20, k_13, k_31, k_23, k_32 =
        @variables Ca, k01, k10, k02, k20, k13, k31, k23, k32;

    eqs = [
        0     ~ -k01 * Ca * N0 + k10 * N1 +
                -k02 * Ca * N0 + k20 * N2
        0     ~  k01 * Ca * N0 - k10 * N1 +
                -k13 * Ca * N1 + k31 * N3  
        0     ~  k02 * Ca * N0 - k20 * N2 +
                -k23 * Ca * N2 + k32 * N3
        0     ~  k13 * Ca * N1 - k31 * N3 +
                 k23 * Ca * N2 - k32 * N3
        CaM_T ~ N0 + N1 + N2 + N3
    ] 
    
    return eqs, vars
end


function get_Faas_lobe_eqs(Ca_free, Kd_T, Kd_R, CaM_t, CaM_equil)
        
    if CaM_equil
        c1  = (2 * Ca_free) / Kd_T
        c2  = (2 * Kd_R) / Ca_free
        tt₀ = -(c2 * CaM_t) / (1 - (c2+1)*(c1 + 1))
        tr₀ = c1 * tt₀
        rr₀ = tr₀ / c2
        return [tt₀; tr₀;rr₀]
    else
        return [CaM_t; 0.0; 0.0]
    end
end



function get_Stefan_TR_eqs()

    T0, T_A, T_B, T_C, T_D, T_AB, T_AC, T_AD, T_BC, T_BD, T_CD, 
    T_ABC, T_ABD, T_ACD, T_BCD, T_ABCD = @variables T0, T_A, T_B, T_C, T_D, T_AB, T_AC, T_AD, T_BC, T_BD, T_CD, 
    T_ABC, T_ABD, T_ACD, T_BCD, T_ABCD

    R0, R_A, R_B, R_C, R_D, R_AB, R_AC, R_AD, R_BC, R_BD, R_CD, 
    R_ABC, R_ABD, R_ACD, R_BCD, R_ABCD = @variables R0, R_A, R_B, R_C, R_D, R_AB, R_AC, R_AD, R_BC, R_BD, R_CD, 
    R_ABC, R_ABD, R_ACD, R_BCD, R_ABCD

    kf_AT, kb_AT, kf_BT, kb_BT, kf_CT, kb_CT, kf_DT, kb_DT = @variables kf_AT, kb_AT, kf_BT, kb_BT, kf_CT, kb_CT, kf_DT, kb_DT
    kf_AR, kb_AR, kf_BR, kb_BR, kf_CR, kb_CR, kf_DR, kb_DR = @variables kf_AR, kb_AR, kf_BR, kb_BR, kf_CR, kb_CR, kf_DR, kb_DR

    k0_TR, k0_RT, k1_TR, k1_RT, k2_TR, k2_RT, 
    k3_TR, k3_RT, k4_TR, k4_RT = @variables k0_TR, k0_RT, k1_TR, k1_RT, k2_TR, k2_RT, 
    k3_TR, k3_RT, k4_TR, k4_RT

    Ca, CaM_T = @variables Ca, CaM_T

    eqs = [
    0 ~ -kf_AT * T0 * Ca + kb_AT * T_A +  ### T0 ### start of T states
    -kf_BT * T0 * Ca + kb_BT * T_B +
    -kf_CT * T0 * Ca + kb_CT * T_C +
    -kf_DT * T0 * Ca + kb_DT * T_D +
    -T0 * k0_TR + R0 * k0_RT,
    0 ~  kf_AT * T0 * Ca - kb_AT * T_A +  ### T_A
    -kf_BT * T_A * Ca + kb_BT * T_AB +
    -kf_CT * T_A * Ca + kb_CT * T_AC +
    -kf_DT * T_A * Ca + kb_DT * T_AD +
    -T_A * k1_TR + R_A * k1_RT,
    0 ~  kf_BT * T0 * Ca - kb_BT * T_B +  ### T_B
    -kf_AT * T_B * Ca + kb_AT * T_AB +
    -kf_CT * T_B * Ca + kb_CT * T_BC +
    -kf_DT * T_B * Ca + kb_DT * T_BD +
    -T_B * k1_TR + R_B * k1_RT,
    0 ~  kf_CT * T0 * Ca - kb_CT * T_C +  ### T_C
    -kf_AT * T_C * Ca + kb_AT * T_AC +
    -kf_BT * T_C * Ca + kb_BT * T_BC +
    -kf_DT * T_C * Ca + kb_DT * T_CD +
    -T_C * k1_TR + R_C * k1_RT,
    0 ~  kf_DT * T0 * Ca - kb_DT * T_D +  ### T_D
    -kf_AT * T_D * Ca + kb_AT * T_AD +
    -kf_BT * T_D * Ca + kb_BT * T_BD +
    -kf_CT * T_D * Ca + kb_CT * T_CD +
    -T_D * k1_TR + R_D * k1_RT,
    0 ~  kf_AT * T_B * Ca - kb_AT * T_AB + ### T_AB
     kf_BT * T_A * Ca - kb_BT * T_AB +
    -kf_CT * T_AB * Ca + kb_CT * T_ABC +
    -kf_DT * T_AB * Ca + kb_DT * T_ABD +
    -T_AB * k2_TR + R_AB * k2_RT,
    0 ~  kf_AT * T_C * Ca - kb_AT * T_AC + ### T_AC
     kf_CT * T_A * Ca - kb_CT * T_AC +
    -kf_BT * T_AC * Ca + kb_BT * T_ABC +
    -kf_DT * T_AC * Ca + kb_DT * T_ACD +
    -T_AC * k2_TR + R_AC * k2_RT,
    0 ~  kf_AT * T_D * Ca - kb_AT * T_AD + ### T_AD
     kf_DT * T_A * Ca - kb_DT * T_AD +
    -kf_BT * T_AD * Ca + kb_BT * T_ABC +
    -kf_CT * T_AD * Ca + kb_CT * T_ACD +
    -T_AD * k2_TR + R_AD * k2_RT,
    0 ~  kf_BT * T_C * Ca - kb_BT * T_BC + ### T_BC
     kf_CT * T_B * Ca - kb_CT * T_BC +
    -kf_AT * T_BC * Ca + kb_AT * T_ABC +
    -kf_DT * T_BC * Ca + kb_DT * T_BCD +
    -T_BC * k2_TR + R_BC * k2_RT,
    0 ~  kf_BT * T_D * Ca - kb_BT * T_BD + ### T_BD
     kf_DT * T_B * Ca - kb_DT * T_BD +
    -kf_AT * T_BC * Ca + kb_AT * T_ABC +
    -kf_DT * T_BC * Ca + kb_DT * T_BCD +
    -T_BD * k2_TR + R_BD * k2_RT,
    0 ~  kf_CT * T_D * Ca - kb_CT * T_CD + ### T_CD
     kf_DT * T_C * Ca - kb_DT * T_CD +
    -kf_AT * T_BC * Ca + kb_AT * T_ABC +
    -kf_DT * T_BC * Ca + kb_DT * T_BCD +
    -T_CD * k2_TR + R_CD * k2_RT,
    0 ~  kf_AT * T_BC * Ca - kb_AT * T_ABC + ### T_ABC
     kf_BT * T_AC * Ca - kb_BT * T_ABC +
     kf_CT * T_AB * Ca - kb_CT * T_ABC +
    -kf_DT * T_ABC * Ca + kb_DT * T_ABCD +
    -T_ABC * k3_TR + R_ABC * k1_RT,
    0 ~  kf_AT * T_BD * Ca - kb_AT * T_ABD + ### T_ABD
     kf_BT * T_AD * Ca - kb_BT * T_ABD +
     kf_DT * T_AB * Ca - kb_DT * T_ABD +
    -kf_CT * T_ABD * Ca + kb_CT * T_ABCD +
    -T_ABD * k3_TR + R_ABD * k3_RT,
    0 ~  kf_AT * T_CD * Ca - kb_AT * T_ACD + ### T_ACD
     kf_CT * T_AD * Ca - kb_CT * T_ACD +
     kf_DT * T_AC * Ca - kb_DT * T_ACD +
    -kf_BT * T_ACD * Ca + kb_BT * T_ABCD +
    -T_ACD * k3_TR + R_ACD * k3_RT,
    0 ~  kf_BT * T_CD * Ca - kb_BT * T_ACD + ### T_BCD
     kf_CT * T_BD * Ca - kb_CT * T_ACD +
     kf_DT * T_BC * Ca - kb_DT * T_ACD +
    -kf_AT * T_BCD * Ca + kb_AT * T_ABCD +
    -T_BCD * k3_TR + R_BCD * k3_RT,
    0 ~  kf_AT * T_BCD * Ca - kb_AT * T_ABCD + ### T_ABCD
     kf_BT * T_ACD * Ca - kb_BT * T_ABCD +
     kf_CT * T_ABD * Ca - kb_CT * T_ABCD +
     kf_DT * T_ABC * Ca - kb_DT * T_ABCD +
    -T_ABCD * k4_TR + R_ABCD * k4_RT,
    0 ~ -kf_AR * R0 * Ca + kb_AR * R_A + ### R0        #### start of R states
    -kf_BR * R0 * Ca + kb_BR * R_B +
    -kf_CR * R0 * Ca + kb_CR * R_C +
    -kf_DR * R0 * Ca + kb_DR * R_D +
     T0 * k0_TR - R0 * k0_RT,
    0 ~  kf_AR * R0 * Ca - kb_AR * R_A + ### R_A
    -kf_BR * R_A * Ca + kb_BR * R_AB +
    -kf_CR * R_A * Ca + kb_CR * R_AC +
    -kf_DR * R_A * Ca + kb_DR * R_AD +
     T_A * k1_TR - R_A * k1_RT,
    0 ~  kf_BR * R0 * Ca - kb_BR * R_B + ### R_B
    -kf_AR * R_B * Ca + kb_AR * R_AB +
    -kf_CR * R_B * Ca + kb_CR * R_BC +
    -kf_DR * R_B * Ca + kb_DR * R_BD +
     T_B * k1_TR - R_B * k1_RT,
    0 ~  kf_CR * R0 * Ca - kb_CR * R_C + ### R_C
    -kf_AR * R_C * Ca + kb_AR * R_AC +
    -kf_BR * R_C * Ca + kb_BR * R_BC +
    -kf_DR * R_C * Ca + kb_DR * R_CD +
     T_C * k1_TR - R_C * k1_RT,
    0 ~  kf_DR * R0 * Ca - kb_DR * R_D + ### R_D
    -kf_AR * R_D * Ca + kb_AR * R_AD +
    -kf_BR * R_D * Ca + kb_BR * R_BD +
    -kf_CR * R_D * Ca + kb_CR * R_CD +
     T_D * k1_TR - R_D * k1_RT,
    0 ~  kf_AR * R_B * Ca - kb_AR * R_AB + ### R_AB
     kf_BR * R_A * Ca - kb_BR * R_AB +
    -kf_CR * R_AB * Ca + kb_CR * R_ABC +
    -kf_DR * R_AB * Ca + kb_DR * R_ABD +
     T_AB * k2_TR - R_AB * k2_RT,
    0 ~  kf_AR * R_C * Ca - kb_AR * R_AC + ### R_AC
     kf_CR * R_A * Ca - kb_CR * R_AC +
    -kf_BR * R_AC * Ca + kb_BR * R_ABC +
    -kf_DR * R_AC * Ca + kb_DR * R_ACD +
     T_AC * k2_TR - R_AC * k2_RT,
    0 ~  kf_AR * R_D * Ca - kb_AR * R_AD + ### R_AD
     kf_DR * R_A * Ca - kb_DR * R_AD +
    -kf_BR * R_AD * Ca + kb_BR * R_ABC +
    -kf_CR * R_AD * Ca + kb_CR * R_ACD +
     T_AD * k2_TR - R_AD * k2_RT,
    0 ~  kf_BR * R_C * Ca - kb_BR * R_BC + ### R_BC
     kf_CR * R_B * Ca - kb_CR * R_BC +
    -kf_AR * R_BC * Ca + kb_AR * R_ABC +
    -kf_DR * R_BC * Ca + kb_DR * R_BCD +
     T_BC * k2_TR - R_BC * k2_RT,
    0 ~  kf_BR * R_D * Ca - kb_BR * R_BD + ### R_BD
     kf_DR * R_B * Ca - kb_DR * R_BD +
    -kf_AR * R_BC * Ca + kb_AR * R_ABC +
    -kf_DR * R_BC * Ca + kb_DR * R_BCD +
     T_BD * k2_TR - R_BD * k2_RT,
    0 ~  kf_CR * R_D * Ca - kb_CR * R_CD + ### R_CD
     kf_DR * R_C * Ca - kb_DR * R_CD +
    -kf_AR * R_BC * Ca + kb_AR * R_ABC +
    -kf_DR * R_BC * Ca + kb_DR * R_BCD +
     T_CD * k2_TR - R_CD * k2_RT,
    0 ~  kf_AR * R_BC * Ca - kb_AR * R_ABC + ### R_ABC
     kf_BR * R_AC * Ca - kb_BR * R_ABC +
     kf_CR * R_AB * Ca - kb_CR * R_ABC +
    -kf_DR * R_ABC * Ca + kb_DR * R_ABCD +
     T_ABC * k3_TR - R_ABC * k3_RT,
    0 ~  kf_AR * R_BD * Ca - kb_AR * R_ABD + ### R_ABD
     kf_BR * R_AD * Ca - kb_BR * R_ABD +
     kf_DR * R_AB * Ca - kb_DR * R_ABD +
    -kf_CR * R_ABD * Ca + kb_CR * R_ABCD +
     T_ABD * k3_TR - R_ABD * k3_RT,
    0 ~  kf_AR * R_CD * Ca - kb_AR * R_ACD + ### R_ACD
     kf_CR * R_AD * Ca - kb_CR * R_ACD +
     kf_DR * R_AC * Ca - kb_DR * R_ACD +
    -kf_BR * R_ACD * Ca + kb_BR * R_ABCD +
     T_ACD * k3_TR - R_ACD * k3_RT,
    0 ~  kf_BR * R_CD * Ca - kb_BR * R_ACD + ### R_BCD
     kf_CR * R_BD * Ca - kb_CR * R_ACD +
     kf_DR * R_BC * Ca - kb_DR * R_ACD +
    -kf_AR * R_BCD * Ca + kb_AR * R_ABCD +
     T_BCD * k3_TR - R_BCD * k3_RT,
    0 ~  kf_AR * R_BCD * Ca - kb_AR * R_ABCD + ### R_ABCD
     kf_BR * R_ACD * Ca - kb_BR * R_ABCD +
     kf_CR * R_ABD * Ca - kb_CR * R_ABCD +
     kf_DR * R_ABC * Ca - kb_DR * R_ABCD +
     T_ABCD * k4_TR + R_ABCD * k4_RT,
    CaM_T ~ T0 + T_A + T_B + T_C + T_D + T_AB + T_AC + T_AD + T_BC + T_BD + T_CD + T_ABC + T_ABD + T_ACD + T_BCD + T_ABCD + 
            R0 + R_A + R_B + R_C + R_D + R_AB + R_AC + R_AD + R_BC + R_BD + R_CD + R_ABC + R_ABD + R_ACD + R_BCD + R_ABCD
    ]
    vars = [T0, T_A, T_B, T_C, T_D, T_AB, T_AC, T_AD, T_BC, T_BD, T_CD, T_ABC, T_ABD, T_ACD, T_BCD, T_ABCD,
            R0, R_A, R_B, R_C, R_D, R_AB, R_AC, R_AD, R_BC, R_BD, R_CD, R_ABC, R_ABD, R_ACD, R_BCD, R_ABCD]
    return eqs, vars
end



function Stefan_TR_model_equilibriums(m_params)
    outs = zeros(16, size(m_params)[2])
    A = zeros(16, 16)
    A[16,:] .= 1
    b = zeros(16)
    for k in 1:size(m_params)[2]
        kf_AT, kb_AT, kf_BT, kb_BT, kf_CT, kb_CT, kf_DT, kb_DT, Ca, CaM_t = m_params[:,k]

        A[1, 1] = -Ca*kf_AT - Ca*kf_BT - Ca*kf_CT - Ca*kf_DT
        A[1, 2] = kb_AT 
        A[1, 3] = kb_BT 
        A[1, 4] = kb_CT 
        A[1, 5] = kb_DT 

        A[2, 1] = Ca*kf_AT
        A[2, 2] = -kb_AT - Ca*kf_BT - Ca*kf_CT - Ca*kf_DT
        A[2, 6] = kb_BT
        A[2, 7] = kb_CT
        A[2, 8] = kb_DT

        A[3, 1] = Ca*kf_BT
        A[3, 3] = -kb_BT - Ca*kf_AT - Ca*kf_CT - Ca*kf_DT
        A[3, 6] = kb_AT
        A[3, 9] = kb_CT
        A[3,10] = kb_DT

        A[4, 1] = Ca*kf_CT
        A[4, 4] = -kb_CT - Ca*kf_AT - Ca*kf_BT - Ca*kf_DT
        A[4, 7] = kb_AT
        A[4, 9] = kb_BT
        A[4,11] = kb_DT

        A[5, 1] = Ca*kf_DT
        A[5, 5] = -kb_DT - Ca*kf_AT - Ca*kf_BT - Ca*kf_CT
        A[5, 8] = kb_AT
        A[5,10] = kb_BT
        A[5,11] = kb_CT

        A[6, 2] = Ca*kf_BT
        A[6, 3] = Ca*kf_AT
        A[6, 6] = -kb_AT - kb_BT - Ca*kf_CT - Ca*kf_DT
        A[6,12] = kb_CT
        A[6,13] = kb_DT

        A[7, 2] = Ca*kf_CT
        A[7, 4] = Ca*kf_AT
        A[7, 7] = -kb_AT - kb_CT - Ca*kf_BT - Ca*kf_DT
        A[7,12] = kb_BT
        A[7,14] = kb_DT

        A[8, 2] = Ca*kf_DT
        A[8, 5] = Ca*kf_AT
        A[8, 8] = -kb_AT - kb_DT - Ca*kf_BT - Ca*kf_CT
        A[8,12] = kb_BT
        A[8,14] = kb_CT

        A[9, 3] = Ca*kf_CT
        A[9, 4] = Ca*kf_BT
        A[9, 9] = -kb_BT - kb_CT - Ca*kf_AT - Ca*kf_DT
        A[9,12] = kb_AT
        A[9,15] = kb_DT

        A[10, 3] = Ca*kf_DT
        A[10, 5] = Ca*kf_BT
        A[10, 9] = -Ca*kf_AT - Ca*kf_DT
        A[10,10] = -kb_BT - kb_DT
        A[10,12] = kb_AT
        A[10,15] = kb_DT

        A[11, 4] = Ca*kf_DT
        A[11, 5] = Ca*kf_CT
        A[11, 9] = -Ca*kf_AT - Ca*kf_DT
        A[11,11] = -kb_CT - kb_DT
        A[11,12] = kb_AT
        A[11,15] = kb_DT

        A[12, 6] = Ca*kf_CT
        A[12, 7] = Ca*kf_BT
        A[12, 9] = Ca*kf_AT
        A[12,12] = -kb_AT - kb_BT - kb_CT - Ca*kf_DT
        A[12,16] = kb_DT

        A[13, 6] = Ca*kf_DT
        A[13, 8] = Ca*kf_BT
        A[13,10] = Ca*kf_AT
        A[13,13] = -kb_AT - kb_BT - kb_DT - Ca*kf_CT
        A[13,16] = kb_CT

        A[14, 7] = Ca*kf_DT
        A[14, 8] = Ca*kf_CT
        A[14,11] = Ca*kf_AT
        A[14,14] = -kb_AT - kb_CT - kb_DT - Ca*kf_BT
        A[14,16] = kb_BT

        A[15, 9] = Ca*kf_DT
        A[15,10] = Ca*kf_CT
        A[15,11] = Ca*kf_BT
        A[15,14] = -kb_BT - kb_CT - kb_DT
        A[15,15] = -Ca*kf_AT
        A[15,16] = kb_AT

        b[end] = CaM_t
        sol = A\b
        outs[:,k] .= sol
    end
    return outs
end


function solve_Stefan_TR_eq(Kd_AT, kf_AT, kb_AT, Kd_BT, kf_BT, kb_BT, Kd_CT, kf_CT, kb_CT, Kd_DT, kf_DT, kb_DT, 
                            Kd_AR, kf_AR, kb_AR, Kd_BR, kf_BR, kb_BR, Kd_CR, kf_CR, kb_CR, Kd_DR, kf_DR, kb_DR, 
                            Kd_0_TR, k0_TR, k0_RT, Kd_1_TR, k1_TR, k1_RT, Kd_2_TR, k2_TR, k2_RT, 
                            Kd_3_TR, k3_TR, k3_RT, Kd_4_TR, k4_TR, k4_RT, Ca_free, CaM_t)

    A = [-k0_TR - Ca_free*kf_AT - Ca_free*kf_BT - Ca_free*kf_CT - Ca_free*kf_DT kb_AT kb_BT kb_CT kb_DT 0 0 0 0 0 0 0 0 0 0 0 k0_RT 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0; 
    Ca_free*kf_AT -k1_TR - kb_AT - Ca_free*kf_BT - Ca_free*kf_CT - Ca_free*kf_DT 0 0 0 kb_BT kb_CT kb_DT 0 0 0 0 0 0 0 0 0 k1_RT 0 0 0 0 0 0 0 0 0 0 0 0 0 0; 
    Ca_free*kf_BT 0 -k1_TR - kb_BT - Ca_free*kf_AT - Ca_free*kf_CT - Ca_free*kf_DT 0 0 kb_AT 0 0 kb_CT kb_DT 0 0 0 0 0 0 0 0 k1_RT 0 0 0 0 0 0 0 0 0 0 0 0 0; 
    Ca_free*kf_CT 0 0 -k1_TR - kb_CT - Ca_free*kf_AT - Ca_free*kf_BT - Ca_free*kf_DT 0 0 kb_AT 0 kb_BT 0 kb_DT 0 0 0 0 0 0 0 0 k1_RT 0 0 0 0 0 0 0 0 0 0 0 0; 
    Ca_free*kf_DT 0 0 0 -k1_TR - kb_DT - Ca_free*kf_AT - Ca_free*kf_BT - Ca_free*kf_CT 0 0 kb_AT 0 kb_BT kb_CT 0 0 0 0 0 0 0 0 0 k1_RT 0 0 0 0 0 0 0 0 0 0 0; 
    0 Ca_free*kf_BT Ca_free*kf_AT 0 0 -k2_TR - kb_AT - kb_BT - Ca_free*kf_CT - Ca_free*kf_DT 0 0 0 0 0 kb_CT kb_DT 0 0 0 0 0 0 0 0 k2_RT 0 0 0 0 0 0 0 0 0 0; 
    0 Ca_free*kf_CT 0 Ca_free*kf_AT 0 0 -k2_TR - kb_AT - kb_CT - Ca_free*kf_BT - Ca_free*kf_DT 0 0 0 0 kb_BT 0 kb_DT 0 0 0 0 0 0 0 0 k2_RT 0 0 0 0 0 0 0 0 0; 
    0 Ca_free*kf_DT 0 0 Ca_free*kf_AT 0 0 -k2_TR - kb_AT - kb_DT - Ca_free*kf_BT - Ca_free*kf_CT 0 0 0 kb_BT 0 kb_CT 0 0 0 0 0 0 0 0 0 k2_RT 0 0 0 0 0 0 0 0; 
    0 0 Ca_free*kf_CT Ca_free*kf_BT 0 0 0 0 -k2_TR - kb_BT - kb_CT - Ca_free*kf_AT - Ca_free*kf_DT 0 0 kb_AT 0 0 kb_DT 0 0 0 0 0 0 0 0 0 k2_RT 0 0 0 0 0 0 0; 
    0 0 Ca_free*kf_DT 0 Ca_free*kf_BT 0 0 0 -Ca_free*kf_AT - Ca_free*kf_DT -k2_TR - kb_BT - kb_DT 0 kb_AT 0 0 kb_DT 0 0 0 0 0 0 0 0 0 0 k2_RT 0 0 0 0 0 0; 
    0 0 0 Ca_free*kf_DT Ca_free*kf_CT 0 0 0 -Ca_free*kf_AT - Ca_free*kf_DT 0 -k2_TR - kb_CT - kb_DT kb_AT 0 0 kb_DT 0 0 0 0 0 0 0 0 0 0 0 k2_RT 0 0 0 0 0; 
    0 0 0 0 0 Ca_free*kf_CT Ca_free*kf_BT 0 Ca_free*kf_AT 0 0 -k3_TR - kb_AT - kb_BT - kb_CT - Ca_free*kf_DT 0 0 0 kb_DT 0 0 0 0 0 0 0 0 0 0 0 k1_RT 0 0 0 0; 
    0 0 0 0 0 Ca_free*kf_DT 0 Ca_free*kf_BT 0 Ca_free*kf_AT 0 0 -k3_TR - kb_AT - kb_BT - kb_DT - Ca_free*kf_CT 0 0 kb_CT 0 0 0 0 0 0 0 0 0 0 0 0 k3_RT 0 0 0; 
    0 0 0 0 0 0 Ca_free*kf_DT Ca_free*kf_CT 0 0 Ca_free*kf_AT 0 0 -k3_TR - kb_AT - kb_CT - kb_DT - Ca_free*kf_BT 0 kb_BT 0 0 0 0 0 0 0 0 0 0 0 0 0 k3_RT 0 0; 
    0 0 0 0 0 0 0 0 Ca_free*kf_DT Ca_free*kf_CT Ca_free*kf_BT 0 0 -kb_BT - kb_CT - kb_DT -k3_TR - Ca_free*kf_AT kb_AT 0 0 0 0 0 0 0 0 0 0 0 0 0 0 k3_RT 0; 
    0 0 0 0 0 0 0 0 0 0 0 Ca_free*kf_DT Ca_free*kf_CT Ca_free*kf_BT Ca_free*kf_AT -k4_TR - kb_AT - kb_BT - kb_CT - kb_DT 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 k4_RT; 
    k0_TR 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -k0_RT - Ca_free*kf_AR - Ca_free*kf_BR - Ca_free*kf_CR - Ca_free*kf_DR kb_AR kb_BR kb_CR kb_DR 0 0 0 0 0 0 0 0 0 0 0; 
    0 k1_TR 0 0 0 0 0 0 0 0 0 0 0 0 0 0 Ca_free*kf_AR -k1_RT - kb_AR - Ca_free*kf_BR - Ca_free*kf_CR - Ca_free*kf_DR 0 0 0 kb_BR kb_CR kb_DR 0 0 0 0 0 0 0 0; 
    0 0 k1_TR 0 0 0 0 0 0 0 0 0 0 0 0 0 Ca_free*kf_BR 0 -k1_RT - kb_BR - Ca_free*kf_AR - Ca_free*kf_CR - Ca_free*kf_DR 0 0 kb_AR 0 0 kb_CR kb_DR 0 0 0 0 0 0; 
    0 0 0 k1_TR 0 0 0 0 0 0 0 0 0 0 0 0 Ca_free*kf_CR 0 0 -k1_RT - kb_CR - Ca_free*kf_AR - Ca_free*kf_BR - Ca_free*kf_DR 0 0 kb_AR 0 kb_BR 0 kb_DR 0 0 0 0 0; 
    0 0 0 0 k1_TR 0 0 0 0 0 0 0 0 0 0 0 Ca_free*kf_DR 0 0 0 -k1_RT - kb_DR - Ca_free*kf_AR - Ca_free*kf_BR - Ca_free*kf_CR 0 0 kb_AR 0 kb_BR kb_CR 0 0 0 0 0; 
    0 0 0 0 0 k2_TR 0 0 0 0 0 0 0 0 0 0 0 Ca_free*kf_BR Ca_free*kf_AR 0 0 -k2_RT - kb_AR - kb_BR - Ca_free*kf_CR - Ca_free*kf_DR 0 0 0 0 0 kb_CR kb_DR 0 0 0; 
    0 0 0 0 0 0 k2_TR 0 0 0 0 0 0 0 0 0 0 Ca_free*kf_CR 0 Ca_free*kf_AR 0 0 -k2_RT - kb_AR - kb_CR - Ca_free*kf_BR - Ca_free*kf_DR 0 0 0 0 kb_BR 0 kb_DR 0 0; 
    0 0 0 0 0 0 0 k2_TR 0 0 0 0 0 0 0 0 0 Ca_free*kf_DR 0 0 Ca_free*kf_AR 0 0 -k2_RT - kb_AR - kb_DR - Ca_free*kf_BR - Ca_free*kf_CR 0 0 0 kb_BR 0 kb_CR 0 0; 
    0 0 0 0 0 0 0 0 k2_TR 0 0 0 0 0 0 0 0 0 Ca_free*kf_CR Ca_free*kf_BR 0 0 0 0 -k2_RT - kb_BR - kb_CR - Ca_free*kf_AR - Ca_free*kf_DR 0 0 kb_AR 0 0 kb_DR 0; 
    0 0 0 0 0 0 0 0 0 k2_TR 0 0 0 0 0 0 0 0 Ca_free*kf_DR 0 Ca_free*kf_BR 0 0 0 -Ca_free*kf_AR - Ca_free*kf_DR -k2_RT - kb_BR - kb_DR 0 kb_AR 0 0 kb_DR 0; 
    0 0 0 0 0 0 0 0 0 0 k2_TR 0 0 0 0 0 0 0 0 Ca_free*kf_DR Ca_free*kf_CR 0 0 0 -Ca_free*kf_AR - Ca_free*kf_DR 0 -k2_RT - kb_CR - kb_DR kb_AR 0 0 kb_DR 0; 
    0 0 0 0 0 0 0 0 0 0 0 k3_TR 0 0 0 0 0 0 0 0 0 Ca_free*kf_CR Ca_free*kf_BR 0 Ca_free*kf_AR 0 0 -k3_RT - kb_AR - kb_BR - kb_CR - Ca_free*kf_DR 0 0 0 kb_DR;
    0 0 0 0 0 0 0 0 0 0 0 0 k3_TR 0 0 0 0 0 0 0 0 Ca_free*kf_DR 0 Ca_free*kf_BR 0 Ca_free*kf_AR 0 0 -k3_RT - kb_AR - kb_BR - kb_DR - Ca_free*kf_CR 0 0 kb_CR; 
    0 0 0 0 0 0 0 0 0 0 0 0 0 k3_TR 0 0 0 0 0 0 0 0 Ca_free*kf_DR Ca_free*kf_CR 0 0 Ca_free*kf_AR 0 0 -k3_RT - kb_AR - kb_CR - kb_DR - Ca_free*kf_BR 0 kb_BR; 
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 k3_TR 0 0 0 0 0 0 0 0 0 Ca_free*kf_DR Ca_free*kf_CR Ca_free*kf_BR 0 0 -kb_BR - kb_CR - kb_DR -k3_RT - Ca_free*kf_AR kb_AR;
###              0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 k4_TR 0 0 0 0 0 0 0 0 0 0 0 Ca*kf_DR Ca*kf_CR Ca*kf_BR Ca*kf_AR k4_RT - kb_AR - kb_BR - kb_CR - kb_DR;
    1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1] 

    b = [zeros(31); CaM_t]

    sol = A\b
    return sol
end


function get_Pepke_TR_eqs()

    x = @variables T0, T1C, T1N, T2C, T1N1C, T2N, T1N2C, T2N1C, T4, R0, R1C, R1N, R2C, R1N1C, R2N, R1N2C, R2N1C, R4;
    Ca, CaM_T = @variables Ca, CaM_T;

    k_off_C1_T, k_on_C1_T, k_off_N1_T, k_on_N1_T, k_off_C2_T, k_on_C2_T, k_off_N2_T, k_on_N2_T =
    @variables k_off_C1_T, k_on_C1_T, k_off_N1_T, k_on_N1_T, k_off_C2_T, k_on_C2_T, k_off_N2_T, k_on_N2_T;

    k_off_C1_R, k_on_C1_R, k_off_N1_R, k_on_N1_R, k_off_C2_R, k_on_C2_R, k_off_N2_R, k_on_N2_R =
    @variables k_off_C1_R, k_on_C1_R, k_off_N1_R, k_on_N1_R, k_off_C2_R, k_on_C2_R, k_off_N2_R, k_on_N2_R;

    k0_TR, k0_RT, k1_TR, k1_RT, k2_TR, k2_RT, 
    k3_TR, k3_RT, k4_TR, k4_RT = @variables k0_TR, k0_RT, k1_TR, k1_RT, k2_TR, k2_RT, 
    k3_TR, k3_RT, k4_TR, k4_RT

    eqs =[
0       ~ -k_on_C1_T * T0 * Ca    + k_off_C1_T * T1C +    ### start of T states
    -k_on_N1_T * T0 * Ca    + k_off_N1_T * T1N +
    -k0_TR * T0     + k0_RT*R0,
0      ~  k_on_C1_T * T0 * Ca    - k_off_C1_T * T1C +
    -k_on_C2_T * T1C * Ca   + k_off_C2_T * T2C +
    -k_on_N1_T * T1C * Ca   + k_off_N1_T * T1N1C +
    -k1_TR * T1C    + k1_RT*R1C,
0      ~  k_on_N1_T * T0 * Ca    - k_off_N1_T * T1N +
    -k_on_N2_T * T1N * Ca   + k_off_N2_T * T2N +
    -k_on_C1_T * T1N * Ca   + k_off_C1_T * T1N1C +
    -k1_TR * T1N    + k1_RT*R1N,
0      ~  k_on_C2_T * T1C * Ca   - k_off_C2_T * T2C +
    -k_on_N1_T * T2C * Ca   + k_off_N1_T * T1N2C +
    -k2_TR * T2C    + k2_RT*R2C,
0    ~  k_on_N1_T * T1C * Ca   - k_off_N1_T * T1N1C +
     k_on_C1_T * T1N * Ca   - k_off_C1_T * T1N1C +
    -k_on_C2_T * T1N1C * Ca + k_off_C2_T * T1N2C +
    -k_on_N2_T * T1N1C * Ca + k_off_N2_T * T2N1C +
    -k2_TR * T1N1C  + k2_RT*R1N1C,
0      ~  k_on_N2_T * T1N * Ca   - k_off_N2_T * T2N +
    -k_on_C1_T * T2N * Ca   + k_off_C1_T * T2N1C +
    -k2_TR * T2N  + k2_RT*R2N,
0    ~  k_on_N1_T * T2C * Ca   - k_off_N1_T * T1N2C +
     k_on_C2_T * T1N1C * Ca - k_off_C2_T * T1N2C +
    -k_on_N2_T * T1N2C * Ca + k_off_N2_T * T4 +
    -k3_TR * T1N2C  + k3_RT*R1N2C,
0    ~  k_on_C1_T * T2N * Ca   - k_off_C1_T * T2N1C +
     k_on_N2_T * T1N1C * Ca - k_off_N2_T * T2N1C +
    -k_on_C2_T * T2N1C * Ca + k_off_C2_T * T4 +
    -k3_TR * T2N1C  + k3_RT*R2N1C,
0       ~  k_on_N2_T * T1N2C * Ca - k_off_N2_T * T4 +
     k_on_C2_T * T2N1C * Ca - k_off_C2_T * T4 +
    -k4_TR * T4     + k4_RT*R4,
0       ~ -k_on_C1_R * R0 * Ca    + k_off_C1_R * R1C +            ### start of R states
    -k_on_N1_R * R0 * Ca    + k_off_N1_R * R1N +
     k0_TR * T0     - k0_RT*R0,
0      ~  k_on_C1_R * R0 * Ca    - k_off_C1_R * R1C +
    -k_on_C2_R * R1C * Ca   + k_off_C2_R * R2C +
    -k_on_N1_R * R1C * Ca   + k_off_N1_R * R1N1C +
     k1_TR * T1C    - k1_RT*R1C,
0      ~  k_on_N1_R * R0 * Ca    - k_off_N1_R * R1N +
    -k_on_N2_R * R1N * Ca   + k_off_N2_R * R2N +
    -k_on_C1_R * R1N * Ca   + k_off_C1_R * R1N1C +
     k1_TR * T1N    - k1_RT*R1N,
0      ~  k_on_C2_R * R1C * Ca   - k_off_C2_R * R2C +
    -k_on_N1_R * R2C * Ca   + k_off_N1_R * R1N2C +
     k2_TR * T2C    - k2_RT*R2C,
0    ~  k_on_N1_R * R1C * Ca   - k_off_N1_R * R1N1C +
     k_on_C1_R * R1N * Ca   - k_off_C1_R * R1N1C +
    -k_on_C2_R * R1N1C * Ca + k_off_C2_R * R1N2C +
    -k_on_N2_R * R1N1C * Ca + k_off_N2_R * R2N1C +
     k2_TR * T1N1C  - k2_RT*R1N1C,
0      ~  k_on_N2_R * R1N * Ca   - k_off_N2_R * R2N +
    -k_on_C1_R * R2N * Ca   + k_off_C1_R * R2N1C +
     k2_TR * T2N  - k2_RT*R2N,
0    ~  k_on_N1_R * R2C * Ca   - k_off_N1_R * R1N2C +
     k_on_C2_R * R1N1C * Ca - k_off_C2_R * R1N2C +
    -k_on_N2_R * R1N2C * Ca + k_off_N2_R * R4 +
     k3_TR * T1N2C  - k3_RT*R1N2C,
0    ~  k_on_C1_R * R2N * Ca   - k_off_C1_R * R2N1C +
     k_on_N2_R * R1N1C * Ca - k_off_N2_R * R2N1C +
    -k_on_C2_R * R2N1C * Ca + k_off_C2_R * R4 +
     k3_TR * T2N1C  - k3_RT*R2N1C,
0       ~  k_on_N2_R * R1N2C * Ca - k_off_N2_R * R4 +
     k_on_C2_R * R2N1C * Ca - k_off_C2_R * R4 +
     k4_TR * T4     - k4_RT*R4,
    CaM_T ~ T0 + T1C + T1N + T2C + T1N1C + T2N + T1N2C + T2N1C + T4 + R0 + R1C + R1N + R2C + R1N1C + R2N + R1N2C + R2N1C + R4
    ]

    return eqs, x
end


function solve_Pepke_TR_eq(k_on_C1_T, k_off_C1_T, k_on_C2_T, k_off_C2_T, k_on_N1_T, k_off_N1_T, k_on_N2_T, k_off_N2_T, 
                           k_on_C1_R, k_off_C1_R, k_on_C2_R, k_off_C2_R, k_on_N1_R, k_off_N1_R, k_on_N2_R, k_off_N2_R,
                           Kd_0_TR, k0_TR, k0_RT, Kd_1_TR, k1_TR, k1_RT, Kd_2_TR, k2_TR, k2_RT, 
                           Kd_3_TR, k3_TR, k3_RT, Kd_4_TR, k4_TR, k4_RT, Ca, CaM_t)
    A = [
         -k0_TR - Ca*k_on_C1_T - Ca*k_on_N1_T k_off_C1_T k_off_N1_T 0 0 0 0 0 0 k0_RT 0 0 0 0 0 0 0 0; 
          Ca*k_on_C1_T -k1_TR - k_off_C1_T - Ca*k_on_C2_T - Ca*k_on_N1_T 0 k_off_C2_T k_off_N1_T 0 0 0 0 0 k1_RT 0 0 0 0 0 0 0; 
          Ca*k_on_N1_T 0 -k1_TR - k_off_N1_T - Ca*k_on_C1_T - Ca*k_on_N2_T 0 k_off_C1_T k_off_N2_T 0 0 0 0 0 k1_RT 0 0 0 0 0 0; 
          0 Ca*k_on_C2_T 0 -k2_TR - k_off_C2_T - Ca*k_on_N1_T 0 0 k_off_N1_T 0 0 0 0 0 k2_RT 0 0 0 0 0; 
          0 Ca*k_on_N1_T Ca*k_on_C1_T 0 -k2_TR - k_off_C1_T - k_off_N1_T - Ca*k_on_C2_T - Ca*k_on_N2_T 0 k_off_C2_T k_off_N2_T 0 0 0 0 0 k2_RT 0 0 0 0; 
          0 0 Ca*k_on_N2_T 0 0 -k2_TR - k_off_N2_T - Ca*k_on_C1_T 0 k_off_C1_T 0 0 0 0 0 0 k2_RT 0 0 0; 
          0 0 0 Ca*k_on_N1_T Ca*k_on_C2_T 0 -k3_TR - k_off_C2_T - k_off_N1_T - Ca*k_on_N2_T 0 k_off_N2_T 0 0 0 0 0 0 k3_RT 0 0; 
          0 0 0 0 Ca*k_on_N2_T Ca*k_on_C1_T 0 -k3_TR - k_off_C1_T - k_off_N2_T - Ca*k_on_C2_T k_off_C2_T 0 0 0 0 0 0 0 k3_RT 0; 
          0 0 0 0 0 0 Ca*k_on_N2_T Ca*k_on_C2_T -k4_TR - k_off_C2_T - k_off_N2_T 0 0 0 0 0 0 0 0 k4_RT; 
          k0_TR 0 0 0 0 0 0 0 0 -k0_RT - Ca*k_on_C1_R - Ca*k_on_N1_R k_off_C1_R k_off_N1_R 0 0 0 0 0 0; 
          0 k1_TR 0 0 0 0 0 0 0 Ca*k_on_C1_R -k1_RT - k_off_C1_R - Ca*k_on_C2_R - Ca*k_on_N1_R 0 k_off_C2_R k_off_N1_R 0 0 0 0; 
          0 0 k1_TR 0 0 0 0 0 0 Ca*k_on_N1_R 0 -k1_RT - k_off_N1_R - Ca*k_on_C1_R - Ca*k_on_N2_R 0 k_off_C1_R k_off_N2_R 0 0 0; 
          0 0 0 k2_TR 0 0 0 0 0 0 Ca*k_on_C2_R 0 -k2_RT - k_off_C2_R - Ca*k_on_N1_R 0 0 k_off_N1_R 0 0; 
          0 0 0 0 k2_TR 0 0 0 0 0 Ca*k_on_N1_R Ca*k_on_C1_R 0 -k2_RT - k_off_C1_R - k_off_N1_R - Ca*k_on_C2_R - Ca*k_on_N2_R 0 k_off_C2_R k_off_N2_R 0; 
          0 0 0 0 0 k2_TR 0 0 0 0 0 Ca*k_on_N2_R 0 0 -k2_RT - k_off_N2_R - Ca*k_on_C1_R 0 k_off_C1_R 0; 
          0 0 0 0 0 0 k3_TR 0 0 0 0 0 Ca*k_on_N1_R Ca*k_on_C2_R 0 -k3_RT - k_off_C2_R - k_off_N1_R - Ca*k_on_N2_R 0 k_off_N2_R; 
          0 0 0 0 0 0 0 k3_TR 0 0 0 0 0 Ca*k_on_N2_R Ca*k_on_C1_R 0 -k3_RT - k_off_C1_R - k_off_N2_R - Ca*k_on_C2_R k_off_C2_R; 
###          0 0 0 0 0 0 0 0 k4_TR 0 0 0 0 0 0 Ca*k_on_N2_R Ca*k_on_C2_R -k4_RT - k_off_C2_R - k_off_N2_R;
          1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1]




    b = [zeros(17); CaM_t]
    sol = A\b
    return sol
end



function solve_Byrne_eq(k01, k10, k02, k20, k13, k31, k23, k32, Ca, CaM_T, CaM_equil)

    if CaM_equil

        A = [
            -Ca*k01 - Ca*k02 k10 k20 0; 
             Ca*k01 -k10 - Ca*k13 0 k31; 
             Ca*k02 0 -k20 - Ca*k23 k32; 
             ###0 Ca*k13 Ca*k23 -k31 - k32; 
             1 1 1 1] 

        b = [zeros(3); CaM_T]

        sol = A\b

        return sol
    
    else

        return [CaM_T; zeros(3)]
    end
end



function solve_Bwell_eq(Kd_1, Kd_2, Ca_free, CaM_t, CaM_equil)

    if CaM_equil
        c1      = 1 + Kd_2/Ca_free^2
        c2      = 1 + Ca_free^2/Kd_1
        CaM₀    = (c1 - 1) / (c2 * c1 - 1) * CaM_t
        CaM2Ca₀ = CaM₀ * Ca_free^2 / Kd_1
        CaM4Ca₀ = (c2 - 1) / (c2 * c1 - 1) * CaM_t
        return [CaM₀; CaM2Ca₀; CaM4Ca₀]
    else
        return [CaM_t; 0.0; 0.0]
    end
end