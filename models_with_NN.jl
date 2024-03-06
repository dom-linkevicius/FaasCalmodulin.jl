Faas_params_NN = (;
    tv_kon_TN = 5.9,  
    tv_Kd_TN  = -3.7,
    tv_kon_RN = 7.5, 
    tv_Kd_RN  = -6.1,
    tv_kon_TC = 4.9, 
    tv_Kd_TC  = -4.6,
    tv_kon_RC = 4.4, 
    tv_Kd_RC  = -6.6,
);

Blackwell_params_NN = (;
    tv_kon_1 = 3.78,
    tv_Kd_1  = -6.77,
    tv_kon_2 = 5.0,
    tv_Kd_2  = -4.0,
)

Byrne_params_NN = (;
    tv_k01_N    = log10(275e3),
    tv_K01d_N   = log10(33e-6),
    tv_k02_N    = log10(275e3),
    tv_K02d_N   = log10(229.88e-6),
    tv_k13_N    = log10(507.2e3),
    tv_K13d_N   = log10(3.45e-6),
    tv_k23_N    = log10(500e3),
    tv_K23d_N   = log10(0.5e-6),
    tv_k01_C    = log10(2.75e5),
    tv_K01d_C   = log10(18.5e-6),
    tv_k02_C    = log10(2.75e5),
    tv_K02d_C   = log10(116e-6),
    tv_k13_C    = log10(3.71e3),
    tv_K13d_C   = log10(0.38e-6),
    tv_k23_C    = log10(1.18e5),
    tv_K23d_C   = log10(0.06e-6),
)



Blackwell_scheme_NN_uncaging = @model begin ### Simplest scheme
    @param begin
        tv_kon_1     ∈ Uniform(2,  12)
        tv_Kd_1      ∈ Uniform(-12, -2)
        tv_kon_2     ∈ Uniform(2,  12)
        tv_Kd_2      ∈ Uniform(-12, -2)

        A11           ∈ MvNormal(5, 1.0)
        A12           ∈ MvNormal(5, 1.0)
        b1            ∈ MvNormal(5, 1.0)
        Ae            ∈ MvNormal(5, 1.0)
        be            ∈ Normal(1.0)
        Ω             ∈ RealDomain(; lower=1e-10, init=1, upper=3)
        σ_vec         ∈ VectorDomain(3; lower=1e-10, upper=3, init=1)
    end

    @random begin
        η ~ Normal(Ω)
    end
    
    @covariates begin
        cov_f_frac  
        cov_τ_f  
        cov_τ_s  
        cov_kon_DMn  
        cov_koff_DMn  
        cov_Kd_DMn  
        cov_koff_PP  
        cov_DMn_t  
        cov_CaM_t  
        cov_OGB5_t  
        cov_Ca_free
        PCD
        batch_id
    end



    @pre begin
        σ         = σ_vec[batch_id]

        i1 = tanh.(hcat(A11, A12) * [(PCD - 431)/41.5; η[1]] .+ b1)
        U         = sigmoid(i1' * Ae + be)

        Fmax_Fmin = 39.364

        f_frac   = cov_f_frac 
        τ_f      = 1/cov_τ_f
        τ_s      = 1/cov_τ_s
        kon_DMn  = cov_kon_DMn
        koff_DMn = cov_koff_DMn
        Kd_DMn   = cov_Kd_DMn
        kon_PP   = cov_kon_DMn
        koff_PP  = cov_koff_PP

        DMn_t   = cov_DMn_t
        CaM_t   = cov_CaM_t
        OGB5_t  = cov_OGB5_t
        Ca_free = cov_Ca_free

        kon_D    = 8.77e5         ### in M^-1 ms^-1
        koff_D   = 33.2          ### in ms^-1 
        Kd_D     = koff_D / kon_D
    

        Kd_1  = 10^tv_Kd_1
        kon_1 = 10^tv_kon_1
        koff_1= Kd_1 * kon_1

        Kd_2  = 10^tv_Kd_2 
        kon_2 = 10^tv_kon_2
        koff_2= Kd_2 * kon_2 

        DMn₀     = (Kd_DMn * DMn_t)  / (Ca_free + Kd_DMn)
        CaDMn₀   = DMn_t - DMn₀
        OGB5₀    = (Kd_D   * OGB5_t) / (Ca_free + Kd_D)
        CaOGB5₀  = OGB5_t - OGB5₀

        c1      = 1 + Kd_2/Ca_free^2
        c2      = 1 + Ca_free^2/Kd_1
        CaM₀    = (c1 - 1) / (c2 * c1 - 1) * CaM_t
        CaM2Ca₀ = CaM₀ * Ca_free^2 / Kd_1
        CaM4Ca₀ = (c2 - 1) / (c2 * c1 - 1) * CaM_t

        DMn_s₀_u   = (1-f_frac) * DMn₀   * U# * 0
        CaDMn_s₀_u = (1-f_frac) * CaDMn₀ * U# * 0
        DMn_f₀_u   =    f_frac  * DMn₀   * U# * 0
        CaDMn_f₀_u =    f_frac  * CaDMn₀ * U# * 0
    end

    @init begin
        DMn_s    = DMn_s₀_u 
        CaDMn_s  = CaDMn_s₀_u
        DMn_f    = DMn_f₀_u
        CaDMn_f  = CaDMn_f₀_u
        PP       = 0.0
        CaPP     = 0.0
        OGB5     = OGB5₀
        CaOGB5   = CaOGB5₀
        CaM      = CaM₀
        CaM2Ca   = CaM2Ca₀
        CaM4Ca   = CaM4Ca₀
        Ca       = Ca_free
    end

    @dynamics begin
        DMn_s'    = -kon_DMn * DMn_s * Ca + koff_DMn * CaDMn_s -   DMn_s * τ_s
        CaDMn_s'  =  kon_DMn * DMn_s * Ca - koff_DMn * CaDMn_s - CaDMn_s * τ_s
        DMn_f'    = -kon_DMn * DMn_f * Ca + koff_DMn * CaDMn_f -   DMn_f * τ_f
        CaDMn_f'  =  kon_DMn * DMn_f * Ca - koff_DMn * CaDMn_f - CaDMn_f * τ_f
        PP'       = -kon_PP  * PP    * Ca + koff_PP  * CaPP    +  2DMn_s * τ_s +  2DMn_f * τ_f + CaDMn_s * τ_s + CaDMn_f * τ_f
        CaPP'     =  kon_PP  * PP    * Ca - koff_PP  * CaPP    + CaDMn_s * τ_s + CaDMn_f * τ_f
        OGB5'     = -kon_D   * OGB5  * Ca + koff_D   * CaOGB5
        CaOGB5'   =  kon_D   * OGB5  * Ca - koff_D   * CaOGB5
        CaM'      = -kon_1 * CaM * Ca^2 + koff_1 * CaM2Ca 
        CaM2Ca'   =  kon_1 * CaM * Ca^2 - koff_1 * CaM2Ca -kon_2 * CaM2Ca * Ca^2 + koff_2 * CaM4Ca 
        CaM4Ca'   =                                        kon_2 * CaM2Ca * Ca^2 - koff_2 * CaM4Ca
        Ca'       = -kon_DMn * DMn_s * Ca + koff_DMn * CaDMn_s +
                    -kon_DMn * DMn_f * Ca + koff_DMn * CaDMn_f +
                    -kon_DMn * PP    * Ca + koff_PP  * CaPP    + 
                    -kon_D   * OGB5  * Ca + koff_D   * CaOGB5  +
                    -2kon_1   * CaM   * Ca^2 + 2koff_1 * CaM2Ca  +
                    -2kon_2   * CaM2Ca* Ca^2 + 2koff_2 * CaM4Ca
    end

    @derived begin
        ff0  = @. (OGB5 + Fmax_Fmin * CaOGB5) / (OGB5₀ + Fmax_Fmin * CaOGB5₀)
        F_F0 ~ @. Normal(ff0, σ)
    end

    @observed begin
        DMn_s   = DMn_s
        CaDMn_s = CaDMn_s
        DMn_f   = DMn_f
        CaDMn_f = CaDMn_f
        PP      = PP
        CaPP    = CaPP
        OGB5    = OGB5
        CaOGB5  = CaOGB5
        CaM     = CaM
        CaM2Ca  = CaM2Ca
        CaM4Ca  = CaM4Ca
        Ca      = Ca
    end
end;


Blackwell_scheme_CN = @model begin ### Simplest scheme with unique lobes
    @param begin

        tv_kon_2C     ∈ Uniform(2,  12)
        tv_Kd_2C      ∈ Uniform(-12, -2)
        tv_kon_2N     ∈ Uniform(2,  12)
        tv_Kd_2N      ∈ Uniform(-12, -2)

        tv_mα         ∈ RealDomain(; lower=1e-10, init= 0.0011, upper=1)
        tv_α₀         ∈ RealDomain(;              init = -0.39)

        Ω             ∈ Constrained(Normal(0.1), lower=1e-10) ###, upper=1)
        σ             ∈ RealDomain(; lower=1e-10, init=10)
    end

    @random begin
        η ~ Normal(Ω)
    end
    
    @covariates begin
        cov_f_frac  
        cov_τ_f  
        cov_τ_s  
        cov_kon_DMn  
        cov_koff_DMn  
        cov_Kd_DMn  
        cov_koff_PP  
        cov_DMn_t  
        cov_CaM_t  
        cov_OGB5_t  
        cov_Ca_free
        PCD
    end



    @pre begin
        U         = (tv_mα*PCD + tv_α₀) * (1 + η[1])
        Fmax_Fmin = 39.364

        f_frac   = cov_f_frac 
        τ_f      = 1/cov_τ_f
        τ_s      = 1/cov_τ_s
        kon_DMn  = cov_kon_DMn
        koff_DMn = cov_koff_DMn
        Kd_DMn   = cov_Kd_DMn
        kon_PP   = cov_kon_DMn
        koff_PP  = cov_koff_PP

        DMn_t   = cov_DMn_t
        CaM_t   = cov_CaM_t
        OGB5_t  = cov_OGB5_t
        Ca_free = cov_Ca_free

        kon_D    = 8.77e5         ### in M^-1 ms^-1
        koff_D   = 33.2          ### in ms^-1 
        Kd_D     = koff_D / kon_D

        Kd_2C  = 10^tv_Kd_2C
        kon_2C = 10^tv_kon_2C
        koff_2C= Kd_2C * kon_2C

        Kd_2N  = 10^tv_Kd_2N
        kon_2N = 10^tv_kon_2N
        koff_2N= Kd_2N * kon_2N

        DMn₀     = (Kd_DMn * DMn_t)  / (Ca_free + Kd_DMn)
        CaDMn₀   = DMn_t - DMn₀
        OGB5₀    = (Kd_D   * OGB5_t) / (Ca_free + Kd_D)
        CaOGB5₀  = OGB5_t - OGB5₀

        c0       = solve_Blackwell_CN_eqs(kon_2C, koff_2C, kon_2N, koff_2N, Ca_free, CaM_t)

        DMn_s₀_u   = (1-f_frac) * DMn₀   * U# * 0
        CaDMn_s₀_u = (1-f_frac) * CaDMn₀ * U# * 0
        DMn_f₀_u   =    f_frac  * DMn₀   * U# * 0
        CaDMn_f₀_u =    f_frac  * CaDMn₀ * U# * 0
    end

    @init begin
        DMn_s    = DMn_s₀_u 
        CaDMn_s  = CaDMn_s₀_u
        DMn_f    = DMn_f₀_u
        CaDMn_f  = CaDMn_f₀_u
        PP       = 0.0
        CaPP     = 0.0
        OGB5     = OGB5₀
        CaOGB5   = CaOGB5₀
        CaM0     = c0[1]
        CaM2C    = c0[2]
        CaM2N    = c0[3]
        CaM4     = c0[4]
        Ca       = Ca_free
    end

    @dynamics begin
        DMn_s'    = -kon_DMn * DMn_s * Ca + koff_DMn * CaDMn_s -   DMn_s * τ_s
        CaDMn_s'  =  kon_DMn * DMn_s * Ca - koff_DMn * CaDMn_s - CaDMn_s * τ_s
        DMn_f'    = -kon_DMn * DMn_f * Ca + koff_DMn * CaDMn_f -   DMn_f * τ_f
        CaDMn_f'  =  kon_DMn * DMn_f * Ca - koff_DMn * CaDMn_f - CaDMn_f * τ_f
        PP'       = -kon_PP  * PP    * Ca + koff_PP  * CaPP    +  2DMn_s * τ_s +  2DMn_f * τ_f + CaDMn_s * τ_s + CaDMn_f * τ_f
        CaPP'     =  kon_PP  * PP    * Ca - koff_PP  * CaPP    + CaDMn_s * τ_s + CaDMn_f * τ_f
        OGB5'     = -kon_D   * OGB5  * Ca + koff_D   * CaOGB5
        CaOGB5'   =  kon_D   * OGB5  * Ca - koff_D   * CaOGB5
        CaM0'     = -kon_2C * CaM0  * Ca^2 + koff_2C * CaM2C +
                    -kon_2N * CaM0  * Ca^2 + koff_2N * CaM2N 
        CaM2C'    =  kon_2C * CaM0  * Ca^2 - koff_2C * CaM2C +
                    -kon_2N * CaM2C * Ca^2 + koff_2N * CaM4 
        CaM2N'    =  kon_2N * CaM0  * Ca^2 - koff_2N * CaM2N +
                    -kon_2C * CaM2N * Ca^2 + koff_2C * CaM4 
        CaM4'     =  kon_2N * CaM2C * Ca^2 - koff_2N * CaM4 +
                     kon_2C * CaM2N * Ca^2 - koff_2C * CaM4 
        Ca'       = -kon_DMn * DMn_s * Ca + koff_DMn * CaDMn_s +
                    -kon_DMn * DMn_f * Ca + koff_DMn * CaDMn_f +
                    -kon_DMn * PP    * Ca + koff_PP  * CaPP    + 
                    -kon_D   * OGB5  * Ca + koff_D   * CaOGB5  +
                    -kon_2C * CaM0  * Ca^2 + koff_2C * CaM2C +
                    -kon_2N * CaM0  * Ca^2 + koff_2N * CaM2N +
                    -kon_2N * CaM2C * Ca^2 + koff_2N * CaM4 +
                    -kon_2C * CaM2N * Ca^2 + koff_2C * CaM4 
    end

    @derived begin
        F_F0 = @. Normal((OGB5 + Fmax_Fmin * CaOGB5) / (OGB5₀ + Fmax_Fmin * CaOGB5₀), σ)
    end

    @observed begin
        DMn_s   = DMn_s
        CaDMn_s = CaDMn_s
        DMn_f   = DMn_f
        CaDMn_f = CaDMn_f
        PP      = PP
        CaPP    = CaPP
        OGB5    = OGB5
        CaOGB5  = CaOGB5
        CaM0    = CaM0
        CaM2C   = CaM2C
        CaM2N   = CaM2N
        CaM4    = CaM4
        Ca      = Ca
    end
end;


Blackwell_scheme_TR = @model begin ### Simplest scheme with tense/relaxed
    @param begin

        tv_kon_1_T     ∈ Uniform(2,   12)
        tv_Kd_1_T      ∈ Uniform(-12, -2)
        tv_kon_2_T     ∈ Uniform(2,   12)
        tv_Kd_2_T      ∈ Uniform(-12, -2)

        tv_kon_1_R     ∈ Uniform(2,   12)
        tv_Kd_1_R      ∈ Uniform(-12, -2)
        tv_kon_2_R     ∈ Uniform(2,   12)
        tv_Kd_2_R      ∈ Uniform(-12, -2)

        tv_k_TR_0      ∈ Uniform(2,   12)
        tv_Kd_0        ∈ Uniform(-12, -2)
        tv_k_TR_2      ∈ Uniform(2,   12)
        tv_Kd_2        ∈ Uniform(-12, -2)
        tv_k_TR_4      ∈ Uniform(2,   12)
        tv_Kd_4        ∈ Uniform(-12, -2)

        tv_mα         ∈ RealDomain(; lower=1e-10, init= 0.0011, upper=1)
        tv_α₀         ∈ RealDomain(;              init = -0.39)

        Ω             ∈ Constrained(Normal(0.1), lower=1e-10) ###, upper=1)
        σ             ∈ RealDomain(; lower=1e-10, init=10)
    end

    @random begin
        η ~ Normal(Ω)
    end
    
    @covariates begin
        cov_f_frac  
        cov_τ_f  
        cov_τ_s  
        cov_kon_DMn  
        cov_koff_DMn  
        cov_Kd_DMn  
        cov_koff_PP  
        cov_DMn_t  
        cov_CaM_t  
        cov_OGB5_t  
        cov_Ca_free
        PCD
    end



    @pre begin
        U         = (tv_mα*PCD + tv_α₀) * (1 + η[1])
        Fmax_Fmin = 39.364

        f_frac   = cov_f_frac 
        τ_f      = 1/cov_τ_f
        τ_s      = 1/cov_τ_s
        kon_DMn  = cov_kon_DMn
        koff_DMn = cov_koff_DMn
        Kd_DMn   = cov_Kd_DMn
        kon_PP   = cov_kon_DMn
        koff_PP  = cov_koff_PP

        DMn_t   = cov_DMn_t
        CaM_t   = cov_CaM_t
        OGB5_t  = cov_OGB5_t
        Ca_free = cov_Ca_free

        kon_D    = 8.77e5         ### in M^-1 ms^-1
        koff_D   = 33.2          ### in ms^-1 
        Kd_D     = koff_D / kon_D

        Kd_1_T  = 10^tv_Kd_1_T
        kon_1_T = 10^tv_kon_1_T
        koff_1_T= Kd_1_T * kon_1_T

        Kd_2_T  = 10^tv_Kd_2_T
        kon_2_T = 10^tv_kon_2_T
        koff_2_T= Kd_2_T * kon_2_T

        Kd_1_R  = 10^tv_Kd_1_R
        kon_1_R = 10^tv_kon_1_R
        koff_1_R= Kd_1_R * kon_1_R

        Kd_2_R  = 10^tv_Kd_2_R
        kon_2_R = 10^tv_kon_2_R
        koff_2_R= Kd_2_R * kon_2_R

        Kd_0   = 10^tv_Kd_0
        k_TR_0 = 10^tv_k_TR_0
        k_RT_0 = Kd_0 * k_TR_0

        Kd_2   = 10^tv_Kd_2
        k_TR_2 = 10^tv_k_TR_2
        k_RT_2 = Kd_2 * k_TR_2

        Kd_4   = 10^tv_Kd_4
        k_TR_4 = 10^tv_k_TR_4
        k_RT_4 = Kd_4 * k_TR_4

        DMn₀     = (Kd_DMn * DMn_t)  / (Ca_free + Kd_DMn)
        CaDMn₀   = DMn_t - DMn₀
        OGB5₀    = (Kd_D   * OGB5_t) / (Ca_free + Kd_D)
        CaOGB5₀  = OGB5_t - OGB5₀

        c0       = solve_Blackwell_TR_eqs(kon_1_T, koff_1_T, kon_2_T, koff_2_T, 
                                          kon_1_R, koff_1_R, kon_2_R, koff_2_R, 
                                          k_TR_0, k_RT_0, k_TR_2, k_RT_2, k_TR_4, k_RT_4,
                                          Ca_free, CaM_t)

        DMn_s₀_u   = (1-f_frac) * DMn₀   * U# * 0
        CaDMn_s₀_u = (1-f_frac) * CaDMn₀ * U# * 0
        DMn_f₀_u   =    f_frac  * DMn₀   * U# * 0
        CaDMn_f₀_u =    f_frac  * CaDMn₀ * U# * 0
    end

    @init begin
        DMn_s    = DMn_s₀_u 
        CaDMn_s  = CaDMn_s₀_u
        DMn_f    = DMn_f₀_u
        CaDMn_f  = CaDMn_f₀_u
        PP       = 0.0
        CaPP     = 0.0
        OGB5     = OGB5₀
        CaOGB5   = CaOGB5₀
        CaM0_R   = c0[1]
        CaM2_R   = c0[2]
        CaM4_R   = c0[3]
        CaM0_T   = c0[4]
        CaM2_T   = c0[5]
        CaM4_T   = c0[6]
        Ca       = Ca_free
    end

    @dynamics begin
        DMn_s'    = -kon_DMn * DMn_s * Ca + koff_DMn * CaDMn_s -   DMn_s * τ_s
        CaDMn_s'  =  kon_DMn * DMn_s * Ca - koff_DMn * CaDMn_s - CaDMn_s * τ_s
        DMn_f'    = -kon_DMn * DMn_f * Ca + koff_DMn * CaDMn_f -   DMn_f * τ_f
        CaDMn_f'  =  kon_DMn * DMn_f * Ca - koff_DMn * CaDMn_f - CaDMn_f * τ_f
        PP'       = -kon_PP  * PP    * Ca + koff_PP  * CaPP    +  2DMn_s * τ_s +  2DMn_f * τ_f + CaDMn_s * τ_s + CaDMn_f * τ_f
        CaPP'     =  kon_PP  * PP    * Ca - koff_PP  * CaPP    + CaDMn_s * τ_s + CaDMn_f * τ_f
        OGB5'     = -kon_D   * OGB5  * Ca + koff_D   * CaOGB5
        CaOGB5'   =  kon_D   * OGB5  * Ca - koff_D   * CaOGB5
        CaM0_R'   = -kon_1_R * CaM0_R * Ca^2 + koff_1_R * CaM2_R + 
                    -k_RT_0  * CaM0_R + k_TR_0 * CaM0_T 
        CaM2_R'   =  kon_1_R * CaM0_R * Ca^2 - koff_1_R * CaM2_R +
                    -kon_2_R * CaM2_R * Ca^2 + koff_2_R * CaM4_R +
                    -k_RT_2  * CaM2_R + k_TR_2 * CaM2_T 
        CaM4_R'   =  kon_2_R * CaM2_R * Ca^2 - koff_2_R * CaM4_R +
                    -k_RT_4  * CaM4_R + k_TR_4 * CaM4_T 
        CaM0_T'   = -kon_1_T * CaM0_T * Ca^2 + koff_1_T * CaM2_T +
                     k_RT_0  * CaM0_R - k_TR_0 * CaM0_T 
        CaM2_T'   =  kon_1_T * CaM0_T * Ca^2 - koff_1_T * CaM2_T +
                    -kon_2_T * CaM2_T * Ca^2 + koff_2_T * CaM4_T +
                     k_RT_2  * CaM2_R - k_TR_2 * CaM2_T 
        CaM4_T'   =  kon_2_T * CaM2_T * Ca^2 - koff_2_T * CaM4_T +
                     k_RT_4  * CaM4_R - k_TR_4 * CaM4_T
        Ca'       = -kon_DMn * DMn_s * Ca + koff_DMn * CaDMn_s +
                    -kon_DMn * DMn_f * Ca + koff_DMn * CaDMn_f +
                    -kon_DMn * PP    * Ca + koff_PP  * CaPP    + 
                    -kon_D   * OGB5  * Ca + koff_D   * CaOGB5  +
                    -kon_1_R * CaM0_R * Ca^2 + koff_1_R * CaM2_R + 
                    -kon_2_R * CaM2_R * Ca^2 + koff_2_R * CaM4_R +
                    -kon_1_T * CaM0_T * Ca^2 + koff_1_T * CaM2_T +
                    -kon_2_T * CaM2_T * Ca^2 + koff_2_T * CaM4_T
    end

    @derived begin
        F_F0 = @. Normal((OGB5 + Fmax_Fmin * CaOGB5) / (OGB5₀ + Fmax_Fmin * CaOGB5₀), σ)
    end

    @observed begin
        DMn_s   = DMn_s
        CaDMn_s = CaDMn_s
        DMn_f   = DMn_f
        CaDMn_f = CaDMn_f
        PP      = PP
        CaPP    = CaPP
        OGB5    = OGB5
        CaOGB5  = CaOGB5
        CaM0_R  = CaM0_R
        CaM2_R  = CaM2_R
        CaM4_R  = CaM4_R
        CaM0_T  = CaM0_T
        CaM2_T  = CaM2_T
        CaM4_T  = CaM4_T
        Ca      = Ca
    end
end;


Faas_scheme_NN_uncaging = @model begin ### Faas scheme
    @param begin
        tv_kon_TN     ∈ Uniform(2, 12)
        tv_Kd_TN      ∈ Uniform(-12, -2) 
        tv_kon_RN     ∈ Uniform(2, 12)
        tv_Kd_RN      ∈ Uniform(-12, -2)

        tv_kon_TC     ∈ Uniform(2, 12)
        tv_Kd_TC      ∈ Uniform(-12, -2)
        tv_kon_RC     ∈ Uniform(2, 12)
        tv_Kd_RC      ∈ Uniform(-12, -2)

        A11           ∈ MvNormal(5, 1.0)
        A12           ∈ MvNormal(5, 1.0)
        b1            ∈ MvNormal(5, 1.0)
        Ae            ∈ MvNormal(5, 1.0)
        be            ∈ Normal(1.0)
        Ω             ∈ RealDomain(; lower=1e-10, init=1, upper=3)
        σ_vec         ∈ VectorDomain(3; lower=1e-10, upper=3, init=1)
    end

    @random begin
        η     ~ Normal(Ω)
    end
    
    @covariates begin
        cov_f_frac  
        cov_τ_f  
        cov_τ_s  
        cov_kon_DMn  
        cov_koff_DMn  
        cov_Kd_DMn  
        cov_koff_PP  
        cov_DMn_t  
        cov_CaM_t  
        cov_OGB5_t  
        cov_Ca_free
        PCD
        batch_id
    end



    @pre begin
        σ         = σ_vec[batch_id]

        i1 = tanh.(hcat(A11, A12) * [(PCD - 431)/41.5; η[1]] .+ b1)
        U         = sigmoid(i1' * Ae + be)

        Fmax_Fmin = 39.364

        f_frac   = cov_f_frac   
        τ_f      = 1/cov_τ_f    
        τ_s      = 1/cov_τ_s    

        kon_DMn  = cov_kon_DMn  
        Kd_DMn   = cov_Kd_DMn   
        koff_DMn = Kd_DMn * kon_DMn

        kon_PP   = cov_kon_DMn  
        koff_PP  = cov_koff_PP  

        DMn_t   = cov_DMn_t
        CaM_t   = cov_CaM_t
        OGB5_t  = cov_OGB5_t
        Ca_free = cov_Ca_free

        kon_D    = 8.77e5         ### in M^-1 ms^-1
        koff_D   = 33.2          ### in ms^-1 
        Kd_D     = koff_D / kon_D
    

        Kd_TN    = 10^tv_Kd_TN 
        kon_TN   = 10^tv_kon_TN
        koff_TN  = Kd_TN * kon_TN 

        Kd_RN    = 10^tv_Kd_RN 
        kon_RN   = 10^tv_kon_RN
        koff_RN  = Kd_RN * kon_RN 

        Kd_TC    = 10^tv_Kd_TC 
        kon_TC   = 10^tv_kon_TC
        koff_TC  = Kd_TC * kon_TC 

        Kd_RC    = 10^tv_Kd_RC 
        kon_RC   = 10^tv_kon_RC
        koff_RC  = Kd_RC * kon_RC 


        DMn₀     = (Kd_DMn * DMn_t)  / (Ca_free + Kd_DMn)
        CaDMn₀   = DMn_t - DMn₀
        OGB5₀    = (Kd_D   * OGB5_t) / (Ca_free + Kd_D)
        CaOGB5₀  = OGB5_t - OGB5₀

        c1_N  = (2 * Ca_free) / Kd_TN
        c2_N  = (2 * Kd_RN) / Ca_free
        NtNt₀ = -(c2_N * CaM_t) / (1 - (c2_N+1)*(c1_N + 1))
        NtNr₀ = c1_N * NtNt₀
        NrNr₀ = NtNr₀ / c2_N

        c1_C  = (2 * Ca_free) / Kd_TC
        c2_C  = (2 * Kd_RC) / Ca_free
        CtCt₀ = -(c2_C * CaM_t) / (1 - (c2_C+1)*(c1_C + 1))
        CtCr₀ = c1_C * CtCt₀
        CrCr₀ = CtCr₀ / c2_C

        DMn_s₀_u   = (1-f_frac) * DMn₀   * U
        CaDMn_s₀_u = (1-f_frac) * CaDMn₀ * U
        DMn_f₀_u   =    f_frac  * DMn₀   * U
        CaDMn_f₀_u =    f_frac  * CaDMn₀ * U
    end

    @init begin
        DMn_s    = DMn_s₀_u 
        CaDMn_s  = CaDMn_s₀_u
        DMn_f    = DMn_f₀_u
        CaDMn_f  = CaDMn_f₀_u
        PP       = 0.0
        CaPP     = 0.0
        OGB5     = OGB5₀
        CaOGB5   = CaOGB5₀
        NtNt     = NtNt₀
        NtNr     = NtNr₀
        NrNr     = NrNr₀
        CtCt     = CtCt₀
        CtCr     = CtCr₀
        CrCr     = CrCr₀
        Ca       = Ca_free
    end

    @dynamics begin
        DMn_s'    = -kon_DMn * DMn_s * Ca + koff_DMn * CaDMn_s -   DMn_s * τ_s
        CaDMn_s'  =  kon_DMn * DMn_s * Ca - koff_DMn * CaDMn_s - CaDMn_s * τ_s
        DMn_f'    = -kon_DMn * DMn_f * Ca + koff_DMn * CaDMn_f -   DMn_f * τ_f
        CaDMn_f'  =  kon_DMn * DMn_f * Ca - koff_DMn * CaDMn_f - CaDMn_f * τ_f
        PP'       = -kon_PP  * PP    * Ca + koff_PP  * CaPP    +  2DMn_s * τ_s +  2DMn_f * τ_f + CaDMn_s * τ_s + CaDMn_f * τ_f
        CaPP'     =  kon_PP  * PP    * Ca - koff_PP  * CaPP    + CaDMn_s * τ_s + CaDMn_f * τ_f
        OGB5'     = -kon_D   * OGB5  * Ca + koff_D   * CaOGB5
        CaOGB5'   =  kon_D   * OGB5  * Ca - koff_D   * CaOGB5
        NtNt'     = -2kon_TN * NtNt  * Ca + koff_TN  * NtNr
        NtNr'     =  2kon_TN * NtNt  * Ca - koff_TN  * NtNr - kon_RN * NtNr * Ca + 2koff_RN  * NrNr
        NrNr'     =                                         + kon_RN * NtNr * Ca - 2koff_RN  * NrNr
        CtCt'     = -2kon_TC * CtCt  * Ca + koff_TC  * CtCr
        CtCr'     =  2kon_TC * CtCt  * Ca - koff_TC  * CtCr - kon_RC * CtCr * Ca + 2koff_RC  * CrCr
        CrCr'     =                                         + kon_RC * CtCr * Ca - 2koff_RC  * CrCr
        Ca'       = -kon_DMn * DMn_s * Ca + koff_DMn * CaDMn_s +
                    -kon_DMn * DMn_f * Ca + koff_DMn * CaDMn_f +
                    -kon_DMn * PP    * Ca + koff_PP  * CaPP    + 
                    -kon_D   * OGB5  * Ca + koff_D   * CaOGB5  +
                    -2kon_TN * NtNt  * Ca + koff_TN  * NtNr    +
                    -2kon_TC * CtCt  * Ca + koff_TC  * CtCr    +     
                    -kon_RN  * NtNr  * Ca + 2koff_RN * NrNr    +
                    -kon_RC  * CtCr  * Ca + 2koff_RC * CrCr
    end

    @derived begin
        signal = @. (OGB5 + Fmax_Fmin * CaOGB5) / (OGB5₀ + Fmax_Fmin * CaOGB5₀)
        F_F0 = @. Normal(signal, σ)
    end

    @observed begin
        DMn_s   = DMn_s
        CaDMn_s = CaDMn_s
        DMn_f   = DMn_f
        CaDMn_f = CaDMn_f
        PP      = PP
        CaPP    = CaPP
        OGB5    = OGB5
        CaOGB5  = CaOGB5
        NtNt    = NtNt
        NtNr    = NtNr
        NrNr    = NrNr
        CtCt    = CtCt
        CtCr    = CtCr
        CrCr    = CrCr
        Ca      = Ca
    end

end;



Byrne_scheme = @model begin ### Byrne scheme
    @param begin
        tv_k01_N      ∈ Uniform(2,   12)
        tv_K01d_N     ∈ Uniform(-12, -2)
        tv_k02_N      ∈ Uniform(2,   12)
        tv_K02d_N     ∈ Uniform(-12, -2)
        tv_k13_N      ∈ Uniform(2,   12)
        tv_K13d_N     ∈ Uniform(-12, -2)
        tv_k23_N      ∈ Uniform(2,   12)
        tv_K23d_N     ∈ Uniform(-12, -2)

        tv_k01_C      ∈ Uniform(2,   12)
        tv_K01d_C     ∈ Uniform(-12, -2)
        tv_k02_C      ∈ Uniform(2,   12)
        tv_K02d_C     ∈ Uniform(-12, -2)
        tv_k13_C      ∈ Uniform(2,   12)
        tv_K13d_C     ∈ Uniform(-12, -2)
        tv_k23_C      ∈ Uniform(2,   12)
        tv_K23d_C     ∈ Uniform(-12, -2)

        tv_mα         ∈ RealDomain(; lower=1e-10, init= 0.001, upper=1)
        tv_α₀         ∈ RealDomain(;              init = -0.39)

        Ω             ∈ Constrained(Normal(0.1), lower=1e-10) ###, upper=1)
        σ             ∈ RealDomain(; lower=1e-10, init=10)
    end

    @random begin
        η ~ Normal(Ω)
    end
    
    @covariates begin
        cov_f_frac  
        cov_τ_f  
        cov_τ_s  
        cov_kon_DMn  
        cov_koff_DMn  
        cov_Kd_DMn  
        cov_koff_PP  
        cov_DMn_t  
        cov_CaM_t  
        cov_OGB5_t  
        cov_Ca_free
        PCD
    end



    @pre begin
        U         = (tv_mα*PCD + tv_α₀) * (1 + η[1])
        Fmax_Fmin = 39.364

        f_frac   = cov_f_frac 
        τ_f      = 1/cov_τ_f
        τ_s      = 1/cov_τ_s
        kon_DMn  = cov_kon_DMn
        koff_DMn = cov_koff_DMn
        Kd_DMn   = cov_Kd_DMn
        kon_PP   = cov_kon_DMn
        koff_PP  = cov_koff_PP

        DMn_t   = cov_DMn_t
        CaM_t   = cov_CaM_t
        OGB5_t  = cov_OGB5_t
        Ca_free = cov_Ca_free

        kon_D    = 8.77e5         ### in M^-1 ms^-1
        koff_D   = 33.2          ### in ms^-1 
        Kd_D     = koff_D / kon_D

        K01d_N  = 10^tv_K01d_N 
        k01_N   = 10^tv_k01_N
        k10_N   = K01d_N * k01_N 

        K02d_N  = 10^tv_K02d_N 
        k02_N   = 10^tv_k02_N
        k20_N   = K02d_N * k02_N 

        K13d_N  = 10^tv_K13d_N 
        k13_N   = 10^tv_k13_N
        k31_N   = K13d_N * k13_N 

        K23d_N  = 10^tv_K23d_N 
        k23_N   = 10^tv_k23_N
        k32_N   = K23d_N * k23_N 

        K01d_C  = 10^tv_K01d_C 
        k01_C   = 10^tv_k01_C
        k10_C   = K01d_C * k01_C 

        K02d_C  = 10^tv_K02d_C 
        k02_C   = 10^tv_k02_C
        k20_C   = K02d_C * k02_C 

        K13d_C  = 10^tv_K13d_C 
        k13_C   = 10^tv_k13_C
        k31_C   = K13d_C * k13_C 

        K23d_C  = 10^tv_K23d_C 
        k23_C   = 10^tv_k23_C
        k32_C   = K23d_C * k23_C 


        DMn₀     = (Kd_DMn * DMn_t)  / (Ca_free + Kd_DMn)
        CaDMn₀   = DMn_t - DMn₀
        OGB5₀    = (Kd_D   * OGB5_t) / (Ca_free + Kd_D)
        CaOGB5₀  = OGB5_t - OGB5₀


        N₀ = solve_Byrne_eq(k01_N, k10_N, k02_N, k20_N, k13_N, k31_N, k23_N, k32_N, Ca_free, CaM_t)
        C₀ = solve_Byrne_eq(k01_C, k10_C, k02_C, k20_C, k13_C, k31_C, k23_C, k32_C, Ca_free, CaM_t)
        

        DMn_s₀_u   = (1-f_frac) * DMn₀   * U * 1
        CaDMn_s₀_u = (1-f_frac) * CaDMn₀ * U * 1
        DMn_f₀_u   =    f_frac  * DMn₀   * U * 1
        CaDMn_f₀_u =    f_frac  * CaDMn₀ * U * 1
    end

    @init begin
        DMn_s    = DMn_s₀_u 
        CaDMn_s  = CaDMn_s₀_u
        DMn_f    = DMn_f₀_u
        CaDMn_f  = CaDMn_f₀_u
        PP       = 0.0
        CaPP     = 0.0
        OGB5     = OGB5₀
        CaOGB5   = CaOGB5₀
        N0       = N₀[1]
        N1       = N₀[2]
        N2       = N₀[3]
        N3       = N₀[4]
        C0       = C₀[1]
        C1       = C₀[2]
        C2       = C₀[3]
        C3       = C₀[4]
        Ca       = Ca_free
    end

    @dynamics begin
        DMn_s'    = -kon_DMn * DMn_s * Ca + koff_DMn * CaDMn_s -   DMn_s * τ_s
        CaDMn_s'  =  kon_DMn * DMn_s * Ca - koff_DMn * CaDMn_s - CaDMn_s * τ_s
        DMn_f'    = -kon_DMn * DMn_f * Ca + koff_DMn * CaDMn_f -   DMn_f * τ_f
        CaDMn_f'  =  kon_DMn * DMn_f * Ca - koff_DMn * CaDMn_f - CaDMn_f * τ_f
        PP'       = -kon_PP  * PP    * Ca + koff_PP  * CaPP    +  2DMn_s * τ_s +  2DMn_f * τ_f + CaDMn_s * τ_s + CaDMn_f * τ_f
        CaPP'     =  kon_PP  * PP    * Ca - koff_PP  * CaPP    + CaDMn_s * τ_s + CaDMn_f * τ_f
        OGB5'     = -kon_D   * OGB5  * Ca + koff_D   * CaOGB5
        CaOGB5'   =  kon_D   * OGB5  * Ca - koff_D   * CaOGB5
        N0'       = -k01_N * Ca * N0 + k10_N * N1 +
                    -k02_N * Ca * N0 + k20_N * N2
        N1'       =  k01_N * Ca * N0 - k10_N * N1 +
                    -k13_N * Ca * N1 + k31_N * N3  
        N2'       =  k02_N * Ca * N0 - k20_N * N2 +
                    -k23_N * Ca * N2 + k32_N * N3
        N3'       =  k13_N * Ca * N1 - k31_N * N3 +
                     k23_N * Ca * N2 - k32_N * N3
        C0'       = -k01_C * Ca * C0 + k10_C * C1 +
                    -k02_C * Ca * C0 + k20_C * C2
        C1'       =  k01_C * Ca * C0 - k10_C * C1 +
                    -k13_C * Ca * C1 + k31_C * C3  
        C2'       =  k02_C * Ca * C0 - k20_C * C2 +
                    -k23_C * Ca * C2 + k32_C * C3
        C3'       =  k13_C * Ca * C1 - k31_C * C3 +
                     k23_C * Ca * C2 - k32_C * C3
        Ca'       = -kon_DMn * DMn_s * Ca + koff_DMn * CaDMn_s +
                    -kon_DMn * DMn_f * Ca + koff_DMn * CaDMn_f +
                    -kon_DMn * PP    * Ca + koff_PP  * CaPP    + 
                    -kon_D   * OGB5  * Ca + koff_D   * CaOGB5  +
                    -k01_N * Ca * N0 + k10_N * N1 +
                    -k02_N * Ca * N0 + k20_N * N2 +
                    -k13_N * Ca * N1 + k31_N * N3 +
                    -k23_N * Ca * N2 + k32_N * N3 +
                    -k01_C * Ca * C0 + k10_C * C1 +
                    -k02_C * Ca * C0 + k20_C * C2 +
                    -k13_C * Ca * C1 + k31_C * C3 + 
                    -k23_C * Ca * C2 + k32_C * C3
    end

    @derived begin
        F_F0 = @. Normal((OGB5 + Fmax_Fmin * CaOGB5) / (OGB5₀ + Fmax_Fmin * CaOGB5₀), σ)
    end

    @observed begin
        DMn_s   = DMn_s
        CaDMn_s = CaDMn_s
        DMn_f   = DMn_f
        CaDMn_f = CaDMn_f
        PP      = PP
        CaPP    = CaPP
        OGB5    = OGB5
        CaOGB5  = CaOGB5
        N0      = N0
        N1      = N1
        N2      = N2
        N3      = N3
        C0      = C0
        C1      = C1
        C2      = C2
        C3      = C3
        Ca      = Ca
    end

end;



###Pepke_TR_model = @model begin
###
###    @param begin
###        tv_k_on_C1_T  ∈ Uniform(2, 12)
###        tv_Kd_C1_T    ∈ Uniform(-8, -2)
###        tv_k_on_C2_T  ∈ Uniform(2, 12)
###        tv_Kd_C2_T    ∈ Uniform(-8, -2)
###        tv_k_on_N1_T  ∈ Uniform(2, 12)
###        tv_Kd_N1_T    ∈ Uniform(-8, -2)
###        tv_k_on_N2_T  ∈ Uniform(2, 12)
###        tv_Kd_N2_T    ∈ Uniform(-8, -2)
###
###        tv_k_on_C1_R  ∈ Uniform(2, 12)
###        tv_Kd_C1_R    ∈ Uniform(-8, -2)
###        tv_k_on_C2_R  ∈ Uniform(2, 12)
###        tv_Kd_C2_R    ∈ Uniform(-8, -2)
###        tv_k_on_N1_R  ∈ Uniform(2, 12)
###        tv_Kd_N1_R    ∈ Uniform(-8, -2)
###        tv_k_on_N2_R  ∈ Uniform(2, 12)
###        tv_Kd_N2_R    ∈ Uniform(-8, -2)
###
###
###        tv_k0_TR      ∈ Uniform(2, 12)
###        tv_Kd_0_TR    ∈ Uniform(-8, -2)
###        tv_k1_TR      ∈ Uniform(2, 12)
###        tv_Kd_1_TR    ∈ Uniform(-8, -2)
###        tv_k2_TR      ∈ Uniform(2, 12)
###        tv_Kd_2_TR    ∈ Uniform(-8, -2)
###        tv_k3_TR      ∈ Uniform(2, 12)
###        tv_Kd_3_TR    ∈ Uniform(-8, -2)
###        tv_k4_TR      ∈ Uniform(2, 12)
###        tv_Kd_4_TR    ∈ Uniform(-8, -2)
###
###        tv_mα         ∈ RealDomain(; lower=1e-10, init= 0.0011, upper=1)
###        tv_α₀         ∈ RealDomain(;              init = -0.39)
###
###        Ω             ∈ Constrained(Normal(0.1), lower=1e-10) ###, upper=1)
###        σ             ∈ RealDomain(; lower=1e-10, init=10)
###    end
###
###    @random begin
###        η ~ Normal(Ω)
###    end
###    
###    @covariates begin
###        cov_f_frac  
###        cov_τ_f  
###        cov_τ_s  
###        cov_kon_DMn  
###        cov_koff_DMn  
###        cov_Kd_DMn  
###        cov_koff_PP  
###        cov_DMn_t  
###        cov_CaM_t  
###        cov_OGB5_t  
###        cov_Ca_free
###        PCD
###    end
###
###    @pre begin
###        U         = (tv_mα*PCD + tv_α₀) * (1 + η[1])
###        Fmax_Fmin = 39.364
###
###        f_frac   = cov_f_frac 
###        τ_f      = 1/cov_τ_f
###        τ_s      = 1/cov_τ_s
###        kon_DMn  = cov_kon_DMn
###        koff_DMn = cov_koff_DMn
###        Kd_DMn   = cov_Kd_DMn
###        kon_PP   = cov_kon_DMn
###        koff_PP  = cov_koff_PP
###
###        DMn_t   = cov_DMn_t
###        CaM_t   = cov_CaM_t
###        OGB5_t  = cov_OGB5_t
###        Ca_free = cov_Ca_free
###
###        kon_D    = 8.77e5         ### in M^-1 ms^-1
###        koff_D   = 33.2          ### in ms^-1 
###        Kd_D     = koff_D / kon_D
###
###        Kd_C1_T    = 10^tv_Kd_C1_T 
###        k_on_C1_T  = 10^tv_k_on_C1_T
###        k_off_C1_T = Kd_C1_T * k_on_C1_T 
###
###        Kd_C2_T    = 10^tv_Kd_C2_T 
###        k_on_C2_T  = 10^tv_k_on_C2_T
###        k_off_C2_T = Kd_C2_T * k_on_C2_T
###
###        Kd_N1_T    = 10^tv_Kd_N1_T 
###        k_on_N1_T  = 10^tv_k_on_N1_T
###        k_off_N1_T = Kd_N1_T * k_on_N1_T
###
###        Kd_N2_T    = 10^tv_Kd_N2_T 
###        k_on_N2_T  = 10^tv_k_on_N2_T
###        k_off_N2_T = Kd_N2_T * k_on_N2_T
###
###        Kd_C1_R    = 10^tv_Kd_C1_R 
###        k_on_C1_R  = 10^tv_k_on_C1_R
###        k_off_C1_R = Kd_C1_R * k_on_C1_R 
###
###        Kd_C2_R    = 10^tv_Kd_C2_R 
###        k_on_C2_R  = 10^tv_k_on_C2_R
###        k_off_C2_R = Kd_C2_R * k_on_C2_R
###
###        Kd_N1_R    = 10^tv_Kd_N1_R 
###        k_on_N1_R  = 10^tv_k_on_N1_R
###        k_off_N1_R = Kd_N1_R * k_on_N1_R
###
###        Kd_N2_R    = 10^tv_Kd_N2_R 
###        k_on_N2_R  = 10^tv_k_on_N2_R
###        k_off_N2_R = Kd_N2_R * k_on_N2_R
###
###        Kd_0_TR = 10^tv_Kd_0_TR 
###        k0_TR   = 10^tv_k0_TR
###        k0_RT   = Kd_0_TR * k0_TR 
###
###        Kd_1_TR = 10^tv_Kd_1_TR 
###        k1_TR   = 10^tv_k1_TR
###        k1_RT   = Kd_1_TR * k1_TR
###
###        Kd_2_TR = 10^tv_Kd_2_TR 
###        k2_TR   = 10^tv_k2_TR
###        k2_RT   = Kd_2_TR * k2_TR
###
###        Kd_3_TR = 10^tv_Kd_3_TR 
###        k3_TR   = 10^tv_k3_TR
###        k3_RT   = Kd_3_TR * k3_TR
###
###        Kd_4_TR = 10^tv_Kd_4_TR 
###        k4_TR   = 10^tv_k4_TR
###        k4_RT   = Kd_4_TR * k4_TR
###
###        DMn₀     = (Kd_DMn * DMn_t)  / (Ca_free + Kd_DMn)
###        CaDMn₀   = DMn_t - DMn₀
###        OGB5₀    = (Kd_D   * OGB5_t) / (Ca_free + Kd_D)
###        CaOGB5₀  = OGB5_t - OGB5₀
###
###
###        c₀ = solve_Pepke_TR_eq(k_on_C1_T, k_off_C1_T, k_on_C2_T, k_off_C2_T, k_on_N1_T, k_off_N1_T, k_on_N2_T, k_off_N2_T, 
###                               k_on_C1_R, k_off_C1_R, k_on_C2_R, k_off_C2_R, k_on_N1_R, k_off_N1_R, k_on_N2_R, k_off_N2_R,
###                               Kd_0_TR, k0_TR, k0_RT, Kd_1_TR, k1_TR, k1_RT, Kd_2_TR, k2_TR, k2_RT, 
###                               Kd_3_TR, k3_TR, k3_RT, Kd_4_TR, k4_TR, k4_RT, Ca_free, CaM_t)
###
###        DMn_s₀_u   = (1-f_frac) * DMn₀   * U * 1
###        CaDMn_s₀_u = (1-f_frac) * CaDMn₀ * U * 1
###        DMn_f₀_u   =    f_frac  * DMn₀   * U * 1
###        CaDMn_f₀_u =    f_frac  * CaDMn₀ * U * 1
###    end
###
###    @init begin
###        DMn_s    = DMn_s₀_u 
###        CaDMn_s  = CaDMn_s₀_u
###        DMn_f    = DMn_f₀_u
###        CaDMn_f  = CaDMn_f₀_u
###        PP       = 0.0
###        CaPP     = 0.0
###        OGB5     = OGB5₀
###        CaOGB5   = CaOGB5₀
###        T0       = c₀[1]
###        T1C      = c₀[2] 
###        T1N      = c₀[3] 
###        T2C      = c₀[4] 
###        T1N1C    = c₀[5] 
###        T2N      = c₀[6] 
###        T1N2C    = c₀[7] 
###        T2N1C    = c₀[8] 
###        T4       = c₀[9] 
###        R0       = c₀[10]
###        R1C      = c₀[11] 
###        R1N      = c₀[12] 
###        R2C      = c₀[13] 
###        R1N1C    = c₀[14] 
###        R2N      = c₀[15] 
###        R1N2C    = c₀[16] 
###        R2N1C    = c₀[17] 
###        R4       = c₀[18]
###        Ca     = Ca_free
###    end
###
###    @dynamics begin
###        DMn_s'    = -kon_DMn * DMn_s * Ca + koff_DMn * CaDMn_s -   DMn_s * τ_s
###        CaDMn_s'  =  kon_DMn * DMn_s * Ca - koff_DMn * CaDMn_s - CaDMn_s * τ_s
###        DMn_f'    = -kon_DMn * DMn_f * Ca + koff_DMn * CaDMn_f -   DMn_f * τ_f
###        CaDMn_f'  =  kon_DMn * DMn_f * Ca - koff_DMn * CaDMn_f - CaDMn_f * τ_f
###        PP'       = -kon_PP  * PP    * Ca + koff_PP  * CaPP    +  2DMn_s * τ_s +  2DMn_f * τ_f + CaDMn_s * τ_s + CaDMn_f * τ_f
###        CaPP'     =  kon_PP  * PP    * Ca - koff_PP  * CaPP    + CaDMn_s * τ_s + CaDMn_f * τ_f
###        OGB5'     = -kon_D   * OGB5  * Ca + koff_D   * CaOGB5
###        CaOGB5'   =  kon_D   * OGB5  * Ca - koff_D   * CaOGB5
###        T0'       = -k_on_C1_T * T0 * Ca    + k_off_C1_T * T1C +    ### start of T states
###                    -k_on_N1_T * T0 * Ca    + k_off_N1_T * T1N +
###                    -k0_TR * T0     + k0_RT*R0
###        T1C'      =  k_on_C1_T * T0 * Ca    - k_off_C1_T * T1C +
###                    -k_on_C2_T * T1C * Ca   + k_off_C2_T * T2C +
###                    -k_on_N1_T * T1C * Ca   + k_off_N1_T * T1N1C +
###                    -k1_TR * T1C    + k1_RT*R1C
###        T1N'      =  k_on_N1_T * T0 * Ca    - k_off_N1_T * T1N +
###                    -k_on_N2_T * T1N * Ca   + k_off_N2_T * T2N +
###                    -k_on_C1_T * T1N * Ca   + k_off_C1_T * T1N1C +
###                    -k1_TR * T1N    + k1_RT*R1N
###        T2C'      =  k_on_C2_T * T1C * Ca   - k_off_C2_T * T2C +
###                    -k_on_N1_T * T2C * Ca   + k_off_N1_T * T1N2C +
###                    -k2_TR * T2C    + k2_RT*R2C
###        T1N1C'    =  k_on_N1_T * T1C * Ca   - k_off_N1_T * T1N1C +
###                     k_on_C1_T * T1N * Ca   - k_off_C1_T * T1N1C +
###                    -k_on_C2_T * T1N1C * Ca + k_off_C2_T * T1N2C +
###                    -k_on_N2_T * T1N1C * Ca + k_off_N2_T * T2N1C +
###                    -k2_TR * T1N1C  + k2_RT*R1N1C
###        T2N'      =  k_on_N2_T * T1N * Ca   - k_off_N2_T * T2N +
###                    -k_on_C1_T * T2N * Ca   + k_off_C1_T * T2N1C +
###                    -k2_TR * T2N  + k2_RT*R2N
###        T1N2C'    =  k_on_N1_T * T2C * Ca   - k_off_N1_T * T1N2C +
###                     k_on_C2_T * T1N1C * Ca - k_off_C2_T * T1N2C +
###                    -k_on_N2_T * T1N2C * Ca + k_off_N2_T * T4 +
###                    -k3_TR * T1N2C  + k3_RT*R1N2C
###        T2N1C'    =  k_on_C1_T * T2N * Ca   - k_off_C1_T * T2N1C +
###                     k_on_N2_T * T1N1C * Ca - k_off_N2_T * T2N1C +
###                    -k_on_C2_T * T2N1C * Ca + k_off_C2_T * T4 +
###                    -k3_TR * T2N1C  + k3_RT*R2N1C
###        T4'       =  k_on_N2_T * T1N2C * Ca - k_off_N2_T * T4 +
###                     k_on_C2_T * T2N1C * Ca - k_off_C2_T * T4 +
###                    -k4_TR * T4     + k4_RT*R4
###        R0'       = -k_on_C1_R * R0 * Ca    + k_off_C1_R * R1C +            ### start of R states
###                    -k_on_N1_R * R0 * Ca    + k_off_N1_R * R1N +
###                     k0_TR * T0     - k0_RT*R0
###        R1C'      =  k_on_C1_R * R0 * Ca    - k_off_C1_R * R1C +
###                    -k_on_C2_R * R1C * Ca   + k_off_C2_R * R2C +
###                    -k_on_N1_R * R1C * Ca   + k_off_N1_R * R1N1C +
###                     k1_TR * T1C    - k1_RT*R1C
###        R1N'      =  k_on_N1_R * R0 * Ca    - k_off_N1_R * R1N +
###                    -k_on_N2_R * R1N * Ca   + k_off_N2_R * R2N +
###                    -k_on_C1_R * R1N * Ca   + k_off_C1_R * R1N1C +
###                     k1_TR * T1N    - k1_RT*R1N
###        R2C'      =  k_on_C2_R * R1C * Ca   - k_off_C2_R * R2C +
###                    -k_on_N1_R * R2C * Ca   + k_off_N1_R * R1N2C +
###                     k2_TR * T2C    - k2_RT*R2C
###        R1N1C'    =  k_on_N1_R * R1C * Ca   - k_off_N1_R * R1N1C +
###                     k_on_C1_R * R1N * Ca   - k_off_C1_R * R1N1C +
###                    -k_on_C2_R * R1N1C * Ca + k_off_C2_R * R1N2C +
###                    -k_on_N2_R * R1N1C * Ca + k_off_N2_R * R2N1C +
###                     k2_TR * T1N1C  - k2_RT*R1N1C
###        R2N'      =  k_on_N2_R * R1N * Ca   - k_off_N2_R * R2N +
###                    -k_on_C1_R * R2N * Ca   + k_off_C1_R * R2N1C +
###                     k2_TR * T2N  - k2_RT*R2N
###        R1N2C'    =  k_on_N1_R * R2C * Ca   - k_off_N1_R * R1N2C +
###                     k_on_C2_R * R1N1C * Ca - k_off_C2_R * R1N2C +
###                    -k_on_N2_R * R1N2C * Ca + k_off_N2_R * R4 +
###                     k3_TR * T1N2C  - k3_RT*R1N2C
###        R2N1C'    =  k_on_C1_R * R2N * Ca   - k_off_C1_R * R2N1C +
###                     k_on_N2_R * R1N1C * Ca - k_off_N2_R * R2N1C +
###                    -k_on_C2_R * R2N1C * Ca + k_off_C2_R * R4 +
###                     k3_TR * T2N1C  - k3_RT*R2N1C
###        R4'       =  k_on_N2_R * R1N2C * Ca - k_off_N2_R * R4 +
###                     k_on_C2_R * R2N1C * Ca - k_off_C2_R * R4 +
###                     k4_TR * T4     - k4_RT*R4   
###        Ca'       = -kon_DMn * DMn_s * Ca + koff_DMn * CaDMn_s +
###                    -kon_DMn * DMn_f * Ca + koff_DMn * CaDMn_f +
###                    -kon_DMn * PP    * Ca + koff_PP  * CaPP    + 
###                    -kon_D   * OGB5  * Ca + koff_D   * CaOGB5  +
###                    -k_on_C1_T * T0 * Ca    + k_off_C1_T * T1C +    ### start of T states
###                    -k_on_N1_T * T0 * Ca    + k_off_N1_T * T1N +
###                    -k_on_C2_T * T1C * Ca   + k_off_C2_T * T2C +
###                    -k_on_N1_T * T1C * Ca   + k_off_C2_T * T1N1C +
###                    -k_on_N2_T * T1N * Ca   + k_off_N2_T * T2N +
###                    -k_on_C1_T * T1N * Ca   + k_off_C1_T * T1N1C +
###                    -k_on_N1_T * T2C * Ca   + k_off_N1_T * T1N2C +
###                    -k_on_C2_T * T1N1C * Ca + k_off_C2_T * T1N2C +
###                    -k_on_N2_T * T1N1C * Ca + k_off_N2_T * T2N1C +
###                    -k_on_C1_T * T2N * Ca   + k_off_C1_T * T2N1C +
###                    -k_on_N2_T * T1N2C * Ca + k_off_N2_T * T4 +
###                    -k_on_C2_T * T2N1C * Ca + k_off_C2_T * T4 +
###                    -k_on_C1_R * R0 * Ca    + k_off_C1_R * R1C +            ### start of R states
###                    -k_on_N1_R * R0 * Ca    + k_off_N1_R * R1N +
###                    -k_on_C2_R * R1C * Ca   + k_off_C2_R * R2C +
###                    -k_on_N1_R * R1C * Ca   + k_off_C2_R * R1N1C +
###                    -k_on_N2_R * R1N * Ca   + k_off_N2_R * R2N +
###                    -k_on_C1_R * R1N * Ca   + k_off_C1_R * R1N1C +
###                    -k_on_N1_R * R2C * Ca   + k_off_N1_R * R1N2C +
###                    -k_on_C2_R * R1N1C * Ca + k_off_C2_R * R1N2C +
###                    -k_on_N2_R * R1N1C * Ca + k_off_N2_R * R2N1C +
###                    -k_on_C1_R * R2N * Ca   + k_off_C1_R * R2N1C +
###                    -k_on_N2_R * R1N2C * Ca + k_off_N2_R * R4 +
###                    -k_on_C2_R * R2N1C * Ca + k_off_C2_R * R4
###
###
###    end
###
###    @derived begin
###        F_F0 = @. Normal((OGB5 + Fmax_Fmin * CaOGB5) / (OGB5₀ + Fmax_Fmin * CaOGB5₀), σ)
###    end
###
###    @observed begin
###        DMn_s   = DMn_s
###        CaDMn_s = CaDMn_s
###        DMn_f   = DMn_f
###        CaDMn_f = CaDMn_f
###        PP      = PP
###        CaPP    = CaPP
###        OGB5    = OGB5
###        CaOGB5  = CaOGB5
###        T0      = T0
###        T1C     = T1C
###        T1N     = T1N
###        T2C     = T2C
###        T1N1C   = T1N1C
###        T2N     = T2N
###        T1N2C   = T1N2C
###        T2N1C   = T2N1C
###        T4      = T4
###        R0      = R0
###        R1C     = R1C
###        R1N     = R1N
###        R2C     = R2C
###        R1N1C   = R1N1C
###        R2N     = R2N
###        R1N2C   = R1N2C
###        R2N1C   = R2N1C
###        R4      = R4
###        Ca      = Ca
###    end
###
###end;


###Stefan_model = @model begin
###    @options begin
###        inplace = true
###    end
###
###    @param begin
###        tv_kf_AT      ∈ RealDomain(; lower=2, init= 2.1, upper=12)
###        tv_Kd_AT      ∈ RealDomain(; lower=-8, init= -2.2, upper=-2)
###        tv_kf_BT      ∈ RealDomain(; lower=2, init= 2.3, upper=12)
###        tv_Kd_BT      ∈ RealDomain(; lower=-8, init= -2.4, upper=-2)
###        tv_kf_CT      ∈ RealDomain(; lower=2, init= 2.5, upper=12)
###        tv_Kd_CT      ∈ RealDomain(; lower=-8, init= -2.6, upper=-2)
###        tv_kf_DT      ∈ RealDomain(; lower=2, init= 2.7, upper=12)
###        tv_Kd_DT      ∈ RealDomain(; lower=-8, init= -2.8, upper=-2)
###
###        tv_kf_AR      ∈ RealDomain(; lower=2, init= 2.9, upper=12)
###        tv_Kd_AR      ∈ RealDomain(; lower=-8, init= -2.11, upper=-2)
###        tv_kf_BR      ∈ RealDomain(; lower=2, init= 2.12, upper=12)
###        tv_Kd_BR      ∈ RealDomain(; lower=-8, init= -2.13, upper=-2)
###        tv_kf_CR      ∈ RealDomain(; lower=2, init= 2.14, upper=12)
###        tv_Kd_CR      ∈ RealDomain(; lower=-8, init= -2.15, upper=-2)
###        tv_kf_DR      ∈ RealDomain(; lower=2, init= 2.16, upper=12)
###        tv_Kd_DR      ∈ RealDomain(; lower=-8, init= -2.17, upper=-2)
###
###        tv_k0_TR      ∈ RealDomain(; lower=2, init= 2.18, upper=12)
###        tv_Kd_0_TR    ∈ RealDomain(; lower=-8, init= -2.19, upper=-2)
###        tv_k1_TR      ∈ RealDomain(; lower=2, init= 2.20, upper=12)
###        tv_Kd_1_TR    ∈ RealDomain(; lower=-8, init= -2.21, upper=-2)
###        tv_k2_TR      ∈ RealDomain(; lower=2, init= 2.22, upper=12)
###        tv_Kd_2_TR    ∈ RealDomain(; lower=-8, init= -2.23, upper=-2)
###        tv_k3_TR      ∈ RealDomain(; lower=2, init= 2.24, upper=12)
###        tv_Kd_3_TR    ∈ RealDomain(; lower=-8, init= -2.25, upper=-2)
###        tv_k4_TR      ∈ RealDomain(; lower=2, init= 2.26, upper=12)
###        tv_Kd_4_TR    ∈ RealDomain(; lower=-8, init= -2.27, upper=-2)
###
###        tv_mα         ∈ RealDomain(; lower=1e-10, init= 0.001, upper=1)
###        tv_α₀         ∈ RealDomain(;              init = -0.39)
###
######        Ω             ∈ Constrained(Normal(1), lower=1e-10, upper=1)
###        σ             ∈ RealDomain(; lower=1e-10, init=10)
###    end
###
###    @random begin
###        ###η ~ Normal(Ω)
###        η ~ Normal(0.1)
###    end
###    
###    @covariates begin
###        cov_f_frac  
###        cov_τ_f  
###        cov_τ_s  
###        cov_kon_DMn  
###        cov_koff_DMn  
###        cov_Kd_DMn  
###        cov_koff_PP  
###        cov_DMn_t  
###        cov_CaM_t  
###        cov_OGB5_t  
###        cov_Ca_free
###        PCD
###    end
###
###    @pre begin
###        U         = (tv_mα*PCD + tv_α₀) * (1 + η[1])
###        Fmax_Fmin = 39.364
###
###        f_frac   = cov_f_frac 
###        τ_f      = 1/cov_τ_f
###        τ_s      = 1/cov_τ_s
###        kon_DMn  = cov_kon_DMn
###        koff_DMn = cov_koff_DMn
###        Kd_DMn   = cov_Kd_DMn
###        kon_PP   = cov_kon_DMn
###        koff_PP  = cov_koff_PP
###
###        DMn_t   = cov_DMn_t
###        CaM_t   = cov_CaM_t
###        OGB5_t  = cov_OGB5_t
###        Ca_free = cov_Ca_free
###
###        kon_D    = 8.77e5         ### in M^-1 ms^-1
###        koff_D   = 33.2          ### in ms^-1 
###        Kd_D     = koff_D / kon_D
###
###        Kd_AT  = 10^tv_Kd_AT 
###        kf_AT   = 10^tv_kf_AT
###        kb_AT   = Kd_AT * kf_AT 
###
###        Kd_BT  = 10^tv_Kd_BT 
###        kf_BT   = 10^tv_kf_BT
###        kb_BT   = Kd_BT * kf_BT 
###
###        Kd_CT  = 10^tv_Kd_CT 
###        kf_CT   = 10^tv_kf_CT
###        kb_CT   = Kd_CT * kf_CT 
###
###        Kd_DT  = 10^tv_Kd_DT 
###        kf_DT   = 10^tv_kf_DT
###        kb_DT   = Kd_DT * kf_DT 
###
###        Kd_AR  = 10^tv_Kd_AR 
###        kf_AR   = 10^tv_kf_AR
###        kb_AR   = Kd_AR * kf_AR 
###
###        Kd_BR  = 10^tv_Kd_BR 
###        kf_BR   = 10^tv_kf_BR
###        kb_BR   = Kd_BR * kf_BR 
###
###        Kd_CR  = 10^tv_Kd_CR 
###        kf_CR   = 10^tv_kf_CR
###        kb_CR   = Kd_CR * kf_CR 
###
###        Kd_DR  = 10^tv_Kd_DR 
###        kf_DR   = 10^tv_kf_DR
###        kb_DR   = Kd_DR * kf_DR
###
###        Kd_0_TR = 10^tv_Kd_0_TR 
###        k0_TR   = 10^tv_k0_TR
###        k0_RT   = Kd_0_TR * k0_TR 
###
###        Kd_1_TR = 10^tv_Kd_1_TR 
###        k1_TR   = 10^tv_k1_TR
###        k1_RT   = Kd_1_TR * k1_TR
###
###        Kd_2_TR = 10^tv_Kd_2_TR 
###        k2_TR   = 10^tv_k2_TR
###        k2_RT   = Kd_2_TR * k2_TR
###
###        Kd_3_TR = 10^tv_Kd_3_TR 
###        k3_TR   = 10^tv_k3_TR
###        k3_RT   = Kd_3_TR * k3_TR
###
###        Kd_4_TR = 10^tv_Kd_4_TR 
###        k4_TR   = 10^tv_k4_TR
###        k4_RT   = Kd_4_TR * k4_TR
###
###        DMn₀     = (Kd_DMn * DMn_t)  / (Ca_free + Kd_DMn)
###        CaDMn₀   = DMn_t - DMn₀
###        OGB5₀    = (Kd_D   * OGB5_t) / (Ca_free + Kd_D)
###        CaOGB5₀  = OGB5_t - OGB5₀
###
######        A = [-k0_TR - Ca_free*kf_AT - Ca_free*kf_BT - Ca_free*kf_CT - Ca_free*kf_DT kb_AT kb_BT kb_CT kb_DT 0 0 0 0 0 0 0 0 0 0 0 k0_RT 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0; 
######              Ca_free*kf_AT -k1_TR - kb_AT - Ca_free*kf_BT - Ca_free*kf_CT - Ca_free*kf_DT 0 0 0 kb_BT kb_CT kb_DT 0 0 0 0 0 0 0 0 0 k1_RT 0 0 0 0 0 0 0 0 0 0 0 0 0 0; 
######              Ca_free*kf_BT 0 -k1_TR - kb_BT - Ca_free*kf_AT - Ca_free*kf_CT - Ca_free*kf_DT 0 0 kb_AT 0 0 kb_CT kb_DT 0 0 0 0 0 0 0 0 k1_RT 0 0 0 0 0 0 0 0 0 0 0 0 0; 
######              Ca_free*kf_CT 0 0 -k1_TR - kb_CT - Ca_free*kf_AT - Ca_free*kf_BT - Ca_free*kf_DT 0 0 kb_AT 0 kb_BT 0 kb_DT 0 0 0 0 0 0 0 0 k1_RT 0 0 0 0 0 0 0 0 0 0 0 0; 
######              Ca_free*kf_DT 0 0 0 -k1_TR - kb_DT - Ca_free*kf_AT - Ca_free*kf_BT - Ca_free*kf_CT 0 0 kb_AT 0 kb_BT kb_CT 0 0 0 0 0 0 0 0 0 k1_RT 0 0 0 0 0 0 0 0 0 0 0; 
######              0 Ca_free*kf_BT Ca_free*kf_AT 0 0 -k2_TR - kb_AT - kb_BT - Ca_free*kf_CT - Ca_free*kf_DT 0 0 0 0 0 kb_CT kb_DT 0 0 0 0 0 0 0 0 k2_RT 0 0 0 0 0 0 0 0 0 0; 
######              0 Ca_free*kf_CT 0 Ca_free*kf_AT 0 0 -k2_TR - kb_AT - kb_CT - Ca_free*kf_BT - Ca_free*kf_DT 0 0 0 0 kb_BT 0 kb_DT 0 0 0 0 0 0 0 0 k2_RT 0 0 0 0 0 0 0 0 0; 
######              0 Ca_free*kf_DT 0 0 Ca_free*kf_AT 0 0 -k2_TR - kb_AT - kb_DT - Ca_free*kf_BT - Ca_free*kf_CT 0 0 0 kb_BT 0 kb_CT 0 0 0 0 0 0 0 0 0 k2_RT 0 0 0 0 0 0 0 0; 
######              0 0 Ca_free*kf_CT Ca_free*kf_BT 0 0 0 0 -k2_TR - kb_BT - kb_CT - Ca_free*kf_AT - Ca_free*kf_DT 0 0 kb_AT 0 0 kb_DT 0 0 0 0 0 0 0 0 0 k2_RT 0 0 0 0 0 0 0; 
######              0 0 Ca_free*kf_DT 0 Ca_free*kf_BT 0 0 0 -Ca_free*kf_AT - Ca_free*kf_DT -k2_TR - kb_BT - kb_DT 0 kb_AT 0 0 kb_DT 0 0 0 0 0 0 0 0 0 0 k2_RT 0 0 0 0 0 0; 
######              0 0 0 Ca_free*kf_DT Ca_free*kf_CT 0 0 0 -Ca_free*kf_AT - Ca_free*kf_DT 0 -k2_TR - kb_CT - kb_DT kb_AT 0 0 kb_DT 0 0 0 0 0 0 0 0 0 0 0 k2_RT 0 0 0 0 0; 
######              0 0 0 0 0 Ca_free*kf_CT Ca_free*kf_BT 0 Ca_free*kf_AT 0 0 -k3_TR - kb_AT - kb_BT - kb_CT - Ca_free*kf_DT 0 0 0 kb_DT 0 0 0 0 0 0 0 0 0 0 0 k1_RT 0 0 0 0; 
######              0 0 0 0 0 Ca_free*kf_DT 0 Ca_free*kf_BT 0 Ca_free*kf_AT 0 0 -k3_TR - kb_AT - kb_BT - kb_DT - Ca_free*kf_CT 0 0 kb_CT 0 0 0 0 0 0 0 0 0 0 0 0 k3_RT 0 0 0; 
######              0 0 0 0 0 0 Ca_free*kf_DT Ca_free*kf_CT 0 0 Ca_free*kf_AT 0 0 -k3_TR - kb_AT - kb_CT - kb_DT - Ca_free*kf_BT 0 kb_BT 0 0 0 0 0 0 0 0 0 0 0 0 0 k3_RT 0 0; 
######              0 0 0 0 0 0 0 0 Ca_free*kf_DT Ca_free*kf_CT Ca_free*kf_BT 0 0 -kb_BT - kb_CT - kb_DT -k3_TR - Ca_free*kf_AT kb_AT 0 0 0 0 0 0 0 0 0 0 0 0 0 0 k3_RT 0; 
######              0 0 0 0 0 0 0 0 0 0 0 Ca_free*kf_DT Ca_free*kf_CT Ca_free*kf_BT Ca_free*kf_AT -k4_TR - kb_AT - kb_BT - kb_CT - kb_DT 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 k4_RT; 
######              k0_TR 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -k0_RT - Ca_free*kf_AR - Ca_free*kf_BR - Ca_free*kf_CR - Ca_free*kf_DR kb_AR kb_BR kb_CR kb_DR 0 0 0 0 0 0 0 0 0 0 0; 
######              0 k1_TR 0 0 0 0 0 0 0 0 0 0 0 0 0 0 Ca_free*kf_AR -k1_RT - kb_AR - Ca_free*kf_BR - Ca_free*kf_CR - Ca_free*kf_DR 0 0 0 kb_BR kb_CR kb_DR 0 0 0 0 0 0 0 0; 
######              0 0 k1_TR 0 0 0 0 0 0 0 0 0 0 0 0 0 Ca_free*kf_BR 0 -k1_RT - kb_BR - Ca_free*kf_AR - Ca_free*kf_CR - Ca_free*kf_DR 0 0 kb_AR 0 0 kb_CR kb_DR 0 0 0 0 0 0; 
######              0 0 0 k1_TR 0 0 0 0 0 0 0 0 0 0 0 0 Ca_free*kf_CR 0 0 -k1_RT - kb_CR - Ca_free*kf_AR - Ca_free*kf_BR - Ca_free*kf_DR 0 0 kb_AR 0 kb_BR 0 kb_DR 0 0 0 0 0; 
######              0 0 0 0 k1_TR 0 0 0 0 0 0 0 0 0 0 0 Ca_free*kf_DR 0 0 0 -k1_RT - kb_DR - Ca_free*kf_AR - Ca_free*kf_BR - Ca_free*kf_CR 0 0 kb_AR 0 kb_BR kb_CR 0 0 0 0 0; 
######              0 0 0 0 0 k2_TR 0 0 0 0 0 0 0 0 0 0 0 Ca_free*kf_BR Ca_free*kf_AR 0 0 -k2_RT - kb_AR - kb_BR - Ca_free*kf_CR - Ca_free*kf_DR 0 0 0 0 0 kb_CR kb_DR 0 0 0; 
######              0 0 0 0 0 0 k2_TR 0 0 0 0 0 0 0 0 0 0 Ca_free*kf_CR 0 Ca_free*kf_AR 0 0 -k2_RT - kb_AR - kb_CR - Ca_free*kf_BR - Ca_free*kf_DR 0 0 0 0 kb_BR 0 kb_DR 0 0; 
######              0 0 0 0 0 0 0 k2_TR 0 0 0 0 0 0 0 0 0 Ca_free*kf_DR 0 0 Ca_free*kf_AR 0 0 -k2_RT - kb_AR - kb_DR - Ca_free*kf_BR - Ca_free*kf_CR 0 0 0 kb_BR 0 kb_CR 0 0; 
######              0 0 0 0 0 0 0 0 k2_TR 0 0 0 0 0 0 0 0 0 Ca_free*kf_CR Ca_free*kf_BR 0 0 0 0 -k2_RT - kb_BR - kb_CR - Ca_free*kf_AR - Ca_free*kf_DR 0 0 kb_AR 0 0 kb_DR 0; 
######              0 0 0 0 0 0 0 0 0 k2_TR 0 0 0 0 0 0 0 0 Ca_free*kf_DR 0 Ca_free*kf_BR 0 0 0 -Ca_free*kf_AR - Ca_free*kf_DR -k2_RT - kb_BR - kb_DR 0 kb_AR 0 0 kb_DR 0; 
######              0 0 0 0 0 0 0 0 0 0 k2_TR 0 0 0 0 0 0 0 0 Ca_free*kf_DR Ca_free*kf_CR 0 0 0 -Ca_free*kf_AR - Ca_free*kf_DR 0 -k2_RT - kb_CR - kb_DR kb_AR 0 0 kb_DR 0; 
######              0 0 0 0 0 0 0 0 0 0 0 k3_TR 0 0 0 0 0 0 0 0 0 Ca_free*kf_CR Ca_free*kf_BR 0 Ca_free*kf_AR 0 0 -k3_RT - kb_AR - kb_BR - kb_CR - Ca_free*kf_DR 0 0 0 kb_DR;
######              0 0 0 0 0 0 0 0 0 0 0 0 k3_TR 0 0 0 0 0 0 0 0 Ca_free*kf_DR 0 Ca_free*kf_BR 0 Ca_free*kf_AR 0 0 -k3_RT - kb_AR - kb_BR - kb_DR - Ca_free*kf_CR 0 0 kb_CR; 
######              0 0 0 0 0 0 0 0 0 0 0 0 0 k3_TR 0 0 0 0 0 0 0 0 Ca_free*kf_DR Ca_free*kf_CR 0 0 Ca_free*kf_AR 0 0 -k3_RT - kb_AR - kb_CR - kb_DR - Ca_free*kf_BR 0 kb_BR; 
######              0 0 0 0 0 0 0 0 0 0 0 0 0 0 k3_TR 0 0 0 0 0 0 0 0 0 Ca_free*kf_DR Ca_free*kf_CR Ca_free*kf_BR 0 0 -kb_BR - kb_CR - kb_DR -k3_RT - Ca_free*kf_AR kb_AR;
#########              0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 k4_TR 0 0 0 0 0 0 0 0 0 0 0 Ca*kf_DR Ca*kf_CR Ca*kf_BR Ca*kf_AR k4_RT - kb_AR - kb_BR - kb_CR - kb_DR;
######              1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1] 
######
######        b = [zeros(31); CaM_t]
######
######        c₀ = ones(32) # A\b 
###
###        c₀ = solve_Stefan_TR_eq(Kd_AT, kf_AT, kb_AT, Kd_BT, kf_BT, kb_BT, Kd_CT, kf_CT, kb_CT, Kd_DT, kf_DT, kb_DT, 
###                            Kd_AR, kf_AR, kb_AR, Kd_BR, kf_BR, kb_BR, Kd_CR, kf_CR, kb_CR, Kd_DR, kf_DR, kb_DR, 
###                            Kd_0_TR, k0_TR, k0_RT, Kd_1_TR, k1_TR, k1_RT, Kd_2_TR, k2_TR, k2_RT, 
###                            Kd_3_TR, k3_TR, k3_RT, Kd_4_TR, k4_TR, k4_RT, Ca_free, CaM_t)
###
######        T0₀     = c₀[1]
######        T_A₀    = c₀[2] 
######        T_B₀    = c₀[3] 
######        T_C₀    = c₀[4] 
######        T_D₀    = c₀[5] 
######        T_AB₀   = c₀[6] 
######        T_AC₀   = c₀[7] 
######        T_AD₀   = c₀[8] 
######        T_BC₀   = c₀[9] 
######        T_BD₀   = c₀[10] 
######        T_CD₀   = c₀[11] 
######        T_ABC₀  = c₀[12] 
######        T_ABD₀  = c₀[13] 
######        T_ACD₀  = c₀[14] 
######        T_BCD₀  = c₀[15] 
######        T_ABCD₀ = c₀[16]
######        R0₀     = c₀[17]
######        R_A₀    = c₀[18] 
######        R_B₀    = c₀[19] 
######        R_C₀    = c₀[20] 
######        R_D₀    = c₀[21] 
######        R_AB₀   = c₀[22] 
######        R_AC₀   = c₀[23] 
######        R_AD₀   = c₀[24] 
######        R_BC₀   = c₀[25] 
######        R_BD₀   = c₀[26] 
######        R_CD₀   = c₀[27] 
######        R_ABC₀  = c₀[28] 
######        R_ABD₀  = c₀[29] 
######        R_ACD₀  = c₀[30] 
######        R_BCD₀  = c₀[31] 
######        R_ABCD₀ = c₀[32]
###
###        DMn_s₀_u   = (1-f_frac) * DMn₀   * U * 1
###        CaDMn_s₀_u = (1-f_frac) * CaDMn₀ * U * 1
###        DMn_f₀_u   =    f_frac  * DMn₀   * U * 1
###        CaDMn_f₀_u =    f_frac  * CaDMn₀ * U * 1
###    end
###
###    @init begin
###        DMn_s    = DMn_s₀_u 
###        CaDMn_s  = CaDMn_s₀_u
###        DMn_f    = DMn_f₀_u
###        CaDMn_f  = CaDMn_f₀_u
###        PP       = 0.0
###        CaPP     = 0.0
###        OGB5     = OGB5₀
###        CaOGB5   = CaOGB5₀
####        T0       = T0₀
####        T_A      = T_A₀
####        T_B      = T_B₀
####        T_C      = T_C₀
####        T_D      = T_D₀
####        T_AB     = T_AB₀
####        T_AC     = T_AC₀
####        T_AD     = T_AD₀
####        T_BC     = T_BC₀
####        T_BD     = T_BD₀
####        T_CD     = T_CD₀
####        T_ABC    = T_ABC₀
####        T_ABD    = T_ABD₀
####        T_ACD    = T_ACD₀
####        T_BCD    = T_BCD₀
####        T_ABCD   = T_ABCD₀
####        R0       = R0₀
####        R_A      = R_A₀
####        R_B      = R_B₀
####        R_C      = R_C₀
####        R_D      = R_D₀
####        R_AB     = R_AB₀
####        R_AC     = R_AC₀
####        R_AD     = R_AD₀
####        R_BC     = R_BC₀
####        R_BD     = R_BD₀
####        R_CD     = R_CD₀
####        R_ABC    = R_ABC₀
####        R_ABD    = R_ABD₀
####        R_ACD    = R_ACD₀
####        R_BCD    = R_BCD₀
####        R_ABCD   = R_ABCD₀
###        T0     = c₀[1]
###        T_A    = c₀[2] 
###        T_B    = c₀[3] 
###        T_C    = c₀[4] 
###        T_D    = c₀[5] 
###        T_AB   = c₀[6] 
###        T_AC   = c₀[7] 
###        T_AD   = c₀[8] 
###        T_BC   = c₀[9] 
###        T_BD   = c₀[10] 
###        T_CD   = c₀[11] 
###        T_ABC  = c₀[12] 
###        T_ABD  = c₀[13] 
###        T_ACD  = c₀[14] 
###        T_BCD  = c₀[15] 
###        T_ABCD = c₀[16]
###        R0     = c₀[17]
###        R_A    = c₀[18] 
###        R_B    = c₀[19] 
###        R_C    = c₀[20] 
###        R_D    = c₀[21] 
###        R_AB   = c₀[22] 
###        R_AC   = c₀[23] 
###        R_AD   = c₀[24] 
###        R_BC   = c₀[25] 
###        R_BD   = c₀[26] 
###        R_CD   = c₀[27] 
###        R_ABC  = c₀[28] 
###        R_ABD  = c₀[29] 
###        R_ACD  = c₀[30] 
###        R_BCD  = c₀[31] 
###        R_ABCD = c₀[32]
###        Ca     = Ca_free
###    end
###
###    @dynamics begin
###        DMn_s'    = -kon_DMn * DMn_s * Ca + koff_DMn * CaDMn_s -   DMn_s * τ_s
###        CaDMn_s'  =  kon_DMn * DMn_s * Ca - koff_DMn * CaDMn_s - CaDMn_s * τ_s
###        DMn_f'    = -kon_DMn * DMn_f * Ca + koff_DMn * CaDMn_f -   DMn_f * τ_f
###        CaDMn_f'  =  kon_DMn * DMn_f * Ca - koff_DMn * CaDMn_f - CaDMn_f * τ_f
###        PP'       = -kon_PP  * PP    * Ca + koff_PP  * CaPP    +  2DMn_s * τ_s +  2DMn_f * τ_f + CaDMn_s * τ_s + CaDMn_f * τ_f
###        CaPP'     =  kon_PP  * PP    * Ca - koff_PP  * CaPP    + CaDMn_s * τ_s + CaDMn_f * τ_f
###        OGB5'     = -kon_D   * OGB5  * Ca + koff_D   * CaOGB5
###        CaOGB5'   =  kon_D   * OGB5  * Ca - koff_D   * CaOGB5
###        T0' = -kf_AT * T0 * Ca + kb_AT * T_A +  ### T0 ### start of T states
###        -kf_BT * T0 * Ca + kb_BT * T_B +
###        -kf_CT * T0 * Ca + kb_CT * T_C +
###        -kf_DT * T0 * Ca + kb_DT * T_D +
###        -T0 * k0_TR + R0 * k0_RT
###        T_A' =  kf_AT * T0 * Ca - kb_AT * T_A +  ### T_A
###        -kf_BT * T_A * Ca + kb_BT * T_AB +
###        -kf_CT * T_A * Ca + kb_CT * T_AC +
###        -kf_DT * T_A * Ca + kb_DT * T_AD +
###        -T_A * k1_TR + R_A * k1_RT
###        T_B' =  kf_BT * T0 * Ca - kb_BT * T_B +  ### T_B
###        -kf_AT * T_B * Ca + kb_AT * T_AB +
###        -kf_CT * T_B * Ca + kb_CT * T_BC +
###        -kf_DT * T_B * Ca + kb_DT * T_BD +
###        -T_B * k1_TR + R_B * k1_RT
###        T_C' =  kf_CT * T0 * Ca - kb_CT * T_C +  ### T_C
###        -kf_AT * T_C * Ca + kb_AT * T_AC +
###        -kf_BT * T_C * Ca + kb_BT * T_BC +
###        -kf_DT * T_C * Ca + kb_DT * T_CD +
###        -T_C * k1_TR + R_C * k1_RT
###        T_D' =  kf_DT * T0 * Ca - kb_DT * T_D +  ### T_D
###        -kf_AT * T_D * Ca + kb_AT * T_AD +
###        -kf_BT * T_D * Ca + kb_BT * T_BD +
###        -kf_CT * T_D * Ca + kb_CT * T_CD +
###        -T_D * k1_TR + R_D * k1_RT
###        T_AB' =  kf_AT * T_B * Ca - kb_AT * T_AB + ### T_AB
###         kf_BT * T_A * Ca - kb_BT * T_AB +
###        -kf_CT * T_AB * Ca + kb_CT * T_ABC +
###        -kf_DT * T_AB * Ca + kb_DT * T_ABD +
###        -T_AB * k2_TR + R_AB * k2_RT
###        T_AC' =  kf_AT * T_C * Ca - kb_AT * T_AC + ### T_AC
###         kf_CT * T_A * Ca - kb_CT * T_AC +
###        -kf_BT * T_AC * Ca + kb_BT * T_ABC +
###        -kf_DT * T_AC * Ca + kb_DT * T_ACD +
###        -T_AC * k2_TR + R_AC * k2_RT
###        T_AD' =  kf_AT * T_D * Ca - kb_AT * T_AD + ### T_AD
###         kf_DT * T_A * Ca - kb_DT * T_AD +
###        -kf_BT * T_AD * Ca + kb_BT * T_ABC +
###        -kf_CT * T_AD * Ca + kb_CT * T_ACD +
###        -T_AD * k2_TR + R_AD * k2_RT
###        T_BC' =  kf_BT * T_C * Ca - kb_BT * T_BC + ### T_BC
###         kf_CT * T_B * Ca - kb_CT * T_BC +
###        -kf_AT * T_BC * Ca + kb_AT * T_ABC +
###        -kf_DT * T_BC * Ca + kb_DT * T_BCD +
###        -T_BC * k2_TR + R_BC * k2_RT
###        T_BD' =  kf_BT * T_D * Ca - kb_BT * T_BD + ### T_BD
###         kf_DT * T_B * Ca - kb_DT * T_BD +
###        -kf_AT * T_BC * Ca + kb_AT * T_ABC +
###        -kf_DT * T_BC * Ca + kb_DT * T_BCD +
###        -T_BD * k2_TR + R_BD * k2_RT
###        T_CD' =  kf_CT * T_D * Ca - kb_CT * T_CD + ### T_CD
###         kf_DT * T_C * Ca - kb_DT * T_CD +
###        -kf_AT * T_BC * Ca + kb_AT * T_ABC +
###        -kf_DT * T_BC * Ca + kb_DT * T_BCD +
###        -T_CD * k2_TR + R_CD * k2_RT
###        T_ABC' =  kf_AT * T_BC * Ca - kb_AT * T_ABC + ### T_ABC
###         kf_BT * T_AC * Ca - kb_BT * T_ABC +
###         kf_CT * T_AB * Ca - kb_CT * T_ABC +
###        -kf_DT * T_ABC * Ca + kb_DT * T_ABCD +
###        -T_ABC * k3_TR + R_ABC * k1_RT
###        T_ABD' =  kf_AT * T_BD * Ca - kb_AT * T_ABD + ### T_ABD
###         kf_BT * T_AD * Ca - kb_BT * T_ABD +
###         kf_DT * T_AB * Ca - kb_DT * T_ABD +
###        -kf_CT * T_ABD * Ca + kb_CT * T_ABCD +
###        -T_ABD * k3_TR + R_ABD * k3_RT
###        T_ACD' =  kf_AT * T_CD * Ca - kb_AT * T_ACD + ### T_ACD
###         kf_CT * T_AD * Ca - kb_CT * T_ACD +
###         kf_DT * T_AC * Ca - kb_DT * T_ACD +
###        -kf_BT * T_ACD * Ca + kb_BT * T_ABCD +
###        -T_ACD * k3_TR + R_ACD * k3_RT
###        T_BCD' =  kf_BT * T_CD * Ca - kb_BT * T_ACD + ### T_BCD
###         kf_CT * T_BD * Ca - kb_CT * T_ACD +
###         kf_DT * T_BC * Ca - kb_DT * T_ACD +
###        -kf_AT * T_BCD * Ca + kb_AT * T_ABCD +
###        -T_BCD * k3_TR + R_BCD * k3_RT
###        T_ABCD' =  kf_AT * T_BCD * Ca - kb_AT * T_ABCD + ### T_ABCD
###         kf_BT * T_ACD * Ca - kb_BT * T_ABCD +
###         kf_CT * T_ABD * Ca - kb_CT * T_ABCD +
###         kf_DT * T_ABC * Ca - kb_DT * T_ABCD +
###        -T_ABCD * k4_TR + R_ABCD * k4_RT
###        R0' = -kf_AR * R0 * Ca + kb_AR * R_A + ### R0        #### start of R states
###        -kf_BR * R0 * Ca + kb_BR * R_B +
###        -kf_CR * R0 * Ca + kb_CR * R_C +
###        -kf_DR * R0 * Ca + kb_DR * R_D +
###         T0 * k0_TR - R0 * k0_RT
###        R_A' =  kf_AR * R0 * Ca - kb_AR * R_A + ### R_A
###        -kf_BR * R_A * Ca + kb_BR * R_AB +
###        -kf_CR * R_A * Ca + kb_CR * R_AC +
###        -kf_DR * R_A * Ca + kb_DR * R_AD +
###         T_A * k1_TR - R_A * k1_RT
###        R_B' =  kf_BR * R0 * Ca - kb_BR * R_B + ### R_B
###        -kf_AR * R_B * Ca + kb_AR * R_AB +
###        -kf_CR * R_B * Ca + kb_CR * R_BC +
###        -kf_DR * R_B * Ca + kb_DR * R_BD +
###         T_B * k1_TR - R_B * k1_RT
###        R_C' =  kf_CR * R0 * Ca - kb_CR * R_C + ### R_C
###        -kf_AR * R_C * Ca + kb_AR * R_AC +
###        -kf_BR * R_C * Ca + kb_BR * R_BC +
###        -kf_DR * R_C * Ca + kb_DR * R_CD +
###         T_C * k1_TR - R_C * k1_RT
###        R_D' =  kf_DR * R0 * Ca - kb_DR * R_D + ### R_D
###        -kf_AR * R_D * Ca + kb_AR * R_AD +
###        -kf_BR * R_D * Ca + kb_BR * R_BD +
###        -kf_CR * R_D * Ca + kb_CR * R_CD +
###         T_D * k1_TR - R_D * k1_RT
###        R_AB' =  kf_AR * R_B * Ca - kb_AR * R_AB + ### R_AB
###         kf_BR * R_A * Ca - kb_BR * R_AB +
###        -kf_CR * R_AB * Ca + kb_CR * R_ABC +
###        -kf_DR * R_AB * Ca + kb_DR * R_ABD +
###         T_AB * k2_TR - R_AB * k2_RT
###        R_AC' =  kf_AR * R_C * Ca - kb_AR * R_AC + ### R_AC
###         kf_CR * R_A * Ca - kb_CR * R_AC +
###        -kf_BR * R_AC * Ca + kb_BR * R_ABC +
###        -kf_DR * R_AC * Ca + kb_DR * R_ACD +
###         T_AC * k2_TR - R_AC * k2_RT
###        R_AD' =  kf_AR * R_D * Ca - kb_AR * R_AD + ### R_AD
###         kf_DR * R_A * Ca - kb_DR * R_AD +
###        -kf_BR * R_AD * Ca + kb_BR * R_ABC +
###        -kf_CR * R_AD * Ca + kb_CR * R_ACD +
###         T_AD * k2_TR - R_AD * k2_RT
###        R_BC' =  kf_BR * R_C * Ca - kb_BR * R_BC + ### R_BC
###         kf_CR * R_B * Ca - kb_CR * R_BC +
###        -kf_AR * R_BC * Ca + kb_AR * R_ABC +
###        -kf_DR * R_BC * Ca + kb_DR * R_BCD +
###         T_BC * k2_TR - R_BC * k2_RT
###        R_BD' =  kf_BR * R_D * Ca - kb_BR * R_BD + ### R_BD
###         kf_DR * R_B * Ca - kb_DR * R_BD +
###        -kf_AR * R_BC * Ca + kb_AR * R_ABC +
###        -kf_DR * R_BC * Ca + kb_DR * R_BCD +
###         T_BD * k2_TR - R_BD * k2_RT
###        R_CD' =  kf_CR * R_D * Ca - kb_CR * R_CD + ### R_CD
###         kf_DR * R_C * Ca - kb_DR * R_CD +
###        -kf_AR * R_BC * Ca + kb_AR * R_ABC +
###        -kf_DR * R_BC * Ca + kb_DR * R_BCD +
###         T_CD * k2_TR - R_CD * k2_RT
###        R_ABC' =  kf_AR * R_BC * Ca - kb_AR * R_ABC + ### R_ABC
###         kf_BR * R_AC * Ca - kb_BR * R_ABC +
###         kf_CR * R_AB * Ca - kb_CR * R_ABC +
###        -kf_DR * R_ABC * Ca + kb_DR * R_ABCD +
###         T_ABC * k3_TR - R_ABC * k3_RT
###        R_ABD' =  kf_AR * R_BD * Ca - kb_AR * R_ABD + ### R_ABD
###         kf_BR * R_AD * Ca - kb_BR * R_ABD +
###         kf_DR * R_AB * Ca - kb_DR * R_ABD +
###        -kf_CR * R_ABD * Ca + kb_CR * R_ABCD +
###         T_ABD * k3_TR - R_ABD * k3_RT
###        R_ACD' =  kf_AR * R_CD * Ca - kb_AR * R_ACD + ### R_ACD
###         kf_CR * R_AD * Ca - kb_CR * R_ACD +
###         kf_DR * R_AC * Ca - kb_DR * R_ACD +
###        -kf_BR * R_ACD * Ca + kb_BR * R_ABCD +
###         T_ACD * k3_TR - R_ACD * k3_RT
###        R_BCD' =  kf_BR * R_CD * Ca - kb_BR * R_ACD + ### R_BCD
###         kf_CR * R_BD * Ca - kb_CR * R_ACD +
###         kf_DR * R_BC * Ca - kb_DR * R_ACD +
###        -kf_AR * R_BCD * Ca + kb_AR * R_ABCD +
###         T_BCD * k3_TR - R_BCD * k3_RT
###        R_ABCD' =  kf_AR * R_BCD * Ca - kb_AR * R_ABCD + ### R_ABCD
###         kf_BR * R_ACD * Ca - kb_BR * R_ABCD +
###         kf_CR * R_ABD * Ca - kb_CR * R_ABCD +
###         kf_DR * R_ABC * Ca - kb_DR * R_ABCD +
###         T_ABCD * k4_TR + R_ABCD * k4_RT        
###        Ca'       = -kon_DMn * DMn_s * Ca + koff_DMn * CaDMn_s +
###                    -kon_DMn * DMn_f * Ca + koff_DMn * CaDMn_f +
###                    -kon_DMn * PP    * Ca + koff_PP  * CaPP    + 
###                    -kon_D   * OGB5  * Ca + koff_D   * CaOGB5  +
###                    -kf_AT * T0 * Ca + kb_AT * T_A +        
###                    -kf_BT * T0 * Ca + kb_BT * T_B +
###                    -kf_CT * T0 * Ca + kb_CT * T_C +
###                    -kf_DT * T0 * Ca + kb_DT * T_D +
###                    -kf_BT * T_A * Ca + kb_BT * T_AB +
###                    -kf_CT * T_A * Ca + kb_CT * T_AC +
###                    -kf_DT * T_A * Ca + kb_DT * T_AD +
###                    -kf_AT * T_B * Ca + kb_AT * T_AB +
###                    -kf_CT * T_B * Ca + kb_CT * T_BC +
###                    -kf_DT * T_B * Ca + kb_DT * T_BD +
###                    -kf_AT * T_C * Ca + kb_AT * T_AC +
###                    -kf_BT * T_C * Ca + kb_BT * T_BC +
###                    -kf_DT * T_C * Ca + kb_DT * T_CD +
###                    -kf_AT * T_D * Ca + kb_AT * T_AD +
###                    -kf_BT * T_D * Ca + kb_BT * T_BD +
###                    -kf_CT * T_D * Ca + kb_CT * T_CD +
###                    -kf_CT * T_AB * Ca + kb_CT * T_ABC +
###                    -kf_DT * T_AB * Ca + kb_DT * T_ABD +
###                    -kf_BT * T_AC * Ca + kb_BT * T_ABC +
###                    -kf_DT * T_AC * Ca + kb_DT * T_ACD +
###                    -kf_BT * T_AD * Ca + kb_BT * T_ABC +
###                    -kf_CT * T_AD * Ca + kb_CT * T_ACD +
###                    -kf_AT * T_BC * Ca + kb_AT * T_ABC +
###                    -kf_DT * T_BC * Ca + kb_DT * T_BCD +
###                    -kf_AT * T_BC * Ca + kb_AT * T_ABC +
###                    -kf_DT * T_BC * Ca + kb_DT * T_BCD +
###                    -kf_AT * T_BC * Ca + kb_AT * T_ABC +
###                    -kf_DT * T_BC * Ca + kb_DT * T_BCD +
###                    -kf_DT * T_ABC * Ca + kb_DT * T_ABCD +
###                    -kf_CT * T_ABD * Ca + kb_CT * T_ABCD +
###                    -kf_BT * T_ACD * Ca + kb_BT * T_ABCD +
###                    -kf_AT * T_BCD * Ca + kb_AT * T_ABCD +
###                    -kf_AR * R0 * Ca + kb_AR * R_A +        
###                    -kf_BR * R0 * Ca + kb_BR * R_B +
###                    -kf_CR * R0 * Ca + kb_CR * R_C +
###                    -kf_DR * R0 * Ca + kb_DR * R_D +
###                    -kf_BR * R_A * Ca + kb_BR * R_AB +
###                    -kf_CR * R_A * Ca + kb_CR * R_AC +
###                    -kf_DR * R_A * Ca + kb_DR * R_AD +
###                    -kf_AR * R_B * Ca + kb_AR * R_AB +
###                    -kf_CR * R_B * Ca + kb_CR * R_BC +
###                    -kf_DR * R_B * Ca + kb_DR * R_BD +
###                    -kf_AR * R_C * Ca + kb_AR * R_AC +
###                    -kf_BR * R_C * Ca + kb_BR * R_BC +
###                    -kf_DR * R_C * Ca + kb_DR * R_CD +
###                    -kf_AR * R_D * Ca + kb_AR * R_AD +
###                    -kf_BR * R_D * Ca + kb_BR * R_BD +
###                    -kf_CR * R_D * Ca + kb_CR * R_CD +
###                    -kf_CR * R_AB * Ca + kb_CR * R_ABC +
###                    -kf_DR * R_AB * Ca + kb_DR * R_ABD +
###                    -kf_BR * R_AC * Ca + kb_BR * R_ABC +
###                    -kf_DR * R_AC * Ca + kb_DR * R_ACD +
###                    -kf_BR * R_AD * Ca + kb_BR * R_ABC +
###                    -kf_CR * R_AD * Ca + kb_CR * R_ACD +
###                    -kf_AR * R_BC * Ca + kb_AR * R_ABC +
###                    -kf_DR * R_BC * Ca + kb_DR * R_BCD +
###                    -kf_AR * R_BC * Ca + kb_AR * R_ABC +
###                    -kf_DR * R_BC * Ca + kb_DR * R_BCD +
###                    -kf_AR * R_BC * Ca + kb_AR * R_ABC +
###                    -kf_DR * R_BC * Ca + kb_DR * R_BCD +
###                    -kf_DR * R_ABC * Ca + kb_DR * R_ABCD +
###                    -kf_CR * R_ABD * Ca + kb_CR * R_ABCD +
###                    -kf_BR * R_ACD * Ca + kb_BR * R_ABCD +
###                    -kf_AR * R_BCD * Ca + kb_AR * R_ABCD
###
###    end
###
###    @derived begin
###        F_F0 = @. Normal((OGB5 + Fmax_Fmin * CaOGB5) / (OGB5₀ + Fmax_Fmin * CaOGB5₀), σ)
###    end
###
###    @observed begin
###        DMn_s   = DMn_s
###        CaDMn_s = CaDMn_s
###        DMn_f   = DMn_f
###        CaDMn_f = CaDMn_f
###        PP      = PP
###        CaPP    = CaPP
###        OGB5    = OGB5
###        CaOGB5  = CaOGB5
###        T0      = T0
###        T_A     = T_A
###        T_B     = T_B
###        T_C     = T_C
###        T_D     = T_D
###        T_AB    = T_AB
###        T_AC    = T_AC
###        T_AD    = T_AD
###        T_BC    = T_BC
###        T_BD    = T_BD
###        T_CD    = T_CD
###        T_ABC   = T_ABC
###        T_ABD   = T_ABD
###        T_ACD   = T_ACD
###        T_BCD   = T_BCD
###        T_ABCD  = T_ABCD
###        R0      = R0
###        R_A     = R_A
###        R_B     = R_B
###        R_C     = R_C
###        R_D     = R_D
###        R_AB    = R_AB
###        R_AC    = R_AC
###        R_AD    = R_AD
###        R_BC    = R_BC
###        R_BD    = R_BD
###        R_CD    = R_CD
###        R_ABC   = R_ABC
###        R_ABD   = R_ABD
###        R_ACD   = R_ACD
###        R_BCD   = R_BCD
###        R_ABCD  = R_ABCD
###        Ca      = Ca
###    end
###
###end;