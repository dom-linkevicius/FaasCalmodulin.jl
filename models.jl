Blackwell_const_names = (
    :tv_kon_1,
    :tv_Kd_1,
    :tv_kon_2,
    :tv_Kd_2,
#    :tv_mα,
#    :tv_α₀
)
Pepke_M2_const_names = (
    :tv_kon_2N, 
    :tv_Kd_2N,
    :tv_kon_2C, 
    :tv_Kd_2C,
#    :tv_mα,
#    :tv_α₀
);
Faas_const_names = (
    :tv_kon_TN,
    :tv_Kd_TN,
    :tv_kon_RN,
    :tv_Kd_RN,
    :tv_kon_TC,
    :tv_Kd_TC,
    :tv_kon_RC,
    :tv_Kd_RC,
#    :tv_mα,
#    :tv_α₀
);
Byrne_const_names = (
    :tv_k01_N,
    :tv_K01d_N,
    :tv_k02_N,
    :tv_K02d_N,
    :tv_k13_N,
    :tv_K13d_N,
    :tv_k23_N,
###    :tv_K23d_N,
    :tv_k01_C,
    :tv_K01d_C,
    :tv_k02_C,
    :tv_K02d_C,
    :tv_k13_C,
    :tv_K13d_C,
    :tv_k23_C,
###    :tv_K23d_C,
#    :tv_mα,
#    :tv_α₀
)


Blackwell_params = (;
    tv_kon_1 = 3.78,
    tv_Kd_1  = -6.77,
    tv_kon_2 = 5.0,
    tv_Kd_2  = -4.0,
    tv_mα     = 0.0011,
    tv_α₀     = -0.39,
    Ω         = ones(3),
    σ_vec     = [1; 1; 1]
)
Faas_params = (;
    tv_kon_TN = 5.9,  
    tv_Kd_TN  = -3.7,
    tv_kon_RN = 7.5, 
    tv_Kd_RN  = -6.1,
    tv_kon_TC = 4.9, 
    tv_Kd_TC  = -4.6,
    tv_kon_RC = 4.4, 
    tv_Kd_RC  = -6.6,
    tv_mα     = 0.0011,
    tv_α₀     = -0.39,
    Ω         = ones(3),
    σ_vec     = [1; 1; 1]
);
Pepke_M1_params = (;
    tv_kon_TN = log10(142.5e3),  
    tv_Kd_TN  = log10(27.5e-6),
    tv_kon_RN = log10(175e3), 
    tv_Kd_RN  = log10(6.15e-6),
    tv_kon_TC = log10(5.4e3), 
    tv_Kd_TC  = log10(9.65e-6),
    tv_kon_RC = log10(15e3), 
    tv_Kd_RC  = log10(1.05e-6),
    tv_mα     = 0.0011,
    tv_α₀     = -0.39,
    Ω         = ones(3),
    σ_vec     = [1; 1; 1]
);
Pepke_M2_params = (;
    tv_kon_2N = log10(175e3), 
    tv_Kd_2N  = log10(6.15e-6),
    tv_kon_2C = log10(15e3), 
    tv_Kd_2C  = log10(1.05e-6),
    tv_mα     = 0.0011,
    tv_α₀     = -0.39,
    Ω         = ones(3),
    σ_vec     = [1; 1; 1]
);
Byrne_params = (;
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
    tv_mα       = 0.0011,
    tv_α₀       = -0.39,
    Ω           = ones(3),
    σ_vec       = [1; 1; 1]
)


Blackwell_scheme = @model begin ### Simplest scheme
    @param begin
        tv_kon_1     ∈ Uniform(2,  12)
        tv_Kd_1      ∈ Uniform(-12, -2)
        tv_kon_2     ∈ Uniform(2,  12)
        tv_Kd_2      ∈ Uniform(-12, -2)

        tv_mα         ∈ RealDomain(; lower=1e-10, init= 0.0011, upper=1)
        tv_α₀         ∈ RealDomain(;              init = -0.39)

        Ω             ∈ VectorDomain(3; lower=1e-10, upper=3, init=1)
        σ_vec         ∈ VectorDomain(3; lower=1e-10, upper=3, init=1)

    end

    @random begin
        η     ~ MvNormal(I(3) .* Ω)
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
        cov_τ_decay
        cov_Ca_i
        CaM_equil
        cov_Fluo4FF_t
    end



    @pre begin
        σ         = σ_vec[batch_id]
        ###U         = max(0.0, (tv_mα*PCD + tv_α₀) * (1 + η[batch_id]))
        U         = (tv_mα*PCD + tv_α₀) * (1 + η[batch_id])

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
        Fluo4FF_t = cov_Fluo4FF_t
        Ca_free = cov_Ca_free

        τ_decay = cov_τ_decay
        Ca_i    = cov_Ca_i

        kon_D    = 8.77e5         ### in M^-1 ms^-1
        koff_D   = 33.2          ### in ms^-1 
        Kd_D     = koff_D / kon_D
    
        Kd_Fluo4FF = 23e-6 ## https://pubmed.ncbi.nlm.nih.gov/29604967/
        kon_Fluo4FF  = 1 ## assumed it, it doesn't matter currently
        koff_Fluo4FF = Kd_Fluo4FF * kon_Fluo4FF 

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
        Fluo4FF₀    = (Kd_Fluo4FF * Fluo4FF_t) / (Ca_free + Kd_Fluo4FF)
        CaFluo4FF₀  = Fluo4FF_t - Fluo4FF₀       

        CaM₀_all = solve_Bwell_eq(Kd_1, Kd_2, Ca_free, CaM_t, CaM_equil)

        DMn_s₀_u   = (1-f_frac) * DMn₀   * U# * 0
        CaDMn_s₀_u = (1-f_frac) * CaDMn₀ * U# * 0
        DMn_f₀_u   =    f_frac  * DMn₀   * U# * 0
        CaDMn_f₀_u =    f_frac  * CaDMn₀ * U# * 0
    end

    @init begin
        Ca       = Ca_free
        DMn_s    = DMn_s₀_u 
        CaDMn_s  = CaDMn_s₀_u
        DMn_f    = DMn_f₀_u
        CaDMn_f  = CaDMn_f₀_u
        PP       = 0.0
        CaPP     = 0.0
        OGB5     = OGB5₀
        CaOGB5   = CaOGB5₀
        Fluo4FF  = Fluo4FF₀
        CaFluo4FF= CaFluo4FF₀
        CaM      = CaM₀_all[1]
        CaM2Ca   = CaM₀_all[2]
        CaM4Ca   = CaM₀_all[3]
    end

    @dynamics begin
        Ca'       = -kon_DMn * DMn_s * Ca + koff_DMn * CaDMn_s +
                    -kon_DMn * DMn_f * Ca + koff_DMn * CaDMn_f +
                    -kon_DMn * PP    * Ca + koff_PP  * CaPP    + 
                    -kon_D   * OGB5  * Ca + koff_D   * CaOGB5  +
                    -kon_Fluo4FF * Fluo4FF * Ca + koff_Fluo4FF * CaFluo4FF + 
                    -2kon_1   * CaM   * Ca^2 + 2koff_1 * CaM2Ca +
                    -2kon_2   * CaM2Ca* Ca^2 + 2koff_2 * CaM4Ca +
                    -(Ca - Ca_i) * τ_decay
        DMn_s'    = -kon_DMn * DMn_s * Ca + koff_DMn * CaDMn_s -   DMn_s * τ_s
        CaDMn_s'  =  kon_DMn * DMn_s * Ca - koff_DMn * CaDMn_s - CaDMn_s * τ_s
        DMn_f'    = -kon_DMn * DMn_f * Ca + koff_DMn * CaDMn_f -   DMn_f * τ_f
        CaDMn_f'  =  kon_DMn * DMn_f * Ca - koff_DMn * CaDMn_f - CaDMn_f * τ_f
        PP'       = -kon_PP  * PP    * Ca + koff_PP  * CaPP    +  2DMn_s * τ_s +  2DMn_f * τ_f + CaDMn_s * τ_s + CaDMn_f * τ_f
        CaPP'     =  kon_PP  * PP    * Ca - koff_PP  * CaPP    + CaDMn_s * τ_s + CaDMn_f * τ_f
        OGB5'     = -kon_D   * OGB5  * Ca + koff_D   * CaOGB5
        CaOGB5'   =  kon_D   * OGB5  * Ca - koff_D   * CaOGB5
        Fluo4FF'  = -kon_Fluo4FF   * Fluo4FF  * Ca + koff_Fluo4FF   * CaFluo4FF
        CaFluo4FF'=  kon_Fluo4FF   * Fluo4FF  * Ca - koff_Fluo4FF   * CaFluo4FF
        CaM'      = -kon_1 * CaM * Ca^2 + koff_1 * CaM2Ca 
        CaM2Ca'   =  kon_1 * CaM * Ca^2 - koff_1 * CaM2Ca -kon_2 * CaM2Ca * Ca^2 + koff_2 * CaM4Ca 
        CaM4Ca'   =                                        kon_2 * CaM2Ca * Ca^2 - koff_2 * CaM4Ca
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
        Fluo4FF    = Fluo4FF
        CaFluo4FF  = CaFluo4FF
        CaM     = CaM
        CaM2Ca  = CaM2Ca
        CaM4Ca  = CaM4Ca
        Ca      = Ca
    end
end;


Pepke_m2_scheme = @model begin ### Simplest scheme with unique lobes
    @param begin

        tv_kon_2C     ∈ Uniform(2,  12)
        tv_Kd_2C      ∈ Uniform(-12, -2)
        tv_kon_2N     ∈ Uniform(2,  12)
        tv_Kd_2N      ∈ Uniform(-12, -2)

        tv_mα         ∈ RealDomain(; lower=1e-10, init= 0.0011, upper=1)
        tv_α₀         ∈ RealDomain(;              init = -0.39)

        Ω             ∈ VectorDomain(3; lower=1e-10, upper=3, init=1)
        σ_vec         ∈ VectorDomain(3; lower=1e-10, upper=3, init=1)
    end

    @random begin
        η     ~ MvNormal(I(3) .* Ω)
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
        cov_τ_decay
        cov_Ca_i
        CaM_equil
        cov_Fluo4FF_t
    end



    @pre begin
        σ         = σ_vec[batch_id]
        ###U         = max(0.0, (tv_mα*PCD + tv_α₀) * (1 + η[batch_id]))
        U         = (tv_mα*PCD + tv_α₀) * (1 + η[batch_id])

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
        Fluo4FF_t = cov_Fluo4FF_t
        Ca_free = cov_Ca_free

        τ_decay = cov_τ_decay
        Ca_i    = cov_Ca_i

        kon_D    = 8.77e5         ### in M^-1 ms^-1
        koff_D   = 33.2          ### in ms^-1 
        Kd_D     = koff_D / kon_D

        Kd_Fluo4FF = 23e-6 ## https://pubmed.ncbi.nlm.nih.gov/29604967/
        kon_Fluo4FF  = 1 ## assumed it, it doesn't matter currently
        koff_Fluo4FF = Kd_Fluo4FF * kon_Fluo4FF 

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
        Fluo4FF₀    = (Kd_Fluo4FF * Fluo4FF_t) / (Ca_free + Kd_Fluo4FF)
        CaFluo4FF₀  = Fluo4FF_t - Fluo4FF₀    

        c0       = solve_Blackwell_CN_eqs(kon_2C, koff_2C, kon_2N, koff_2N, Ca_free, CaM_t, CaM_equil)

        DMn_s₀_u   = (1-f_frac) * DMn₀   * U# * 0
        CaDMn_s₀_u = (1-f_frac) * CaDMn₀ * U# * 0
        DMn_f₀_u   =    f_frac  * DMn₀   * U# * 0
        CaDMn_f₀_u =    f_frac  * CaDMn₀ * U# * 0
    end

    @init begin
        Ca       = Ca_free
        DMn_s    = DMn_s₀_u 
        CaDMn_s  = CaDMn_s₀_u
        DMn_f    = DMn_f₀_u
        CaDMn_f  = CaDMn_f₀_u
        PP       = 0.0
        CaPP     = 0.0
        OGB5     = OGB5₀
        CaOGB5   = CaOGB5₀
        Fluo4FF  = Fluo4FF₀
        CaFluo4FF= CaFluo4FF₀
        CaM0     = c0[1]
        CaM2C    = c0[2]
        CaM2N    = c0[3]
        CaM4     = c0[4]
    end

    @dynamics begin
        Ca'       = -kon_DMn * DMn_s * Ca + koff_DMn * CaDMn_s +
                    -kon_DMn * DMn_f * Ca + koff_DMn * CaDMn_f +
                    -kon_DMn * PP    * Ca + koff_PP  * CaPP    + 
                    -kon_D   * OGB5  * Ca + koff_D   * CaOGB5  +
                    -kon_Fluo4FF * Fluo4FF * Ca + koff_Fluo4FF * CaFluo4FF +
                    -kon_2C * CaM0  * Ca^2 + koff_2C * CaM2C +
                    -kon_2N * CaM0  * Ca^2 + koff_2N * CaM2N +
                    -kon_2N * CaM2C * Ca^2 + koff_2N * CaM4 +
                    -kon_2C * CaM2N * Ca^2 + koff_2C * CaM4 +
                    -(Ca - Ca_i) * τ_decay
        DMn_s'    = -kon_DMn * DMn_s * Ca + koff_DMn * CaDMn_s -   DMn_s * τ_s
        CaDMn_s'  =  kon_DMn * DMn_s * Ca - koff_DMn * CaDMn_s - CaDMn_s * τ_s
        DMn_f'    = -kon_DMn * DMn_f * Ca + koff_DMn * CaDMn_f -   DMn_f * τ_f
        CaDMn_f'  =  kon_DMn * DMn_f * Ca - koff_DMn * CaDMn_f - CaDMn_f * τ_f
        PP'       = -kon_PP  * PP    * Ca + koff_PP  * CaPP    +  2DMn_s * τ_s +  2DMn_f * τ_f + CaDMn_s * τ_s + CaDMn_f * τ_f
        CaPP'     =  kon_PP  * PP    * Ca - koff_PP  * CaPP    + CaDMn_s * τ_s + CaDMn_f * τ_f
        OGB5'     = -kon_D   * OGB5  * Ca + koff_D   * CaOGB5
        CaOGB5'   =  kon_D   * OGB5  * Ca - koff_D   * CaOGB5
        Fluo4FF'  = -kon_Fluo4FF   * Fluo4FF  * Ca + koff_Fluo4FF   * CaFluo4FF
        CaFluo4FF'=  kon_Fluo4FF   * Fluo4FF  * Ca - koff_Fluo4FF   * CaFluo4FF
        CaM0'     = -kon_2C * CaM0  * Ca^2 + koff_2C * CaM2C +
                    -kon_2N * CaM0  * Ca^2 + koff_2N * CaM2N 
        CaM2C'    =  kon_2C * CaM0  * Ca^2 - koff_2C * CaM2C +
                    -kon_2N * CaM2C * Ca^2 + koff_2N * CaM4 
        CaM2N'    =  kon_2N * CaM0  * Ca^2 - koff_2N * CaM2N +
                    -kon_2C * CaM2N * Ca^2 + koff_2C * CaM4 
        CaM4'     =  kon_2N * CaM2C * Ca^2 - koff_2N * CaM4 +
                     kon_2C * CaM2N * Ca^2 - koff_2C * CaM4 
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
        Fluo4FF    = Fluo4FF
        CaFluo4FF  = CaFluo4FF
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

        Ω             ∈ VectorDomain(3; lower=1e-10, upper=3, init=1)
        σ_vec         ∈ VectorDomain(3; lower=1e-10, upper=3, init=1)

    end

    @random begin
        η     ~ MvNormal(I(3) .* Ω)
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
        cov_τ_decay
        cov_Ca_i
        CaM_equil
        cov_Fluo4FF_t
    end



    @pre begin
        σ         = σ_vec[batch_id]
        ###U         = max(0.0, (tv_mα*PCD + tv_α₀) * (1 + η[batch_id]))
        U         = (tv_mα*PCD + tv_α₀) * (1 + η[batch_id])

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
        Fluo4FF_t = cov_Fluo4FF_t
        Ca_free = cov_Ca_free

        τ_decay = cov_τ_decay
        Ca_i    = cov_Ca_i

        kon_D    = 8.77e5         ### in M^-1 ms^-1
        koff_D   = 33.2          ### in ms^-1 
        Kd_D     = koff_D / kon_D

        Kd_Fluo4FF = 23e-6 ## https://pubmed.ncbi.nlm.nih.gov/29604967/
        kon_Fluo4FF  = 1 ## assumed it, it doesn't matter currently
        koff_Fluo4FF = Kd_Fluo4FF * kon_Fluo4FF 

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
        Fluo4FF₀    = (Kd_Fluo4FF * Fluo4FF_t) / (Ca_free + Kd_Fluo4FF)
        CaFluo4FF₀  = Fluo4FF_t - Fluo4FF₀

        c0       = solve_Blackwell_TR_eqs(kon_1_T, koff_1_T, kon_2_T, koff_2_T, 
                                          kon_1_R, koff_1_R, kon_2_R, koff_2_R, 
                                          k_TR_0, k_RT_0, k_TR_2, k_RT_2, k_TR_4, k_RT_4,
                                          Ca_free, CaM_t, CaM_equil)

        DMn_s₀_u   = (1-f_frac) * DMn₀   * U# * 0
        CaDMn_s₀_u = (1-f_frac) * CaDMn₀ * U# * 0
        DMn_f₀_u   =    f_frac  * DMn₀   * U# * 0
        CaDMn_f₀_u =    f_frac  * CaDMn₀ * U# * 0
    end

    @init begin
        Ca       = Ca_free
        DMn_s    = DMn_s₀_u 
        CaDMn_s  = CaDMn_s₀_u
        DMn_f    = DMn_f₀_u
        CaDMn_f  = CaDMn_f₀_u
        PP       = 0.0
        CaPP     = 0.0
        OGB5     = OGB5₀
        CaOGB5   = CaOGB5₀
        Fluo4FF  = Fluo4FF₀
        CaFluo4FF= CaFluo4FF₀
        CaM0_R   = c0[1]
        CaM2_R   = c0[2]
        CaM4_R   = c0[3]
        CaM0_T   = c0[4]
        CaM2_T   = c0[5]
        CaM4_T   = c0[6]
    end

    @dynamics begin
        Ca'       = -kon_DMn * DMn_s * Ca + koff_DMn * CaDMn_s +
                    -kon_DMn * DMn_f * Ca + koff_DMn * CaDMn_f +
                    -kon_DMn * PP    * Ca + koff_PP  * CaPP    + 
                    -kon_D   * OGB5  * Ca + koff_D   * CaOGB5  +
                    -kon_Fluo4FF * Fluo4FF * Ca + koff_Fluo4FF * CaFluo4FF + 
                    -kon_1_R * CaM0_R * Ca^2 + koff_1_R * CaM2_R + 
                    -kon_2_R * CaM2_R * Ca^2 + koff_2_R * CaM4_R +
                    -kon_1_T * CaM0_T * Ca^2 + koff_1_T * CaM2_T +
                    -kon_2_T * CaM2_T * Ca^2 + koff_2_T * CaM4_T +
                    -(Ca - Ca_i) * τ_decay
        DMn_s'    = -kon_DMn * DMn_s * Ca + koff_DMn * CaDMn_s -   DMn_s * τ_s
        CaDMn_s'  =  kon_DMn * DMn_s * Ca - koff_DMn * CaDMn_s - CaDMn_s * τ_s
        DMn_f'    = -kon_DMn * DMn_f * Ca + koff_DMn * CaDMn_f -   DMn_f * τ_f
        CaDMn_f'  =  kon_DMn * DMn_f * Ca - koff_DMn * CaDMn_f - CaDMn_f * τ_f
        PP'       = -kon_PP  * PP    * Ca + koff_PP  * CaPP    +  2DMn_s * τ_s +  2DMn_f * τ_f + CaDMn_s * τ_s + CaDMn_f * τ_f
        CaPP'     =  kon_PP  * PP    * Ca - koff_PP  * CaPP    + CaDMn_s * τ_s + CaDMn_f * τ_f
        OGB5'     = -kon_D   * OGB5  * Ca + koff_D   * CaOGB5
        CaOGB5'   =  kon_D   * OGB5  * Ca - koff_D   * CaOGB5
        Fluo4FF'  = -kon_Fluo4FF   * Fluo4FF  * Ca + koff_Fluo4FF   * CaFluo4FF
        CaFluo4FF'=  kon_Fluo4FF   * Fluo4FF  * Ca - koff_Fluo4FF   * CaFluo4FF
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
        Fluo4FF    = Fluo4FF
        CaFluo4FF  = CaFluo4FF
        CaM0_R  = CaM0_R
        CaM2_R  = CaM2_R
        CaM4_R  = CaM4_R
        CaM0_T  = CaM0_T
        CaM2_T  = CaM2_T
        CaM4_T  = CaM4_T
        Ca      = Ca
    end
end;



Faas_scheme = @model begin ### Faas scheme
    @param begin
        tv_kon_TN     ∈ Uniform(2, 12)
        tv_Kd_TN      ∈ Uniform(-12, -2) 
        tv_kon_RN     ∈ Uniform(2, 12)
        tv_Kd_RN      ∈ Uniform(-12, -2)

        tv_kon_TC     ∈ Uniform(2, 12)
        tv_Kd_TC      ∈ Uniform(-12, -2)
        tv_kon_RC     ∈ Uniform(2, 12)
        tv_Kd_RC      ∈ Uniform(-12, -2)

        tv_mα         ∈ RealDomain(; lower=1e-10, init= 0.0011, upper=1)
        tv_α₀         ∈ RealDomain(;              init = -0.39)

        Ω             ∈ VectorDomain(3; lower=1e-10, upper=3, init=1)
        σ_vec         ∈ VectorDomain(3; lower=1e-10, upper=3, init=1)
    end

    @random begin
        η     ~ MvNormal(I(3) .* Ω)
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
        cov_τ_decay
        cov_Ca_i
        CaM_equil
        cov_Fluo4FF_t
    end



    @pre begin
        σ         = σ_vec[batch_id]
        ###U         = max(0.0, (tv_mα*PCD + tv_α₀) * (1 + η[batch_id]))
        U         = (tv_mα*PCD + tv_α₀) * (1 + η[batch_id])

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
        Fluo4FF_t = cov_Fluo4FF_t
        Ca_free = cov_Ca_free

        τ_decay = cov_τ_decay
        Ca_i    = cov_Ca_i

        kon_D    = 8.77e5         ### in M^-1 ms^-1
        koff_D   = 33.2          ### in ms^-1 
        Kd_D     = koff_D / kon_D

        Kd_Fluo4FF = 23e-6 ## https://pubmed.ncbi.nlm.nih.gov/29604967/
        kon_Fluo4FF  = 1 ## assumed it, it doesn't matter currently
        koff_Fluo4FF = Kd_Fluo4FF * kon_Fluo4FF  

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
        Fluo4FF₀    = (Kd_Fluo4FF * Fluo4FF_t) / (Ca_free + Kd_Fluo4FF)
        CaFluo4FF₀  = Fluo4FF_t - Fluo4FF₀  

        N₀_all = get_Faas_lobe_eqs(Ca_free, Kd_TN, Kd_RN, CaM_t, CaM_equil)
        C₀_all = get_Faas_lobe_eqs(Ca_free, Kd_TC, Kd_RC, CaM_t, CaM_equil)

        DMn_s₀_u   = (1-f_frac) * DMn₀   * U
        CaDMn_s₀_u = (1-f_frac) * CaDMn₀ * U
        DMn_f₀_u   =    f_frac  * DMn₀   * U
        CaDMn_f₀_u =    f_frac  * CaDMn₀ * U
    end

    @init begin
        Ca       = Ca_free
        DMn_s    = DMn_s₀_u 
        CaDMn_s  = CaDMn_s₀_u
        DMn_f    = DMn_f₀_u
        CaDMn_f  = CaDMn_f₀_u
        PP       = 0.0
        CaPP     = 0.0
        OGB5     = OGB5₀
        CaOGB5   = CaOGB5₀
        Fluo4FF  = Fluo4FF₀
        CaFluo4FF= CaFluo4FF₀
        NtNt     = N₀_all[1]
        NtNr     = N₀_all[2]
        NrNr     = N₀_all[3]
        CtCt     = C₀_all[1]
        CtCr     = C₀_all[2]
        CrCr     = C₀_all[3]
    end

    @dynamics begin
        Ca'       = -kon_DMn * DMn_s * Ca + koff_DMn * CaDMn_s +
                    -kon_DMn * DMn_f * Ca + koff_DMn * CaDMn_f +
                    -kon_DMn * PP    * Ca + koff_PP  * CaPP    + 
                    -kon_D   * OGB5  * Ca + koff_D   * CaOGB5  +
                    -kon_Fluo4FF * Fluo4FF * Ca + koff_Fluo4FF * CaFluo4FF + 
                    -2kon_TN * NtNt  * Ca + koff_TN  * NtNr    +
                    -2kon_TC * CtCt  * Ca + koff_TC  * CtCr    +     
                    -kon_RN  * NtNr  * Ca + 2koff_RN * NrNr    +
                    -kon_RC  * CtCr  * Ca + 2koff_RC * CrCr    +
                    -(Ca - Ca_i) * τ_decay
        DMn_s'    = -kon_DMn * DMn_s * Ca + koff_DMn * CaDMn_s -   DMn_s * τ_s
        CaDMn_s'  =  kon_DMn * DMn_s * Ca - koff_DMn * CaDMn_s - CaDMn_s * τ_s
        DMn_f'    = -kon_DMn * DMn_f * Ca + koff_DMn * CaDMn_f -   DMn_f * τ_f
        CaDMn_f'  =  kon_DMn * DMn_f * Ca - koff_DMn * CaDMn_f - CaDMn_f * τ_f
        PP'       = -kon_PP  * PP    * Ca + koff_PP  * CaPP    +  2DMn_s * τ_s +  2DMn_f * τ_f + CaDMn_s * τ_s + CaDMn_f * τ_f
        CaPP'     =  kon_PP  * PP    * Ca - koff_PP  * CaPP    + CaDMn_s * τ_s + CaDMn_f * τ_f
        OGB5'     = -kon_D   * OGB5  * Ca + koff_D   * CaOGB5
        CaOGB5'   =  kon_D   * OGB5  * Ca - koff_D   * CaOGB5
        Fluo4FF'  = -kon_Fluo4FF   * Fluo4FF  * Ca + koff_Fluo4FF   * CaFluo4FF
        CaFluo4FF'=  kon_Fluo4FF   * Fluo4FF  * Ca - koff_Fluo4FF   * CaFluo4FF
        NtNt'     = -2kon_TN * NtNt  * Ca + koff_TN  * NtNr
        NtNr'     =  2kon_TN * NtNt  * Ca - koff_TN  * NtNr - kon_RN * NtNr * Ca + 2koff_RN  * NrNr
        NrNr'     =                                         + kon_RN * NtNr * Ca - 2koff_RN  * NrNr
        CtCt'     = -2kon_TC * CtCt  * Ca + koff_TC  * CtCr
        CtCr'     =  2kon_TC * CtCt  * Ca - koff_TC  * CtCr - kon_RC * CtCr * Ca + 2koff_RC  * CrCr
        CrCr'     =                                         + kon_RC * CtCr * Ca - 2koff_RC  * CrCr
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
        Fluo4FF    = Fluo4FF
        CaFluo4FF  = CaFluo4FF
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
###        tv_K23d_N     ∈ Uniform(-12, -2) ### not a free param due to microscopic reversibility

        tv_k01_C      ∈ Uniform(2,   12)
        tv_K01d_C     ∈ Uniform(-12, -2)
        tv_k02_C      ∈ Uniform(2,   12)
        tv_K02d_C     ∈ Uniform(-12, -2)
        tv_k13_C      ∈ Uniform(2,   12)
        tv_K13d_C     ∈ Uniform(-12, -2)
        tv_k23_C      ∈ Uniform(2,   12)
###        tv_K23d_C     ∈ Uniform(-12, -2) ### not a free param due to microscopic reversibility

        tv_mα         ∈ RealDomain(; lower=1e-10, init= 0.001, upper=1)
        tv_α₀         ∈ RealDomain(;              init = -0.39)

        Ω             ∈ VectorDomain(3; lower=1e-10, upper=3, init=1)
        σ_vec         ∈ VectorDomain(3; lower=1e-10, upper=3, init=1)
    end

    @random begin
        η     ~ MvNormal(I(3) .* Ω)
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
        cov_τ_decay
        cov_Ca_i
        CaM_equil
        cov_Fluo4FF_t
    end



    @pre begin
        σ         = σ_vec[batch_id]
        ###U         = max(0.0, (tv_mα*PCD + tv_α₀) * (1 + η[batch_id]))
        U         = (tv_mα*PCD + tv_α₀) * (1 + η[batch_id])

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
        Fluo4FF_t = cov_Fluo4FF_t
        Ca_free = cov_Ca_free

        τ_decay = cov_τ_decay
        Ca_i    = cov_Ca_i

        kon_D    = 8.77e5         ### in M^-1 ms^-1
        koff_D   = 33.2          ### in ms^-1 
        Kd_D     = koff_D / kon_D

        Kd_Fluo4FF = 23e-6 ## https://pubmed.ncbi.nlm.nih.gov/29604967/
        kon_Fluo4FF  = 1 ## assumed it, it doesn't matter currently
        koff_Fluo4FF = Kd_Fluo4FF * kon_Fluo4FF 

        K01d_N  = 10^tv_K01d_N 
        k01_N   = 10^tv_k01_N
        k10_N   = K01d_N * k01_N 

        K02d_N  = 10^tv_K02d_N 
        k02_N   = 10^tv_k02_N
        k20_N   = K02d_N * k02_N 

        K13d_N  = 10^tv_K13d_N 
        k13_N   = 10^tv_k13_N
        k31_N   = K13d_N * k13_N 

        K23d_N  = K01d_N * K13d_N / K02d_N  ###10^tv_K23d_N ##from microscopic reversibility
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

        K23d_C  = K01d_C * K13d_C / K02d_C ###10^tv_K23d_C ##from microscopic reversibility
        k23_C   = 10^tv_k23_C
        k32_C   = K23d_C * k23_C 


        DMn₀     = (Kd_DMn * DMn_t)  / (Ca_free + Kd_DMn)
        CaDMn₀   = DMn_t - DMn₀
        OGB5₀    = (Kd_D   * OGB5_t) / (Ca_free + Kd_D)
        CaOGB5₀  = OGB5_t - OGB5₀
        Fluo4FF₀    = (Kd_Fluo4FF * Fluo4FF_t) / (Ca_free + Kd_Fluo4FF)
        CaFluo4FF₀  = Fluo4FF_t - Fluo4FF₀     


        N₀ = solve_Byrne_eq(k01_N, k10_N, k02_N, k20_N, k13_N, k31_N, k23_N, k32_N, Ca_free, CaM_t, CaM_equil)
        C₀ = solve_Byrne_eq(k01_C, k10_C, k02_C, k20_C, k13_C, k31_C, k23_C, k32_C, Ca_free, CaM_t, CaM_equil)
        

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
        Fluo4FF  = Fluo4FF₀
        CaFluo4FF= CaFluo4FF₀
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
        Ca'       = -kon_DMn * DMn_s * Ca + koff_DMn * CaDMn_s +
                    -kon_DMn * DMn_f * Ca + koff_DMn * CaDMn_f +
                    -kon_DMn * PP    * Ca + koff_PP  * CaPP    + 
                    -kon_D   * OGB5  * Ca + koff_D   * CaOGB5  +
                    -kon_Fluo4FF * Fluo4FF * Ca + koff_Fluo4FF * CaFluo4FF + 
                    -k01_N * Ca * N0 + k10_N * N1 +
                    -k02_N * Ca * N0 + k20_N * N2 +
                    -k13_N * Ca * N1 + k31_N * N3 +
                    -k23_N * Ca * N2 + k32_N * N3 +
                    -k01_C * Ca * C0 + k10_C * C1 +
                    -k02_C * Ca * C0 + k20_C * C2 +
                    -k13_C * Ca * C1 + k31_C * C3 + 
                    -k23_C * Ca * C2 + k32_C * C3 + 
                    -(Ca - Ca_i) * τ_decay
        DMn_s'    = -kon_DMn * DMn_s * Ca + koff_DMn * CaDMn_s -   DMn_s * τ_s
        CaDMn_s'  =  kon_DMn * DMn_s * Ca - koff_DMn * CaDMn_s - CaDMn_s * τ_s
        DMn_f'    = -kon_DMn * DMn_f * Ca + koff_DMn * CaDMn_f -   DMn_f * τ_f
        CaDMn_f'  =  kon_DMn * DMn_f * Ca - koff_DMn * CaDMn_f - CaDMn_f * τ_f
        PP'       = -kon_PP  * PP    * Ca + koff_PP  * CaPP    +  2DMn_s * τ_s +  2DMn_f * τ_f + CaDMn_s * τ_s + CaDMn_f * τ_f
        CaPP'     =  kon_PP  * PP    * Ca - koff_PP  * CaPP    + CaDMn_s * τ_s + CaDMn_f * τ_f
        OGB5'     = -kon_D   * OGB5  * Ca + koff_D   * CaOGB5
        CaOGB5'   =  kon_D   * OGB5  * Ca - koff_D   * CaOGB5
        Fluo4FF'  = -kon_Fluo4FF   * Fluo4FF  * Ca + koff_Fluo4FF   * CaFluo4FF
        CaFluo4FF'=  kon_Fluo4FF   * Fluo4FF  * Ca - koff_Fluo4FF   * CaFluo4FF
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
        Fluo4FF    = Fluo4FF
        CaFluo4FF  = CaFluo4FF
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