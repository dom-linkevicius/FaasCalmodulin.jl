Blackwell_scheme = @model begin 
    @param begin
        tv_kon_1     ∈ Uniform(2,   9)
        tv_Kd_1      ∈ Uniform(-9, -4)
        tv_kon_2     ∈ Uniform(2,   9)    ### THIS IS N-LOBE IN KIM ET AL
        tv_Kd_2      ∈ Uniform(-9, -4)    ### THIS IS N-LOBE IN KIM ET AL

        μ             ∈ RealDomain(init=-1, lower=-5, upper=5)
        Ω             ∈ Constrained(Normal(0, 1), lower=1) 
        σ_faas        ∈ Constrained(Normal(0, 1), lower=1e-10, upper=1)
        σ_shifman     ∈ Constrained(Normal(0, 1), lower=1e-10, upper=1)
    end

    @random begin
        η     ~ Normal(μ, Ω)
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
        cov_τ_decay
        cov_Ca_i
        CaM_equil
        cov_Fluo4FF_t
        cov_kon_D
        cov_koff_D
        cov_Kd_Fluo4FF
        cov_kon_Fluo4FF
    end



    @pre begin
        U         = sigmoid(η)

        Fmax_Fmin = 39.364

        f_frac   = cov_f_frac 
        τ_f      = 1/cov_τ_f
        τ_s      = 1/cov_τ_s
        kon_DMn  = cov_kon_DMn
        koff_DMn = cov_koff_DMn
        Kd_DMn   = cov_Kd_DMn
        kon_PP   = cov_kon_DMn
        koff_PP  = cov_koff_PP

        DMn_t     = cov_DMn_t
        CaM_t     = cov_CaM_t
        OGB5_t    = cov_OGB5_t
        Fluo4FF_t = cov_Fluo4FF_t
        Ca_free   = cov_Ca_free

        τ_decay = cov_τ_decay
        Ca_i    = cov_Ca_i

        kon_D    = cov_kon_D         
        koff_D   = cov_koff_D          
        Kd_D     = koff_D / kon_D
    
        Kd_Fluo4FF = cov_Kd_Fluo4FF
        kon_Fluo4FF  = cov_kon_Fluo4FF
        koff_Fluo4FF = Kd_Fluo4FF * kon_Fluo4FF 

        Kd_1  = 10^tv_Kd_1
        kon_1 = 10^tv_kon_1
        koff_1= Kd_1 * kon_1

        Kd_2  = 10^tv_Kd_2 
        kon_2 = 10^tv_kon_2
        koff_2= Kd_2 * kon_2 

        DMn₀        = (Kd_DMn * DMn_t)  / (Ca_free + Kd_DMn)
        CaDMn₀      = DMn_t - DMn₀
        OGB5₀       = (Kd_D   * OGB5_t) / (Ca_free + Kd_D)
        CaOGB5₀     = OGB5_t - OGB5₀
        Fluo4FF₀    = (Kd_Fluo4FF * Fluo4FF_t) / (Ca_free + Kd_Fluo4FF)
        CaFluo4FF₀  = Fluo4FF_t - Fluo4FF₀       

        CaM₀_all = solve_Bwell_eq(Kd_1, Kd_2, Ca_free, CaM_t, CaM_equil)

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
        F_F0   ~ @. Normal((OGB5 + Fmax_Fmin * CaOGB5) / (OGB5₀ + Fmax_Fmin * CaOGB5₀), σ_faas)
        Ca_CaM ~ @. Normal((2CaM2Ca + 4CaM4Ca) / CaM_t, (2CaM2Ca + 4CaM4Ca) / CaM_t * σ_shifman)
    end
end;



Bhalla_scheme = @model begin 
    @param begin
        tv_kon_1     ∈ Uniform(2,   9)
        tv_Kd_1      ∈ Uniform(-9, -4)
        tv_kon_2     ∈ Uniform(2,   9)    
        tv_Kd_2      ∈ Uniform(-9, -4)    
        tv_kon_3     ∈ Uniform(2,   9)    
        tv_Kd_3      ∈ Uniform(-9, -4)    

        μ             ∈ RealDomain(init=-1, lower=-5, upper=5)
        Ω             ∈ Constrained(Normal(0, 1), lower=1) 
        σ_faas        ∈ Constrained(Normal(0, 1), lower=1e-10, upper=1)
        σ_shifman     ∈ Constrained(Normal(0, 1), lower=1e-10, upper=1)
    end

    @random begin
        η     ~ Normal(μ, Ω)
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
        cov_τ_decay
        cov_Ca_i
        CaM_equil
        cov_Fluo4FF_t
        cov_kon_D
        cov_koff_D
        cov_Kd_Fluo4FF
        cov_kon_Fluo4FF
    end



    @pre begin
        U         = sigmoid(η)

        Fmax_Fmin = 39.364

        f_frac   = cov_f_frac 
        τ_f      = 1/cov_τ_f
        τ_s      = 1/cov_τ_s
        kon_DMn  = cov_kon_DMn
        koff_DMn = cov_koff_DMn
        Kd_DMn   = cov_Kd_DMn
        kon_PP   = cov_kon_DMn
        koff_PP  = cov_koff_PP

        DMn_t     = cov_DMn_t
        CaM_t     = cov_CaM_t
        OGB5_t    = cov_OGB5_t
        Fluo4FF_t = cov_Fluo4FF_t
        Ca_free   = cov_Ca_free

        τ_decay = cov_τ_decay
        Ca_i    = cov_Ca_i

        kon_D    = cov_kon_D         
        koff_D   = cov_koff_D          
        Kd_D     = koff_D / kon_D
    
        Kd_Fluo4FF = cov_Kd_Fluo4FF
        kon_Fluo4FF  = cov_kon_Fluo4FF
        koff_Fluo4FF = Kd_Fluo4FF * kon_Fluo4FF 

        Kd_1  = 10^tv_Kd_1
        kon_1 = 10^tv_kon_1
        koff_1= Kd_1 * kon_1

        Kd_2  = 10^tv_Kd_2 
        kon_2 = 10^tv_kon_2
        koff_2= Kd_2 * kon_2 

        Kd_3  = 10^tv_Kd_3 
        kon_3 = 10^tv_kon_3
        koff_3= Kd_3 * kon_3 

        DMn₀        = (Kd_DMn * DMn_t)  / (Ca_free + Kd_DMn)
        CaDMn₀      = DMn_t - DMn₀
        OGB5₀       = (Kd_D   * OGB5_t) / (Ca_free + Kd_D)
        CaOGB5₀     = OGB5_t - OGB5₀
        Fluo4FF₀    = (Kd_Fluo4FF * Fluo4FF_t) / (Ca_free + Kd_Fluo4FF)
        CaFluo4FF₀  = Fluo4FF_t - Fluo4FF₀       

        CaM₀_all = solve_Bhalla_eq(kon_1, koff_1, kon_2, koff_2, kon_3, koff_3, Ca_free, CaM_t)

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
        CaM      = CaM₀_all[1]
        CaM2Ca   = CaM₀_all[2]
        CaM3Ca   = CaM₀_all[3]
        CaM4Ca   = CaM₀_all[4]
    end

    @dynamics begin
        Ca'       = -kon_DMn * DMn_s * Ca + koff_DMn * CaDMn_s +
                    -kon_DMn * DMn_f * Ca + koff_DMn * CaDMn_f +
                    -kon_DMn * PP    * Ca + koff_PP  * CaPP    + 
                    -kon_D   * OGB5  * Ca + koff_D   * CaOGB5  +
                    -kon_Fluo4FF * Fluo4FF * Ca + koff_Fluo4FF * CaFluo4FF + 
                    -2kon_1   * CaM    * Ca^2 + 2koff_1 * CaM2Ca +
                     -kon_2   * CaM2Ca * Ca   +  koff_2 * CaM3Ca +
                     -kon_3   * CaM3Ca * Ca   +  koff_3 * CaM4Ca +
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
        CaM2Ca'   =  kon_1 * CaM * Ca^2 - koff_1 * CaM2Ca - kon_2 * CaM2Ca * Ca + koff_2 * CaM3Ca 
        CaM3Ca'   =                                         kon_2 * CaM2Ca * Ca - koff_2 * CaM3Ca - kon_3 * CaM3Ca * Ca + koff_3 * CaM4Ca
        CaM4Ca'   =                                                                                 kon_3 * CaM3Ca * Ca - koff_3 * CaM4Ca
    end

    @derived begin
        F_F0   ~ @. Normal((OGB5 + Fmax_Fmin * CaOGB5) / (OGB5₀ + Fmax_Fmin * CaOGB5₀), σ_faas)
        Ca_CaM ~ @. Normal((2CaM2Ca + 3CaM3Ca + 4CaM4Ca) / CaM_t, (2CaM2Ca + 3CaM3Ca + 4CaM4Ca) / CaM_t * σ_shifman)
    end
end;



Shifman_scheme = @model begin 
    @param begin
        tv_kon_1     ∈ Uniform(2,   9)
        tv_Kd_1      ∈ Normal(log10(7.9e-6),  1)
        tv_kon_2     ∈ Uniform(2,   9)    
        tv_Kd_2      ∈ Normal(log10(1.7e-6),  1)
        tv_kon_3     ∈ Uniform(2,   9)    
        tv_Kd_3      ∈ Normal(log10(35.e-6),  1)
        tv_kon_4     ∈ Uniform(2,   9)    
        tv_Kd_4      ∈ Normal(log10(8.9e-6),  1)

        μ             ∈ RealDomain(init=-1, lower=-5, upper=5)
        Ω             ∈ Constrained(Normal(0, 1), lower=1) 
        σ_faas        ∈ Constrained(Normal(0, 1), lower=1e-10, upper=1)
        σ_shifman     ∈ Constrained(Normal(0, 1), lower=1e-10, upper=1)
    end

    @random begin
        η     ~ Normal(μ, Ω)
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
        cov_τ_decay
        cov_Ca_i
        CaM_equil
        cov_Fluo4FF_t
        cov_kon_D
        cov_koff_D
        cov_Kd_Fluo4FF
        cov_kon_Fluo4FF
    end



    @pre begin
        U         = sigmoid(η)

        Fmax_Fmin = 39.364

        f_frac   = cov_f_frac 
        τ_f      = 1/cov_τ_f
        τ_s      = 1/cov_τ_s
        kon_DMn  = cov_kon_DMn
        koff_DMn = cov_koff_DMn
        Kd_DMn   = cov_Kd_DMn
        kon_PP   = cov_kon_DMn
        koff_PP  = cov_koff_PP

        DMn_t     = cov_DMn_t
        CaM_t     = cov_CaM_t
        OGB5_t    = cov_OGB5_t
        Fluo4FF_t = cov_Fluo4FF_t
        Ca_free   = cov_Ca_free

        τ_decay = cov_τ_decay
        Ca_i    = cov_Ca_i

        kon_D    = cov_kon_D         
        koff_D   = cov_koff_D          
        Kd_D     = koff_D / kon_D
    
        Kd_Fluo4FF = cov_Kd_Fluo4FF
        kon_Fluo4FF  = cov_kon_Fluo4FF
        koff_Fluo4FF = Kd_Fluo4FF * kon_Fluo4FF 

        Kd_1  = 10^tv_Kd_1
        kon_1 = 10^tv_kon_1
        koff_1= Kd_1 * kon_1

        Kd_2  = 10^tv_Kd_2 
        kon_2 = 10^tv_kon_2
        koff_2= Kd_2 * kon_2 

        Kd_3  = 10^tv_Kd_3 
        kon_3 = 10^tv_kon_3
        koff_3= Kd_3 * kon_3 

        Kd_4  = 10^tv_Kd_4 
        kon_4 = 10^tv_kon_4
        koff_4= Kd_4 * kon_4 

        DMn₀        = (Kd_DMn * DMn_t)  / (Ca_free + Kd_DMn)
        CaDMn₀      = DMn_t - DMn₀
        OGB5₀       = (Kd_D   * OGB5_t) / (Ca_free + Kd_D)
        CaOGB5₀     = OGB5_t - OGB5₀
        Fluo4FF₀    = (Kd_Fluo4FF * Fluo4FF_t) / (Ca_free + Kd_Fluo4FF)
        CaFluo4FF₀  = Fluo4FF_t - Fluo4FF₀       

        CaM₀_all = solve_Shifman_eq(kon_1, koff_1, kon_2, koff_2, kon_3, koff_3, kon_4, koff_4, Ca_free, CaM_t)

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
        CaM      = CaM₀_all[1]
        CaM1Ca   = CaM₀_all[2]
        CaM2Ca   = CaM₀_all[3]
        CaM3Ca   = CaM₀_all[4]
        CaM4Ca   = CaM₀_all[5]
    end

    @dynamics begin
        Ca'       = -kon_DMn * DMn_s * Ca + koff_DMn * CaDMn_s +
                    -kon_DMn * DMn_f * Ca + koff_DMn * CaDMn_f +
                    -kon_DMn * PP    * Ca + koff_PP  * CaPP    + 
                    -kon_D   * OGB5  * Ca + koff_D   * CaOGB5  +
                    -kon_Fluo4FF * Fluo4FF * Ca + koff_Fluo4FF * CaFluo4FF + 
                    -kon_1 * CaM    * Ca + koff_1 * CaM1Ca +
                    -kon_2 * CaM1Ca * Ca + koff_2 * CaM2Ca +
                    -kon_3 * CaM2Ca * Ca + koff_3 * CaM3Ca +
                    -kon_4 * CaM3Ca * Ca + koff_4 * CaM4Ca +
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
        CaM'      = -kon_1 * CaM * Ca + koff_1 * CaM1Ca 
        CaM1Ca'   =  kon_1 * CaM * Ca - koff_1 * CaM1Ca - kon_2 * CaM1Ca * Ca + koff_2 * CaM2Ca 
        CaM2Ca'   =                                       kon_2 * CaM1Ca * Ca - koff_2 * CaM2Ca - kon_3 * CaM2Ca * Ca + koff_3 * CaM3Ca 
        CaM3Ca'   =                                                                               kon_3 * CaM2Ca * Ca - koff_3 * CaM3Ca - kon_4 * CaM3Ca * Ca + koff_4 * CaM4Ca
        CaM4Ca'   =                                                                                                                       kon_4 * CaM3Ca * Ca - koff_4 * CaM4Ca
    end

    @derived begin
        F_F0   ~ @. Normal((OGB5 + Fmax_Fmin * CaOGB5) / (OGB5₀ + Fmax_Fmin * CaOGB5₀), σ_faas)
        Ca_CaM ~ @. Normal((CaM1Ca + 2CaM2Ca + 3CaM3Ca + 4CaM4Ca) / CaM_t, (CaM1Ca + 2CaM2Ca + 3CaM3Ca + 4CaM4Ca) / CaM_t * σ_shifman)
    end
end;



Pepke_m2_scheme = @model begin 
    @param begin
        tv_kon_TN     ∈ Normal(log10(7.7e5),    1)
        tv_Kd_TN      ∈ Normal(log10(1.93e-6),  1) 
        tv_kon_RN     ∈ Normal(log10(3.2e7),    1)
        tv_Kd_RN      ∈ Normal(log10(0.788e-6), 1)

        tv_kon_TC     ∈ Normal(log10(8.4e4),    1)
        tv_Kd_TC      ∈ Normal(log10(27.8e-6),  1)
        tv_kon_RC     ∈ Normal(log10(2.5e4),    1)
        tv_Kd_RC      ∈ Normal(log10(0.264e-6), 1)

        μ             ∈ RealDomain(init=-1, lower=-5, upper=5)
        Ω             ∈ Constrained(Normal(0, 1), lower=1)
        σ_faas        ∈ Constrained(Normal(0, 1), lower=1e-10, upper=1)
        σ_shifman     ∈ Constrained(Normal(0, 1), lower=1e-10, upper=1)
    end

    @random begin
        η     ~ Normal(μ, Ω)
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
        cov_τ_decay
        cov_Ca_i
        CaM_equil
        cov_Fluo4FF_t
        cov_kon_D
        cov_koff_D
        cov_Kd_Fluo4FF
        cov_kon_Fluo4FF
    end



    @pre begin
        U         = sigmoid(η)

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

        kon_D    = cov_kon_D
        koff_D   = cov_koff_D
        Kd_D     = koff_D / kon_D

        Kd_Fluo4FF   = cov_Kd_Fluo4FF
        kon_Fluo4FF  = cov_kon_Fluo4FF 
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

        kon_2C₀  = (kon_TC * kon_RC)   / (koff_TC + kon_RC * Ca_free)
        koff_2C₀ = (koff_TC * koff_RC) / (koff_TC + kon_RC * Ca_free)

        kon_2N₀  = (kon_TN * kon_RN)   / (koff_TN + kon_RN * Ca_free)
        koff_2N₀ = (koff_TN * koff_RN) / (koff_TN + kon_RN * Ca_free)


        c0       = solve_Blackwell_CN_eqs(kon_2C₀, koff_2C₀, kon_2N₀, koff_2N₀, Ca_free, CaM_t, CaM_equil)

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
        CaM0     = c0[1]
        CaM2C    = c0[2]
        CaM2N    = c0[3]
        CaM4     = c0[4]
    end

    @vars begin
        kon_2C  = (kon_TC * kon_RC)   / (koff_TC + kon_RC * Ca)
        koff_2C = (koff_TC * koff_RC) / (koff_TC + kon_RC * Ca)

        kon_2N  = (kon_TN * kon_RN)   / (koff_TN + kon_RN * Ca)
        koff_2N = (koff_TN * koff_RN) / (koff_TN + kon_RN * Ca)
    end

    @dynamics begin
        Ca'       = -kon_DMn * DMn_s * Ca + koff_DMn * CaDMn_s +
                    -kon_DMn * DMn_f * Ca + koff_DMn * CaDMn_f +
                    -kon_DMn * PP    * Ca + koff_PP  * CaPP    + 
                    -kon_D   * OGB5  * Ca + koff_D   * CaOGB5  +
                    -kon_Fluo4FF * Fluo4FF * Ca + koff_Fluo4FF * CaFluo4FF +
                    -2kon_2C * CaM0  * Ca^2 + 2koff_2C * CaM2C +
                    -2kon_2N * CaM0  * Ca^2 + 2koff_2N * CaM2N +
                    -2kon_2N * CaM2C * Ca^2 + 2koff_2N * CaM4 +
                    -2kon_2C * CaM2N * Ca^2 + 2koff_2C * CaM4 +
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
        F_F0 ~ @. Normal((OGB5 + Fmax_Fmin * CaOGB5) / (OGB5₀ + Fmax_Fmin * CaOGB5₀), σ_faas)
        Ca_CaM ~ @. Normal((2CaM2C + 2CaM2N + 4CaM4) / CaM_t, (2CaM2C + 2CaM2N + 4CaM4) / CaM_t * σ_shifman)
    end
end;



Pepke_m1_scheme = @model begin 
    @param begin
        tv_kon_TN     ∈ Normal(log10(7.7e5),    1)
        tv_Kd_TN      ∈ Normal(log10(1.93e-6),  1) 
        tv_kon_RN     ∈ Normal(log10(3.2e7),    1)
        tv_Kd_RN      ∈ Normal(log10(0.788e-6), 1)

        tv_kon_TC     ∈ Normal(log10(8.4e4),    1)
        tv_Kd_TC      ∈ Normal(log10(27.8e-6),  1)
        tv_kon_RC     ∈ Normal(log10(2.5e4),    1)
        tv_Kd_RC      ∈ Normal(log10(0.264e-6), 1)

        μ             ∈ RealDomain(init=-1, lower=-5, upper=5)
        Ω             ∈ Constrained(Normal(0, 1), lower=1)
        σ_faas        ∈ Constrained(Normal(0, 1), lower=1e-10, upper=1)
        σ_shifman     ∈ Constrained(Normal(0, 1), lower=1e-10, upper=1)
    end

    @random begin
        η     ~ Normal(μ, Ω)
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
        cov_τ_decay
        cov_Ca_i
        CaM_equil
        cov_Fluo4FF_t
        cov_kon_D
        cov_koff_D
        cov_Kd_Fluo4FF
        cov_kon_Fluo4FF
    end



    @pre begin
        U         = sigmoid(η)

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

        kon_D    = cov_kon_D    
        koff_D   = cov_koff_D
        Kd_D     = koff_D / kon_D

        Kd_Fluo4FF = cov_Kd_Fluo4FF
        kon_Fluo4FF  = cov_kon_Fluo4FF
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

        N₀_all = get_Pepke_m1_lobe_eqs(Ca_free, Kd_TN, Kd_RN, CaM_t, CaM_equil)
        C₀_all = get_Pepke_m1_lobe_eqs(Ca_free, Kd_TC, Kd_RC, CaM_t, CaM_equil)

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
        N0       = N₀_all[1]
        N1       = N₀_all[2]
        N2       = N₀_all[3]
        C0       = C₀_all[1]
        C1       = C₀_all[2]
        C2       = C₀_all[3]
    end

    @dynamics begin
        Ca'       = -kon_DMn * DMn_s * Ca + koff_DMn * CaDMn_s +
                    -kon_DMn * DMn_f * Ca + koff_DMn * CaDMn_f +
                    -kon_DMn * PP    * Ca + koff_PP  * CaPP    + 
                    -kon_D   * OGB5  * Ca + koff_D   * CaOGB5  +
                    -kon_Fluo4FF * Fluo4FF * Ca + koff_Fluo4FF * CaFluo4FF + 
                    -kon_TN * N0 * Ca + koff_TN * N1 +
                    -kon_TC * C0 * Ca + koff_TC * C1 +     
                    -kon_RN * N1 * Ca + koff_RN * N2 +
                    -kon_RC * C1 * Ca + koff_RC * C2 +
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
        N0'       = -kon_TN * N0 * Ca + koff_TN * N1
        N1'       =  kon_TN * N0 * Ca - koff_TN * N1 - kon_RN * N1 * Ca + koff_RN * N2
        N2'       =                                  + kon_RN * N1 * Ca - koff_RN * N2
        C0'       = -kon_TC * C0 * Ca + koff_TC * C1
        C1'       =  kon_TC * C0 * Ca - koff_TC * C1 - kon_RC * C1 * Ca + koff_RC * C2
        C2'       =                                  + kon_RC * C1 * Ca - koff_RC * C2
    end

    @derived begin
        signal = @. (OGB5 + Fmax_Fmin * CaOGB5) / (OGB5₀ + Fmax_Fmin * CaOGB5₀)
        F_F0 ~ @. Normal(signal, σ_faas)
        Ca_CaM ~ @. Normal((N1 + 2N2 + C1 + 2C2) / CaM_t, (N1 + 2N2 + C1 + 2C2) / CaM_t * σ_shifman)
    end
end;


Faas_scheme = @model begin 
    @param begin
        tv_kon_TN     ∈ Normal(log10(7.7e5),    1)
        tv_Kd_TN      ∈ Normal(log10(1.93e-6),  1) 
        tv_kon_RN     ∈ Normal(log10(3.2e7),    1)
        tv_Kd_RN      ∈ Normal(log10(0.788e-6), 1)

        tv_kon_TC     ∈ Normal(log10(8.4e4),    1)
        tv_Kd_TC      ∈ Normal(log10(27.8e-6),  1)
        tv_kon_RC     ∈ Normal(log10(2.5e4),    1)
        tv_Kd_RC      ∈ Normal(log10(0.264e-6), 1)

        μ             ∈ RealDomain(init=-1, lower=-5, upper=5)
        Ω             ∈ Constrained(Normal(0, 1), lower=1)
        σ_faas        ∈ Constrained(Normal(0, 1), lower=1e-10, upper=1)
        σ_shifman     ∈ Constrained(Normal(0, 1), lower=1e-10, upper=1)
    end

    @random begin
        η     ~ Normal(μ, Ω)
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
        cov_τ_decay
        cov_Ca_i
        CaM_equil
        cov_Fluo4FF_t
        cov_kon_D
        cov_koff_D
        cov_Kd_Fluo4FF
        cov_kon_Fluo4FF
    end



    @pre begin
        U         = sigmoid(η)

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

        kon_D    = cov_kon_D
        koff_D   = cov_koff_D
        Kd_D     = koff_D / kon_D

        Kd_Fluo4FF   = cov_Kd_Fluo4FF
        kon_Fluo4FF  = cov_kon_Fluo4FF
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
        F_F0 ~ @. Normal(signal, σ_faas)
        Ca_CaM ~ @. Normal((NtNr + 2NrNr + CtCr + 2CrCr) / CaM_t, (NtNr + 2NrNr + CtCr + 2CrCr) / CaM_t * σ_shifman)
    end
end;



Byrne_scheme = @model begin ### Byrne scheme
    @param begin
        tv_k01_N      ∈ Normal(log10(275e3),   1)
        tv_k02_N      ∈ Normal(log10(275e3),   1)
        tv_k13_N      ∈ Normal(log10(507.2e3), 1)
        tv_k23_N      ∈ Normal(log10(500e3),   1)

        tv_K01d_N     ∈ Normal(log10(33e-6),   1)
        tv_K13d_N     ∈ Normal(log10(230e-6),  1)
        tv_K02d_N     ∈ Normal(log10(3.45e-6), 1)

        tv_k01_C      ∈ Normal(log10(275e3),   1)
        tv_k02_C      ∈ Normal(log10(275e3),   1)
        tv_k13_C      ∈ Normal(log10(3.71e3),  1)
        tv_k23_C      ∈ Normal(log10(118e3),   1)

        tv_K01d_C     ∈ Normal(log10(18.5e-6), 1)
        tv_K02d_C     ∈ Normal(log10(116e-6),  1)
        tv_K13d_C     ∈ Normal(log10(0.38e-6), 1)


        μ             ∈ RealDomain(init=-1, lower=-5, upper=5)
        Ω             ∈ Constrained(Normal(0, 1), lower=1)
        σ_faas        ∈ Constrained(Normal(0, 1), lower=1e-10, upper=1)
        σ_shifman     ∈ Constrained(Normal(0, 1), lower=1e-10, upper=1)
    end

    @random begin
        η     ~ Normal(μ, Ω)
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
        cov_τ_decay
        cov_Ca_i
        CaM_equil
        cov_Fluo4FF_t
        cov_kon_D
        cov_koff_D
        cov_Kd_Fluo4FF
        cov_kon_Fluo4FF
    end



    @pre begin
        U         = sigmoid(η)

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

        kon_D    = cov_kon_D
        koff_D   = cov_koff_D
        Kd_D     = koff_D / kon_D

        Kd_Fluo4FF = cov_Kd_Fluo4FF
        kon_Fluo4FF  = cov_kon_Fluo4FF
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
        F_F0 ~ @. Normal((OGB5 + Fmax_Fmin * CaOGB5) / (OGB5₀ + Fmax_Fmin * CaOGB5₀), σ_faas)
        Ca_CaM ~ @. Normal((N1 + N2 + 2N3 + C1 + C2 + 2C3) / CaM_t, (N1 + N2 + 2N3 + C1 + C2 + 2C3) / CaM_t * σ_shifman)
    end
end;


Blackwell_const_names = (
    :tv_kon_1,
    :tv_Kd_1,
    :tv_kon_2,
    :tv_Kd_2,
);
Bhalla_const_names = (
    :tv_kon_1,
    :tv_Kd_1,
    :tv_kon_2,
    :tv_Kd_2,
    :tv_kon_3,
    :tv_Kd_3,
);
Shifman_const_names = (
    :tv_Kd_1,
    :tv_Kd_2,
    :tv_Kd_3,
    :tv_Kd_4,
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
)

Blackwell_params = (;
    tv_kon_1  = log10(4e3),
    tv_Kd_1   = log10(2e-6),
    tv_kon_2  = log10(1e5),
    tv_Kd_2   = log10(1.1e-5),
    μ         = 1,
    Ω         = 2,
    σ_faas    = 0.5,
    σ_shifman = 0.5
);
Bhalla_params = (;
    tv_kon_1  = log10(72e3),
    tv_Kd_1   = log10(1e-6),
    tv_kon_2  = log10(3.6e3),
    tv_Kd_2   = log10(2.8e-6),
    tv_kon_3  = log10(4.65e2),
    tv_Kd_3   = log10(2.1e-5),
    μ         = 1,
    Ω         = 2,
    σ_faas    = 0.5,
    σ_shifman = 0.5
);
Shifman_params = (;
    tv_kon_1  = log10(1e5),
    tv_Kd_1   = log10(7.9e-6),
    tv_kon_2  = log10(1e5),
    tv_Kd_2   = log10(1.7e-6),
    tv_kon_3  = log10(1e5),
    tv_Kd_3   = log10(35e-6),
    tv_kon_4  = log10(1e5),
    tv_Kd_4   = log10(8.9e-6),
    μ         = 1,
    Ω         = 2,
    σ_faas    = 0.5,
    σ_shifman = 0.5
);
Faas_params = (;
    tv_kon_TN = 5.9,  
    tv_Kd_TN  = -3.7,
    tv_kon_RN = 7.5, 
    tv_Kd_RN  = -6.1,
    tv_kon_TC = 4.9, 
    tv_Kd_TC  = -4.6,
    tv_kon_RC = 4.4, 
    tv_Kd_RC  = -6.6,
    μ         = 1,
    Ω         = 2,
    σ_faas    = 0.5,
    σ_shifman = 0.5
);
Pepke_params = (;     ###https://www.ebi.ac.uk/biomodels/MODEL1001150000
    tv_kon_TN = log10(100e3),  
    tv_Kd_TN  = log10(25e-6),
    tv_kon_RN = log10(150e3), 
    tv_Kd_RN  = log10(5e-6),
    tv_kon_TC = log10(4e3), 
    tv_Kd_TC  = log10(10e-6),
    tv_kon_RC = log10(10e3), 
    tv_Kd_RC  = log10(0.925e-6),
    μ         = 1,
    Ω         = 2,
    σ_faas    = 0.5,
    σ_shifman = 0.5
);
Byrne_params = (;
    tv_k01_N    = log10(275e3),
    tv_K01d_N   = log10(33e-6),
    tv_k02_N    = log10(275e3),
    tv_K02d_N   = log10(229.88e-6),
    tv_k13_N    = log10(507.2e3),
    tv_K13d_N   = log10(3.45e-6),
    tv_k23_N    = log10(500e3),
###    tv_K23d_N   = log10(0.5e-6),
    tv_k01_C    = log10(2.75e5),
    tv_K01d_C   = log10(18.5e-6),
    tv_k02_C    = log10(2.75e5),
    tv_K02d_C   = log10(116e-6),
    tv_k13_C    = log10(3.71e3),
    tv_K13d_C   = log10(0.38e-6),
    tv_k23_C    = log10(1.18e5),
###    tv_K23d_C   = log10(0.06e-6),
    μ         = 1,
    Ω         = 2,
    σ_faas    = 0.5,
    σ_shifman = 0.5
)