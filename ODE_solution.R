library(deSolve)  

##############
# User input #
##############

weight <- 70 #in kg
gender <- 0 #0 for male | 1 for female
dose <- 10  # in mg
infusion_time <- 5/60 #infusion time in hours
C0_MU <- 0
C0_AD <- 0
C0_GO <- 0
C0_SK <- 0
C0_HT <- 0
C0_BR <- 0
C0_KI <- 0
C0_RE <- 0
C0_ST <- 0
C0_IN <- 0
C0_LU <- 0
C0_LI <- 0
C0_ART <- 0
C0_VEN <- 0
user_input <-data.frame(weight ,gender, dose, infusion_time, C0_MU, C0_AD, C0_GO, C0_SK, C0_HT, C0_BR, C0_KI,
                      C0_RE, C0_ST, C0_IN, C0_LU, C0_LI, C0_ART, C0_VEN )

#####################
# Compartment names #
#####################

comp_names <- c("MU", "AD", "GO" , "SK", "HT", "BR", "KI", "RE", "ST", "IN","LU","LI","ART","VEN")
      

######################################################
# Calculation or organ flows and volumes #
######################################################

covariates <- function(w,gender){
        k1<-c(9.61e-02,-4.88e-06,3.05e-10,-3.62e-15,1.22e-20,0,0.17)
        k2<-c(3.95e-02,1.59e-05,-6.99e-10,1.09e-14,-5.26e-20,0,0.05)
        k3<-c(1.67e-04,6.2e-10,-6.54e-13,2.48e-17,-2.85e-22,1.03e-27,0.001)
        k4<-c(1.07e-01,-3.26e-06,6.11e-11,-5.43e-16,1.83e-21,0,0.05)
        k5<-c(8.53e-03,-4.07e-07,1.4e-11,-1.9e-16,1.05e-21,-1.94e-27,0.04)
        k6<-c(1.19e-01,-3.51e-06,4.28e-11,-1.82e-16,0,0,0.12)
        k7<-c(7.31e-03,-8.29e-08,5.71e-13,0,0,0,0.19)
        k8<-rep(0,6)
        k9<-c(1.88e-03,8.76e-08,-2.52e-12,1.86e-17,0,0,0.01)
        k10<-c(1.74e-02,-5.3e-07,1.18e-11,-6.74e-17,0,0,0.14)
        k11<-c(1.67e-02,-9.96e-08,-1.09e-13,1.13e-17,0,0,1)
        k12<-c(3.49e-02,-3.23e-07,2.13e-12,0,0,0,0.065)
        k13<-c(3.66e-02,-3.44e-07,5.00e-12, -2.59e-17,0.0,0.0)
        k14<-c(5.49e-02,-5.15e-07, 7.50e-12,-3.87e-17,0.0,0.0)
        
        k1_W<-c(1.17e-01,-3.59e-06,3.19e-10,-3.55e-15,-7.58e-22,0.0)
        k2_W<-c(5.91e-02,1.20e-05,-5.80e-10,1.12e-14,-6.36e-20,0.0)
        k3_W<-c(1.94e-04,-8.32e-09,3.15e-13,0.0,0.0,0.0)
        k4_W<-c(9.54e-02,-1.7e-06,-1.64e-13,2.64e-16,-1.49e-21,0.0)
        k5_W<-c(5.72e-03,-1.02e-07,2.53e-12,-2.71e-17,9.29e-23,0.0)
        k6_W<-c(1.12e-01,-3.33e-06,4.04e-11,-1.70e-16,0.0,0.0)
        k7_W<-c(8.04e-03,-1.38e-07,2.19e-12,-1.34e-17,0.0,0.0)
        k8_W<-rep(0,6)
        k9_W<-c(1.88e-03,8.76e-08,-2.52e-12,1.86e-17,0.0,0.0)
        k10_W<-c(1.89e-02,-6.62e-07,1.56e-11,-9.87e-17,0.0,0.0)
        k11_W<-c(1.74e-02,-7.14e-08,-6.78e-14,0.0,0.0,0.0)
        k12_W<-c(3.59e-02,-4.76e-07,8.50e-12,-5.45e-17,0.0,0.0)
        k13_W<-c(3.66e-02,-3.44e-07,5.00e-12, -2.59e-17,0.0,0.0)
        k14_W<-c(5.49e-02,-5.15e-07, 7.50e-12,-3.87e-17,0.0,0.0)
        
        #density<-c(1.041,0.916,1,1,1.03,1.035,1.05,1,1.05,1.042,1.05,1,1,1)#in kg/L
        density <- rep(1,14)
        flow_frac<-c(0.17,0.05,0.001,0.05,0.04,0.12,0.19,0,0.01,0.14,1,0.065)
        
        const<-list(k1,k2,k3,k4,k5,k6,k7,k8,k9,k10,k11,k12,k13,k14)
        const_W<-list(k1_W,k2_W,k3_W,k4_W,k5_W,k6_W,k7_W,k8_W,k9_W,k10_W,k11_W,k12_W,k13_W,k14_W)
        vol<-rep(0,14)
        flow<-rep(0,12)
        TF<-(187*(w)^0.81)*60/1000 #Total flow
        weight <- w
        
        if (gender==0){
                if (w>75){
                        w<-75*1000
                }
                else{
                        w<-w*1000
                }
                co<-const
                for (i in 1:14){
                        
                        vol[i]<-(co[[i]][1]+co[[i]][2]*w+co[[i]][3]*w^2+co[[i]][4]*w^3+co[[i]][5]*w^4+co[[i]][6]*w^5)
                        vol[i]<-vol[i]*weight/density[i]
                }
        }
        
        else if (gender==1){
                if (w>65){
                        w<-65*1000
                }
                else{
                        w<-w*1000
                }
                co<-const_W
                for (i in 1:14){
                        
                        vol[i]<-(co[[i]][1]+co[[i]][2]*w+co[[i]][3]*w^2+co[[i]][4]*w^3+co[[i]][5]*w^4+co[[i]][6]*w^5)
                        vol[i]<-vol[i]*weight/density[i]
                }
        }
        
        
        
        for (i in 1:12){
                
                flow[i]<-TF*flow_frac[i]
        }
        
        flow[8]<-TF-sum(flow)+TF
        vol[8]<-weight-sum(vol)
        combine<-c(flow,vol)
        
        return(combine)
        
}



#################
# ODEs system   #
#################

odes <- function(time,C,params){
        
        dCdt<-rep(0,14)
        
        Q_MU<-params[1]
        Q_AD<-params[2]
        Q_TE<-params[3]
        Q_SK<-params[4]
        Q_HT<-params[5]
        Q_BR<-params[6]
        Q_KI<-params[7]
        Q_RE<-params[8]
        Q_ST<-params[9]
        Q_SPL<-params[10]
        Q_LU<-params[11]
        Q_LI<-params[12]
        Q_ART<-params[11]
        Q_VEN<-params[11]    
        
        V_MU<-params[13]
        V_AD<-params[14]
        V_TE<-params[15]
        V_SK<-params[16]
        V_HT<-params[17]
        V_BR<-params[18]
        V_KI<-params[19]
        V_RE<-params[20]
        V_ST<-params[21]
        V_SPL<-params[22]
        V_LU<-params[23]
        V_LI<-params[24]
        V_ART<-params[25]
        V_VEN<-params[26]
        T_inf<-params[27]
        dose <- params[28]
        dose <- dose *1e06 #dose in ng
        
        Kp_init <- c(0.56, 3.34, 0.76, 0.46, 0.86, 0.32, 0.73, 0.5, 0.75, 0.53, 0.71, 1.35) 
        a <- 0.12 # scaling factor
        Kp <- Kp_init * a # updated Kp vector
        Kp[2] <- 3.32 # Kp of adipose compartment
        
        
        
        KP_MU<-Kp[1]
        KP_AD<-Kp[2]
        KP_TE<-Kp[3]
        KP_SK<-Kp[4]
        KP_HT<-Kp[5]
        KP_BR<-Kp[6]
        KP_KI<-Kp[7]
        KP_RE<-Kp[8]
        KP_ST<-Kp[9]
        KP_SPL<-Kp[10]
        KP_LU<-Kp[11]
        KP_LI<-Kp[12]
        CL<-64.2
        Q_HA<-params[12]
        Q_HP<-params[9]+params[10]
        Q_H<-Q_HA+Q_HP
        f_UB<-0.015/0.65 
       
        
        dCdt[1]<-(-Q_MU*C[1]/(KP_MU*V_MU))+(Q_MU*C[13]/V_MU)
        dCdt[2]<-(-Q_AD*C[2]/(KP_AD*V_AD))+(Q_AD*C[13]/V_AD)
        dCdt[3]<-(-Q_TE*C[3]/(KP_TE*V_TE))+(Q_TE*C[13]/V_TE)
        dCdt[4]<-(-Q_SK*C[4]/(KP_SK*V_SK))+(Q_SK*C[13]/V_SK)
        dCdt[5]<-(-Q_HT*C[5]/(KP_HT*V_HT))+(Q_HT*C[13]/V_HT)
        dCdt[6]<-(-Q_BR*C[6]/(KP_BR*V_BR))+(Q_BR*C[13]/V_BR)
        dCdt[7]<-(-Q_KI*C[7]/(KP_KI*V_KI))+(Q_KI*C[13]/V_KI)
        dCdt[8]<-(-Q_RE*C[8]/(KP_RE*V_RE))+(Q_RE*C[13]/V_RE)
        dCdt[9]<-(-Q_ST*C[9]/(KP_ST*V_ST))+(Q_ST*C[13]/V_ST)
        dCdt[10]<-(-Q_SPL*C[10]/(KP_SPL*V_SPL))+(Q_SPL*C[13]/V_SPL)
        dCdt[11]<-(-Q_LU*C[11]/(KP_LU*V_LU))+(Q_VEN*C[14]/V_LU)
        dCdt[12]<-(Q_HA*C[13]/V_LI)+(Q_ST*C[9]/(KP_ST*V_LI))+(Q_SPL*C[10]/(KP_SPL*V_LI))-(((Q_H+CL*f_UB)/(KP_LI*V_LI))*C[12])
        dCdt[13]<-(-Q_ART*C[13]/V_ART)+(Q_LU*C[11]/(V_ART*KP_LU))
        
        if (time<=T_inf){
                dCdt[14]=(-Q_LU*C[14]/V_VEN)+(Q_H*C[12]/(KP_LI*V_VEN))+(Q_KI*C[7]/(KP_KI*V_VEN))+
                        (Q_MU*C[1]/(KP_MU*V_VEN))+(Q_AD*C[2]/(KP_AD*V_VEN))+(Q_SK*C[4]/(KP_SK*V_VEN))+
                        (Q_TE*C[3]/(KP_TE*V_VEN))+(Q_HT*C[5]/(KP_HT*V_VEN))+(Q_BR*C[6]/(KP_BR*V_VEN))+
                        (Q_RE*C[8]/(KP_RE*V_VEN))+(dose/(V_VEN*1000*T_inf))  #units are in ng/mL
        }
        else{
                
                dCdt[14]<-(-Q_LU*C[14]/V_VEN)+(Q_H*C[12]/(KP_LI*V_VEN))+(Q_KI*C[7]/(KP_KI*V_VEN))+
                        (Q_MU*C[1]/(KP_MU*V_VEN))+(Q_AD*C[2]/(KP_AD*V_VEN))+(Q_SK*C[4]/(KP_SK*V_VEN))+
                        (Q_TE*C[3]/(KP_TE*V_VEN))+(Q_HT*C[5]/(KP_HT*V_VEN))+(Q_BR*C[6]/(KP_BR*V_VEN))+
                        (Q_RE*C[8]/(KP_RE*V_VEN))
        }
        
        list(dCdt)        
        #return(dCdt)
}

##############################################
initial_concentration<- c(user_input$C0_MU, user_input$C0_AD,user_input$C0_GO, user_input$C0_SK, user_input$C0_HT,
                          user_input$C0_BR, user_input$C0_KI, user_input$C0_RE, user_input$C0_ST,
                          user_input$C0_IN, user_input$C0_LU,
                          user_input$C0_LI, user_input$C0_ART, user_input$C0_VEN)

params<-c(covariates(user_input$weight[1],user_input$gender[1]),user_input$infusion_time[1],
          user_input$dose) 


sample_time <- c(0, 5/60, 0.25, 0.5, 0.75, 1, 1.5, 2, 3, 4, 6, 8, 10, 12, 24, 36, 48, 72) # in hours
solution <- ode(y = initial_concentration, times = sample_time, func = odes, parms = params)
solution

