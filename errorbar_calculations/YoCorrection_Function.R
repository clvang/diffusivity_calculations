# This script calculates changes in the mass fraction between needle retraction
# and the beginning of igniters being energized.  Temperature is assumed constant between
# this time window (e.g. T=298 K).

YoCorrection_Function <- function(N, cf, sol_id, 
                                  dod1_ratio, yo, UYo_sq,
                                  Udo_sq, do_measured ){

  source('/Users/changvang/mygitFiles/diffusivity_calculations/errorbar_calculations/staticProps.R')

  Yo1_Solve <- function(rho_l_1, rho_h_1, Y_ratio_o, d_ratio, err_tol, N){
    Yo1_values <- .Fortran("Yo1_bisectSolve",
          rho_l_1=as.double(rho_l_1),
          rho_h_1=as.double(rho_h_1),
          Y_ratio_o = as.double(Y_ratio_o),
          d_ratio = as.double(d_ratio),
          err_tol = as.double(err_tol),
          N = as.integer(N),
          Yl1_vector = double(N) )
    #there are the returned values
    list(Yo_mcVals=Yo1_values$Yl1_vector)         #values for D
  }

  if (!is.loaded('Yo1_bisectSolve')){
    dyn.load("Yo1_bisectSolve.so")
  }


  # read static thermodynamic properties of species. Can
  # be either heptane/hexadecane or propanol/glycerol based on "sol_id"
    static_properties <- staticProps(sol_id)

  # nominal ratio of diameters
    d_ratio <- (dod1_ratio)^3     

  # nominal initial mass fractions at needle retraction
    Y_lo <- yo          # low volatility component
    Y_ho <- 1-Y_lo                              # high volatility component

    sigma_Y <- sqrt(UYo_sq) / cf
    Ylo_mcN95 <- rnorm(n=N, mean=Y_lo, sd=sigma_Y)
    Yho_mcN95 <- 1 - Ylo_mcN95

  # nominal INITIAL DENSITIES of hight and low volatility components
    # Low volatility component density at 298 K, 1 atm (Hexadecane or Glycerol)
    rho_l_0 <- static_properties$rho_Ao   
    # High volatility component density at 298 K, 1 atm (Heptane or Propanol)
    rho_h_0 <- static_properties$rho_Bo 

    # these uncertainties for static properties of densities 
    # are consistent with those defined in tdelay_tdecay.R, 
    # where static properties are used to calculate viscous decay times
    rho_relative_error <- 0.01                #density 10% relative error for 95% confidnece
    u_rho_l0 <- rho_l_0*rho_relative_error    #relative error U_rho/rho = 1% for 95% confidence
    u_rho_h0 <- rho_h_0*rho_relative_error    #relative error U_rho/rho = 1% for 95% confidence
    sigma_rho_l0 <- u_rho_l0 / cf 
    sigma_rho_h0<- u_rho_h0 / cf 

    rholo_mcN95 <- rnorm(n=N, mean=rho_l_0, sd=sigma_rho_l0)
    rhoho_mcN95 <- rnorm(n=N, mean=rho_h_0, sd= sigma_rho_h0)


  # nominal DENSITY OF HIGH VOLATILITY COMP. AT TEMPERATURE WHERE
  # IGNITER FIRST ENERGIEZES (E.G. 298 K)
    rho_h_1 <- rho_h_0  #kg/m^3

  # nominal DENSITY OF LOW VOLATILITY COMP. AT TEMPERATURE WHERE
  # IGNITER FIRST ENERGIZES (E.G. 298 K)
    rho_l_1 <- rho_l_0  #kg/m^3

    # these uncertainties are also consisten with those
    # used in tdelay_decay.R 
    u_rho_l1 <- rho_l_1*rho_relative_error    #relative error U_rho/rho = 1% for 95% confidence
    u_rho_h0 <- rho_h_1*rho_relative_error    #relative error U_rho/rho = 1% for 95% confidence
    sigma_rho_l1 <- u_rho_l1 / cf 
    sigma_rho_h1<- u_rho_h0 / cf 

    rhol1_mcN95 <- rnorm(n=N, mean=rho_l_1, sd=sigma_rho_l1)
    rhoh1_mcN95 <- rnorm(n=N, mean=rho_h_1, sd= sigma_rho_h1)

  # ---------- generate random variables for each of the variables listed above-------

  # GENERATING RANDOM NUMBERS FOR DIAMETER RATIO
  sigma_do <- sqrt(Udo_sq) / cf
  sigma_d1 <- sigma_do 

  do_mcN95 <- rnorm(n=N, mean=do_measured, sd=sigma_do)
  d1_mcN95 <- rnorm(n=N, mean=(do_measured/d_ratio), sd=sigma_d1)

  # factor <- 2
  # do_temp <- rnorm( n=N*factor, mean=do_measured, sd= sigma_do )
  # d1_mcN95 <- 1
  # i <- 0
  # while(length(d1_mcN95) < N){
  #   i <- i + 1
  #   #set.seed(i)
  #   d1_mcN95 <- rnorm(n=N*factor, mean = (do_measured/d_ratio), sd=sigma_d1 )
  #   index <- which(do_temp >= d1_mcN95)
  #   d1_mcN95 <- d1_mcN95[index]
  #   # print(i)
  # }

  # d1_mcN95 <- d1_mcN95[1:N]   #select only N values from the within range values
  # do_mcN95 <- do_temp[index]
  # do_mcN95 <- do_mcN95[1:N]


  d_ratio_mcN95 <- (do_mcN95/d1_mcN95)^3 #rnorm(n=N, mean = d_ratio, sd= d_ratio*0.01 ) #

  # define terms in equation
  # Y_ratio_o <- (Y_ho - 1)/( (Y_ho/rho_h_0) + (Y_lo/rho_l_0) )
  Y_ratio_o_mcN95 <- (Yho_mcN95 - 1)/( (Yho_mcN95/rhoho_mcN95) + 
                    (Ylo_mcN95/rholo_mcN95) )

  err_tol <- 1e-6
  Yo_mcN95 <- Yo1_Solve(rho_l_1=rhol1_mcN95, 
                            rho_h_1=rhoh1_mcN95, 
                            Y_ratio_o=Y_ratio_o_mcN95, 
                            d_ratio=d_ratio_mcN95, 
                            err_tol=err_tol, 
                            N=N)$Yo_mcVals

  # err_tol <- 1e-6
  # Y_h1_nominal <- Yo1_Solve(rho_l_1=rho_l_1, 
  #                           rho_h_1=rho_h_1, 
  #                           Y_ratio_o=Y_ratio_o, 
  #                           d_ratio=d_ratio, 
  #                           err_tol=err_tol, 
  #                           N=1)$Yo_mcVals
  # Yo_mcN95 <- 1 - Y_h1_mcN95

  list=c(Yo_mcN95=Yo_mcN95) 

}
  