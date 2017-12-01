
functionPlot <- function(tau_vector, LHS_vector, eps_vector, expname)

{
	# functionEval <- function(tau_vector_input, eps_vector_input, LHS_vector_input, N, num_pts, expname){
	functionEval <- function(tau_input, eps_input, LHS_input){
	# SUBROUTINE functional_eval( tau_vector, LHS_vector, eps_vector, N, num_pts, expname )
		# FxValue_Generate <- .Fortran("functional_eval",
		# 					tau_vector = as.double(tau_vector_input),
		# 					LHS_vector = as.double(LHS_vector_input),
		# 					eps_vector = as.double(eps_vector_input),
		# 					N = as.integer(N),
		# 					num_pts = as.integer(num_pts),
		# 					expname = as.character(expname) )
		Fx_value <- .Fortran("Fx_Eval",
					tau=as.double(tau_input),
					eps=as.double(eps_input),
					LHS = as.double(LHS_input),
					Fx = double(1) )
		list(Fx_value_returned=Fx_value$Fx)         #values for D
	}

	if (!is.loaded('Fx_Eval')){
		dyn.load("Fx_Eval.so")
	}
	# if (!is.loaded('functional_eval')){
	# 	dyn.load("functional_eval.so")
	# }	
	print("======== evaluating function to inspect solution space ==========")
 
 	eps_min <- 0.0028 #min(eps_vector)*0.85
 	eps_max <- 0.80   #max(eps_vector)*1.10
 	num_pts <- 20
 	deps <- (eps_max - eps_min)/(num_pts - 1)
 	eps_values <- numeric(num_pts)
 	for( i in 1:num_pts ){
 		eps_values[i] <- eps_min + (i-1)*deps
 	}

 	n <- length(tau_vector)
 	# functionEval(tau_vector_input = tau_vector,
 	# 				eps_vector_input = eps_vector,
 	# 				LHS_vector_input = LHS_vector,
 	# 				N = n, 
 	# 				num_pts = num_pts,
 	# 				expname = expname)

 	Fx <- numeric(num_pts)
 	accumulator <- 0
 	pb <- txtProgressBar(min = 1, max = n*num_pts)
	for(i in 1:n ){
		for(j in 1:num_pts){
			accumulator <- accumulator + 1
			setTxtProgressBar(pb,accumulator)
			Fx[j] <- functionEval(tau_input=tau_vector[i],
					eps_input=eps_values[j],
					LHS_input=LHS_vector[i])$Fx_value_returned
		}
		if( i == 1){
			dev.new()
			plot(eps_values, Fx)
			title(main=expname)
		}else{
			lines(eps_values, Fx)
		}
		# readline("--press enter--")
	}	
	dev.copy(png,paste0(expname,"_functionValues.png"))
	dev.off()

}


