"newtraph" <- 
function(ders, dat, x0, ...)
{

# 	cat("While N-R may be used for either minimization or\nmaximization")
# 	cat("The checks for progress in this function are written")
# 	cat("for maximization.  If you want to minimize, chage your")
# 	cat("derivative calculations (multiply by -1).")
        crit1<-1e-10
        crit2<-1e-06
        crit3<-1e-06
        c1<-0
        c2<-0
        c3<-0
	curnt <- x0
        nump<-length(x0)
	k <- 0
	repeat {
		k <- k + 1
#                 cat(" ",fill=T)
#                 cat(" ",fill=T)
#                 cat("Current estimates beginning iteration ", k, ":", fill = T)
# 		cat(curnt, fill = T)
#                 cat(" ",fill=T)
		int <- ders(curnt, dat, ...)
                logL<-int[[1]]
                gi<-int[[2]]
# 		cat("Log likelihood for these estimates: ", fill = T)
# 		cat(logL, fill = T)
#                 cat("Gradient for these estimates: ",fill=T)
#                 cat(gi,fill=T)
# 		cat(" ", fill = T)
		Gi <- int[[3]]
		GiI <- solve(Gi)
		step <- GiI %*% gi
		new <- curnt - step
		sc <- 1
		repeat {
                        sc <- sc + 1
			check <- ders(new, dat, ...)
			if(check[[1]] < logL) {new <- curnt - (1/sc) * step}
			if(check[[1]] >= logL){newL<-check[[1]]
                                                 newg<-check[[2]]}
                        if(check[[1]] >= logL) break
			if(sc == 10) {cat("Step halving not effective, try new starting values", fill = T)
				     stop()}                        
		}
                dist1<-abs(newL-logL)
		dist2 <- (sum((new - curnt)^2))^0.5
		if(crit1 > dist1){c1<-1
#            cat("Convergence criterion of ",crit1,"met for change in log likelihood",fill=T)
		}
		if(crit2 > dist2){c2<-1
#            cat("Convergence criterion of ",crit2," met for change in estimates",fill=T)
		}
                if(sum(crit3 > abs(newg)) == nump){ c3<-1
#            cat("Convergence criterion of ",crit3," met for sum of derivatives",fill=T)
                }
                if(c1+c2+c3==3)	break
		curnt <- new
	}
# 	cat("", fill = T)
#         cat(" ",fill=T)
	final <- ders(new, dat, ...)
        flogL<-final[[1]]
        fgrad<-final[[2]]
        fInf<--1*solve(final[[3]])
# 	cat("Final Estimates Are: ", new, fill = T)
# 	cat("", fill = T)
# 	cat("Final Log Likelihood: ", flogL, fill = T)
# 	cat("", fill = T)
# 	cat("Value of Gradient at Convergence:", fill = T)
# 	cat(fgrad, fill = T)
# 	cat("", fill = T)
# 	cat("Inverse Observed Information: ", fill = T)
# 	cat("(i.e., Inverse of Negative Hessian)", fill = T)
# 	cat("", fill = T)
# 	print(fInf)
	res<-list(par=new,
            loglik=flogL,
            inv_neg_hess=fInf)
        return(res)
}

