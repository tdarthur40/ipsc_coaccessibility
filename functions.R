# obtained from https://stackoverflow.com/questions/43720260/how-to-extract-p-values-from-lmekin-objects-in-coxme-package
extract_coxme_table <- function (mod){
    beta <- mod$coefficients$fixed
    nvar <- length(beta)
    nfrail <- nrow(mod$var) - nvar
    se <- sqrt(diag(mod$var)[nfrail + 1:nvar])
    z<- beta/se
    p<- pchisq((z)^2, 1, lower.tail = FALSE)
    table=data.frame(cbind(beta,se,z,p))
    return(table)
}

psize = function(height, width) 
{
    options(repr.plot.height=height, repr.plot.width=width)
}

my.invnorm = function(x)
{
    res = rank(x)
    res = qnorm(res/(length(res)+0.5))
    return(res)
}

transform_standard_normal = function(df)
{
    data_valid_expressed_full_qn = normalize.quantiles(as.matrix(df), copy=FALSE)

    input_mat = as.data.frame(t(apply(t(data_valid_expressed_full_qn), 2, my.invnorm)))
    
    return(input_mat)
}


add_rownames = function(x) # add rownames to fread
{
	rownames(x) = x[,1]
	x[,1]       = NULL
	return(x)
}

