k=5  # dimension of the problem
n=75 # sample size
alpha=runif(k) # value of alpha, here chosen at random
p=rgamma(k,alpha) # pre-simulation of the Dirichlet
y=sample(1:k,n,prob=p/sum(p),rep=TRUE) # Multinomial
x=sum(y==1)
for (i in 2:k) x=c(x,sum(y==i))


alpha <- pred_fit$pars %>% 
  filter(grepl("Int", parameter)) %>% 
  pull(estimate) %>% 
  exp() 
p <- rgamma(length(alpha), alpha)
sample(1:length(p), 20, prob = p, replace = T)
