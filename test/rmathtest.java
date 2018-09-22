import java_Rmath.*;
/* lapacketest.java */
public final class rmathtest {
    private rmathtest() {}
    public static void main(String[] args) {

	/* Random Number Generators */

	//norm_rand	
	System.out.println("norm_rand()");	
	System.out.printf("%.4f\t%.4f\t%.4f\t%.4f\n", jniRmath.norm_rand(), jniRmath.norm_rand(), jniRmath.norm_rand(), jniRmath.norm_rand());
        System.out.println();

	//unif_rand	
	System.out.println("unif_rand()");	
	System.out.printf("%.4f\t%.4f\t%.4f\t%.4f\n", jniRmath.unif_rand(), jniRmath.unif_rand(), jniRmath.unif_rand(), jniRmath.unif_rand());
        System.out.println();

	//exp_rand	
	System.out.println("exp_rand()");	
	System.out.printf("%.4f\t%.4f\t%.4f\t%.4f\n", jniRmath.exp_rand(), jniRmath.exp_rand(), jniRmath.exp_rand(), jniRmath.exp_rand());
        System.out.println();

	/* Normal Distribution */
	double x, mu, sigma, p;
	int give_log, lower_tail, log_p, i_tail;
	double[] cum;
	double[] ccum;


	//dnorm
	System.out.println("dnorm");
	x = 1.96; // x is quantiles
	mu = 0; // mu is means
	sigma = 1; // sigma is standard deviation
	give_log = 0; // if give_log = 0, probability p is not given as log(p); otherwise, p is given as log(p)
	System.out.printf("F(%.4f) = %.4f\n", x, jniRmath.dnorm(x, mu, sigma, give_log));
	System.out.println();

	//pnorm
	System.out.println("pnorm");
	x = 1.96; // x is quantiles
	mu = 0; // mu is means
	sigma = 1; // sigma is standard deviation
	lower_tail = 1; // if lower_tail = 0, probability is P(X>x), otherwise, P(X<=x)
	log_p = 0; // if log_p = 0, probability p is not given as log(p); otherwise, p is given as log(p)
	System.out.printf("p(%.4f) = %.4f\n", x, jniRmath.pnorm(x, mu, sigma, lower_tail, log_p));
	System.out.println();

	//qnorm
	System.out.println("qnorm");
	p = 0.95; // p is probability
	mu = 0; // mu is means
	sigma = 1; // sigma is standard deviation
	lower_tail = 1; // if lower_tail = 0, probability is P(X>x), otherwise, P(X<=x)
	log_p = 0; // if log_p = 0, probability p is not given as log(p); otherwise, p is given as log(p)
	System.out.printf("q(%.4f) = %.4f\n", p, jniRmath.qnorm(p, mu, sigma, lower_tail, log_p));
	System.out.println();

	//rnorm
	System.out.println("rnorm");
	mu = 0; // mu is means
	sigma = 1; // sigma is standard deviation
	System.out.printf("rnorm = %.4f\n",jniRmath.rnorm(mu, sigma));
	System.out.println();

	//pnorm_both
	System.out.println("pnorm_both");
	x = 1.96; // x is quantiles
	cum = new double[] {0}; // one-side probability
	ccum = new double[] {0}; // other-side probability
	i_tail = 0; // in {0,1,2}, means "lower", "upper", "both"; if (lower) return *cum = p(X<=x), if (upper) return *ccum = P(X>x)
	log_p = 0; // if log_p = 0, probability p is not given as log(p); otherwise, p is given as log(p)
	jniRmath.pnorm_both (x, cum, ccum, i_tail, log_p);
	System.out.printf("p(%.4f) = %.4f\t1 - p(%.4f) = %.4f\n", x, cum[0], x, ccum[0]);
	System.out.println();

	/* Uniform Distribution */
	double a, b;

	//dunif
	System.out.println("dunif");
	x = 1.96; // x is quantiles
	a = 0; // a is lower limits of the distribution
	b = 2; // b is upper limits of the distribution
	give_log = 0; // if give_log = 0, probability p is not given as log(p); otherwise, p is given as log(p)
	System.out.printf("F(%.4f) = %.4f\n", x, jniRmath.dunif(x, a, b, give_log));
	System.out.println();

	//punif
	System.out.println("punif");
	x = 0.9; // x is quantiles
	a = 0; // a is lower limits of the distribution
	b = 1; // b is upper limits of the distribution
	lower_tail = 1; // if lower_tail = 0, probability is P(X>x), otherwise, P(X<=x)
	log_p = 0; // if log_p = 0, probability p is not given as log(p); otherwise, p is given as log(p)
	System.out.printf("p(%.4f) = %.4f\n", x, jniRmath.punif(x, a, b, lower_tail, log_p));
	System.out.println();

	//qunif
	System.out.println("qunif");
	p = 0.95; // p is quantiles
	a = 0; // a is lower limits of the distribution
	b = 1; // b is upper limits of the distribution
	lower_tail = 1; // if lower_tail = 0, probability is P(X>x), otherwise, P(X<=x)
	log_p = 0; // if log_p = 0, probability p is not given as log(p); otherwise, p is given as log(p)
	System.out.printf("q(%.4f) = %.4f\n", p, jniRmath.qunif(p, a, b, lower_tail, log_p));
	System.out.println();

	//runif
	System.out.println("runif");
	a = 0; // a is lower limits of the distribution
	b = 1; // b is upper limits of the distribution
	System.out.printf("runif = %.4f\n",jniRmath.runif(a, b));
	System.out.println();

	/* Gamma Distribution */
	double shape, scale;

	//dgamma
	System.out.println("dgamma");
	x = 2; // x is quantiles
	shape = 5; // shape is shape parameter
	scale = 1; // scale is scale parameter
	give_log = 0; // if give_log = 0, probability p is not given as log(p); otherwise, p is given as log(p)
	System.out.printf("F(%.4f) = %.4f\n", x, jniRmath.dgamma(x, shape, scale, give_log));
	System.out.println();

	//pgamma
	System.out.println("pgamma");
	x = 0.9; // x is quantiles
	shape = 5; // shape is shape parameter
	scale = 1; // scale is scale parameter
	lower_tail = 1; // if lower_tail = 0, probability is P(X>x), otherwise, P(X<=x)
	log_p = 0; // if log_p = 0, probability p is not given as log(p); otherwise, p is given as log(p)
	System.out.printf("p(%.4f) = %.4f\n", x, jniRmath.pgamma(x, shape, scale, lower_tail, log_p));
	System.out.println();

	//qgamma
	System.out.println("qgamma");
	p = 0.95; // p is quantiles
	shape = 5; // shape is shape parameter
	scale = 1; // scale is scale parameter
	lower_tail = 1; // if lower_tail = 0, probability is P(X>x), otherwise, P(X<=x)
	log_p = 0; // if log_p = 0, probability p is not given as log(p); otherwise, p is given as log(p)
	System.out.printf("q(%.4f) = %.4f\n", p, jniRmath.qgamma(p, shape, scale, lower_tail, log_p));
	System.out.println();

	//rgamma
	System.out.println("rgamma");
	shape = 5; // shape is shape parameter
	scale = 1; // scale is scale parameter
	System.out.printf("rgamma = %.4f\n",jniRmath.rgamma(shape, scale));
	System.out.println();

	/* Beta Distribution */

	//dbeta
	System.out.println("dbeta");
	x = 0.4; // x is quantiles
	a = 2; // a is non-negative parameter of the beta distribution
	b = 2; // b is non-negative parameter of the beta distribution
	give_log = 0; // if give_log = 0, probability p is not given as log(p); otherwise, p is given as log(p)
	System.out.printf("F(%.4f) = %.4f\n", x, jniRmath.dbeta(x, a, b, give_log));
	System.out.println();

	//pbeta
	System.out.println("pbeta");
	x = 0.4; // x is quantiles
	a = 2; // a is non-negative parameter of the beta distribution
	b = 2; // b is non-negative parameter of the beta distribution
	lower_tail = 1; // if lower_tail = 0, probability is P(X>x), otherwise, P(X<=x)
	log_p = 0; // if log_p = 0, probability p is not given as log(p); otherwise, p is given as log(p)
	System.out.printf("p(%.4f) = %.4f\n", x, jniRmath.pbeta(x, a, b, lower_tail, log_p));
	System.out.println();

	//qbeta
	System.out.println("qbeta");
	p = 0.95; // p is quantiles
	a = 2; // a is non-negative parameter of the beta distribution
	b = 2; // b is non-negative parameter of the beta distribution
	lower_tail = 1; // if lower_tail = 0, probability is P(X>x), otherwise, P(X<=x)
	log_p = 0; // if log_p = 0, probability p is not given as log(p); otherwise, p is given as log(p)
	System.out.printf("q(%.4f) = %.4f\n", p, jniRmath.qbeta(p, a, b, lower_tail, log_p));
	System.out.println();

	//rbeta
	System.out.println("rbeta");
	a = 2; // a is non-negative parameter of the beta distribution
	b = 2; // b is non-negative parameter of the beta distribution
	System.out.printf("rbeta = %.4f\n",jniRmath.rbeta(a, b));
	System.out.println();

	/* Lognormal Distribution */
	double meanlog, sdlog;

	//dlnorm
	System.out.println("dlnorm");
	x = 0.4; // x is quantiles
	meanlog = 0; // meanlog is mean of the distribution on the log scale
	sdlog = 1; // sdlog is standard deviation of the distribution on the log scale
	give_log = 0; // if give_log = 0, probability p is not given as log(p); otherwise, p is given as log(p)
	System.out.printf("F(%.4f) = %.4f\n", x, jniRmath.dlnorm(x, meanlog, sdlog, give_log));
	System.out.println();

	//plnorm
	System.out.println("plnorm");
	x = 0.4; // x is quantiles
	meanlog = 0; // meanlog is mean of the distribution on the log scale
	sdlog = 1; // sdlog is standard deviation of the distribution on the log scale	
	lower_tail = 1; // if lower_tail = 0, probability is P(X>x), otherwise, P(X<=x)
	log_p = 0; // if log_p = 0, probability p is not given as log(p); otherwise, p is given as log(p)
	System.out.printf("p(%.4f) = %.4f\n", x, jniRmath.plnorm(x, meanlog, sdlog, lower_tail, log_p));
	System.out.println();

	//qlnorm
	System.out.println("qlnorm");
	p = 0.95; // x is quantiles
	meanlog = 0; // meanlog is mean of the distribution on the log scale
	sdlog = 1; // sdlog is standard deviation of the distribution on the log scale	
	lower_tail = 1; // if lower_tail = 0, probability is P(X>x), otherwise, P(X<=x)
	log_p = 0; // if log_p = 0, probability p is not given as log(p); otherwise, p is given as log(p)
	System.out.printf("q(%.4f) = %.4f\n", p, jniRmath.qlnorm(p, meanlog, sdlog, lower_tail, log_p));
	System.out.println();

	//rlnorm
	System.out.println("rlnorm");
	meanlog = 0; // meanlog is mean of the distribution on the log scale
	sdlog = 1; // sdlog is standard deviation of the distribution on the log scale	
	System.out.printf("rlnorm = %.4f\n",jniRmath.rlnorm(meanlog, sdlog));
	System.out.println();

	/* Chi-squared Distribution */
	double df;

	//dchisq
	System.out.println("dchisq");
	x = 0.4; // x is quantiles
	df = 3;	// df is degrees of freedom
	give_log = 0; // if give_log = 0, probability p is not given as log(p); otherwise, p is given as log(p)
	System.out.printf("F(%.4f) = %.4f\n", x, jniRmath.dchisq(x, df, give_log));
	System.out.println();

	//pchisq
	System.out.println("pchisq");
	x = 0.4; // x is quantiles
	df = 3;	// df is degrees of freedom
	lower_tail = 1; // if lower_tail = 0, probability is P(X>x), otherwise, P(X<=x)
	log_p = 0; // if log_p = 0, probability p is not given as log(p); otherwise, p is given as log(p)
	System.out.printf("p(%.4f) = %.4f\n", x, jniRmath.pchisq(x, df, lower_tail, log_p));
	System.out.println();

	//qchisq
	System.out.println("qchisq");
	p = 0.95; // p is quantiles
	df = 3;	// df is degrees of freedom
	lower_tail = 1; // if lower_tail = 0, probability is P(X>x), otherwise, P(X<=x)
	log_p = 0; // if log_p = 0, probability p is not given as log(p); otherwise, p is given as log(p)
	System.out.printf("q(%.4f) = %.4f\n", p, jniRmath.qchisq(p, df, lower_tail, log_p));
	System.out.println();

	//rchisq
	System.out.println("rchisq");
	df = 3;	// df is degrees of freedom
	System.out.printf("rchisq = %.4f\n",jniRmath.rchisq(df));
	System.out.println();

	/* Non-central Chi-squared Distribution */
	double ncp;

	//dnchisq
	System.out.println("dnchisq");
	x = 0.4; // x is quantiles
	df = 3;	// df is degrees of freedom
	ncp = 1; // ncp is non-centrality parameter (non-negative)
	give_log = 0; // if give_log = 0, probability p is not given as log(p); otherwise, p is given as log(p)
	System.out.printf("F(%.4f) = %.4f\n", x, jniRmath.dnchisq(x, df, ncp, give_log));
	System.out.println();

	//pnchisq
	System.out.println("pnchisq");
	x = 0.4; // x is quantiles
	df = 3;	// df is degrees of freedom
	ncp = 1; // ncp is non-centrality parameter (non-negative)
	lower_tail = 1; // if lower_tail = 0, probability is P(X>x), otherwise, P(X<=x)
	log_p = 0; // if log_p = 0, probability p is not given as log(p); otherwise, p is given as log(p)
	System.out.printf("p(%.4f) = %.4f\n", x, jniRmath.pnchisq(x, df, ncp, lower_tail, log_p));
	System.out.println();

	//qnchisq
	System.out.println("qnchisq");
	p = 0.95; // p is quantiles
	df = 3;	// df is degrees of freedom
	ncp = 1; // ncp is non-centrality parameter (non-negative)
	lower_tail = 1; // if lower_tail = 0, probability is P(X>x), otherwise, P(X<=x)
	log_p = 0; // if log_p = 0, probability p is not given as log(p); otherwise, p is given as log(p)
	System.out.printf("q(%.4f) = %.4f\n", p, jniRmath.qnchisq(p, df, ncp, lower_tail, log_p));
	System.out.println();

	//rnchisq
	System.out.println("rnchisq");
	df = 3;	// df is degrees of freedom
	ncp = 1; // ncp is non-centrality parameter (non-negative)
	System.out.printf("rnchisq = %.4f\n",jniRmath.rnchisq(df, ncp));
	System.out.println();


	/* F Distibution */
	double df1, df2;

	//df
	System.out.println("df");
	x = 0.4; // x is quantiles
	df1 = 3; // df1 is degrees of freedom
	df2 = 1; // df2 is degrees of freedom
	give_log = 0; // if give_log = 0, probability p is not given as log(p); otherwise, p is given as log(p)
	System.out.printf("F(%.4f) = %.4f\n", x, jniRmath.df(x, df1, df2, give_log));
	System.out.println();

	//pf
	System.out.println("pf");
	x = 0.4; // x is quantiles
	df1 = 3; // df1 is degrees of freedom
	df2 = 1; // df2 is degrees of freedom
	lower_tail = 1; // if lower_tail = 0, probability is P(X>x), otherwise, P(X<=x)
	log_p = 0; // if log_p = 0, probability p is not given as log(p); otherwise, p is given as log(p)
	System.out.printf("p(%.4f) = %.4f\n", x, jniRmath.pf(x, df1, df2, lower_tail, log_p));
	System.out.println();

	//qf
	System.out.println("qf");
	p = 0.95; // p is quantiles
	df1 = 3; // df1 is degrees of freedom
	df2 = 1; // df2 is degrees of freedom
	lower_tail = 1; // if lower_tail = 0, probability is P(X>x), otherwise, P(X<=x)
	log_p = 0; // if log_p = 0, probability p is not given as log(p); otherwise, p is given as log(p)
	System.out.printf("q(%.4f) = %.4f\n", p, jniRmath.qf(p, df1, df2, lower_tail, log_p));
	System.out.println();

	//rf
	System.out.println("rf");
	df1 = 3; // df1 is degrees of freedom
	df2 = 1; // df2 is degrees of freedom
	System.out.printf("rf = %.4f\n",jniRmath.rf(df1, df2));
	System.out.println();


	/* Student t Distibution */

	//dt
	System.out.println("dt");
	x = 0.4; // x is quantiles
	df = 3; // df is degrees of freedom
	give_log = 0; // if give_log = 0, probability p is not given as log(p); otherwise, p is given as log(p)
	System.out.printf("F(%.4f) = %.4f\n", x, jniRmath.dt(x, df, give_log));
	System.out.println();

	//pt
	System.out.println("pt");
	x = 0.4; // x is quantiles
	df = 3; // df is degrees of freedom
	lower_tail = 1; // if lower_tail = 0, probability is P(X>x), otherwise, P(X<=x)
	log_p = 0; // if log_p = 0, probability p is not given as log(p); otherwise, p is given as log(p)
	System.out.printf("p(%.4f) = %.4f\n", x, jniRmath.pt(x, df, lower_tail, log_p));
	System.out.println();

	//qt
	System.out.println("qt");
	p = 0.95; // p is quantiles
	df = 3; // df is degrees of freedom
	lower_tail = 1; // if lower_tail = 0, probability is P(X>x), otherwise, P(X<=x)
	log_p = 0; // if log_p = 0, probability p is not given as log(p); otherwise, p is given as log(p)
	System.out.printf("q(%.4f) = %.4f\n", p, jniRmath.qt(p, df, lower_tail, log_p));
	System.out.println();

	//rt
	System.out.println("rt");
	df = 3; // df is degrees of freedom
	System.out.printf("rt = %.4f\n",jniRmath.rt(df));
	System.out.println();


	/* Binomial Distribution */
	double n, pr;

	//dbinom
	System.out.println("dbinom");
	x = 3; // x is quantiles
	n = 5; // n is number of trials (zero or more)
	p = 0.3; // p is probability of success on each trial
	give_log = 0; // if give_log = 0, probability p is not given as log(p); otherwise, p is given as log(p)
	System.out.printf("F(%.4f) = %.4f\n", x, jniRmath.dbinom(x, n, p, give_log));
	System.out.println();

	//pbinom
	System.out.println("pbinom");
	x = 3; // x is quantiles
	n = 5; // n is number of trials (zero or more)
	p = 0.3; // p is probability of success on each trial
	lower_tail = 1; // if lower_tail = 0, probability is P(X>x), otherwise, P(X<=x)
	log_p = 0; // if log_p = 0, probability p is not given as log(p); otherwise, p is given as log(p)
	System.out.printf("p(%.4f) = %.4f\n", x, jniRmath.pbinom(x, n, p, lower_tail, log_p));
	System.out.println();

	//qbinom
	System.out.println("qbinom");
	p = 0.95; // p is quantiles
	n = 5; // n is number of trials (zero or more)
	lower_tail = 1; // if lower_tail = 0, probability is P(X>x), otherwise, P(X<=x)
	log_p = 0; // if log_p = 0, probability p is not given as log(p); otherwise, p is given as log(p)
	System.out.printf("q(%.4f) = %.4f\n", p, jniRmath.qbinom(p, n, p, lower_tail, log_p));
	System.out.println();
	
	//rbinom
	System.out.println("rbinom");
	n = 5; // n is number of trials (zero or more)
	p = 0.3; // p is probability of success on each trial
	System.out.printf("rbinom = %.4f\n",jniRmath.rbinom(n, p));
	System.out.println();


	/* Multnomial Distribution */
	double[] prob;
	int k;
	int[] rN;

	//rmultinom
	System.out.println("rmultinom");
	n = 5; // n is number of trials (zero or more)
	prob = new double[] {0.2,0.2,0.25,0.15,0.2}; // prob specifies the probability for the k classes
	k = 5; // k = length of prob
	rN = new int[]{0,0,0,0,0}; // rN is return vector, length = k, rN[j] ~ Bin (n, prob[j])
	jniRmath.rmultinom(5, prob, k, rN);
	System.out.printf("rmultinom = %d\t%d\t%d\t%d\t%d\n", rN[0], rN[1], rN[2], rN[3], rN[4]);
	System.out.println();


	/* Cauchy Distribution */
	double location;

	//dcauchy
	System.out.println("dcauchy");
	x = 1; // x is quantiles
	location = 0; // location is location parameter
	scale = 1; // scale is scale parameter
	give_log = 0; // if give_log = 0, probability p is not given as log(p); otherwise, p is given as log(p)
	System.out.printf("F(%.4f) = %.4f\n", x, jniRmath.dcauchy(x, location, scale, give_log));
	System.out.println();

	//pcauchy
	System.out.println("pcauchy");
	x = 1; // x is quantiles
	location = 0; // location is location parameter
	scale = 1; // scale is scale parameter
	lower_tail = 1; // if lower_tail = 0, probability is P(X>x), otherwise, P(X<=x)
	log_p = 0; // if log_p = 0, probability p is not given as log(p); otherwise, p is given as log(p)
	System.out.printf("p(%.4f) = %.4f\n", x, jniRmath.pcauchy(x, location, scale, lower_tail, log_p));
	System.out.println();

	//qcauchy
	System.out.println("qcauchy");
	p = 0.95; // p is quantiles
	location = 0; // location is location parameter
	scale = 1; // scale is scale parameter
	lower_tail = 1; // if lower_tail = 0, probability is P(X>x), otherwise, P(X<=x)
	log_p = 0; // if log_p = 0, probability p is not given as log(p); otherwise, p is given as log(p)
	System.out.printf("q(%.4f) = %.4f\n", p, jniRmath.qcauchy(p, location, scale, lower_tail, log_p));
	System.out.println();
	
	//rcauchy
	System.out.println("rcauchy");
	location = 0; // location is location parameter
	scale = 1; // scale is scale parameter
	System.out.printf("rcauchy = %.4f\n",jniRmath.rcauchy(location, scale));
	System.out.println();


	/* Exponential Distribution */

	//dexp
	System.out.println("dexp");
	x = 1; // x is quantiles
	scale = 1; // scale is scale parameter
	give_log = 0; // if give_log = 0, probability p is not given as log(p); otherwise, p is given as log(p)
	System.out.printf("F(%.4f) = %.4f\n", x, jniRmath.dexp(x, scale, give_log));
	System.out.println();

	//pexp
	System.out.println("pexp");
	x = 1; // x is quantiles
	scale = 1; // scale is scale parameter
	lower_tail = 1; // if lower_tail = 0, probability is P(X>x), otherwise, P(X<=x)
	log_p = 0; // if log_p = 0, probability p is not given as log(p); otherwise, p is given as log(p)
	System.out.printf("p(%.4f) = %.4f\n", x, jniRmath.pexp(x, scale, lower_tail, log_p));
	System.out.println();

	//qexp
	System.out.println("qexp");
	p = 0.95; // p is quantiles
	scale = 1; // scale is scale parameter
	lower_tail = 1; // if lower_tail = 0, probability is P(X>x), otherwise, P(X<=x)
	log_p = 0; // if log_p = 0, probability p is not given as log(p); otherwise, p is given as log(p)
	System.out.printf("q(%.4f) = %.4f\n", p, jniRmath.qexp(p, scale, lower_tail, log_p));
	System.out.println();
	
	//rexp
	System.out.println("rexp");
	scale = 1; // scale is scale parameter
	System.out.printf("rexp = %.4f\n",jniRmath.rexp(scale));
	System.out.println();


	/* Geometric Distribution */
	double pro;

	//dgeom
	System.out.println("dgeom");
	x = 1; // x is quantiles
	p = 0.3; // p is probability of success in each trial
	give_log = 0; // if give_log = 0, probability p is not given as log(p); otherwise, p is given as log(p)
	System.out.printf("F(%.4f) = %.4f\n", x, jniRmath.dgeom(x, p, give_log));
	System.out.println();

	//pgeom
	System.out.println("pgeom");
	x = 1; // x is quantiles
	p = 0.3; // p is probability of success in each trial
	lower_tail = 1; // if lower_tail = 1, probability is P(X<=x), otherwise, P(X>x)
	log_p = 0; // if log_p = 1, probability p is given as log(p)
	System.out.printf("p(%.4f) = %.4f\n", x, jniRmath.pgeom(x, p, lower_tail, log_p));
	System.out.println();

	//qgeom
	System.out.println("qgeom");
	p = 0.95; // p is quantiles
	pro = 0.3; // pro is probability of success in each trial
	lower_tail = 1; // if lower_tail = 0, probability is P(X>x), otherwise, P(X<=x)
	log_p = 0; // if log_p = 0, probability p is not given as log(p); otherwise, p is given as log(p)
	System.out.printf("q(%.4f) = %.4f\n", p, jniRmath.qgeom(p, pro, lower_tail, log_p));
	System.out.println();
	
	//rgeom
	System.out.println("rgeom");
	p = 0.3; // p is probability of success in each trial
	System.out.printf("rgeom = %.4f\n",jniRmath.rgeom(p));
	System.out.println();


	/* Hypergeometric Distibution */
	double r;

	//dhyper
	System.out.println("dhyper");
	x = 1; // x is quantiles representing the number of successes in our sample
	r = 4; // r is the number of successes in a sequence
	b = 6; // b is the number of failures in a sequence
	n = 3; // n is the number of items we sample, n < r + b
	give_log = 0; // if give_log = 0, probability p is not given as log(p); otherwise, p is given as log(p)
	System.out.printf("F(%.4f) = %.4f\n", x, jniRmath.dhyper(x, r, b, n, give_log));
	System.out.println();

	//phyper
	System.out.println("phyper");
	x = 1; // x is quantiles
	r = 4; // r is the number of successes in a sequence
	b = 6; // b is the number of failures in a sequence
	n = 3; // n is the number of items in a sequence, n = r + b
	lower_tail = 1; // if lower_tail = 0, probability is P(X>x), otherwise, P(X<=x)
	log_p = 0; // if log_p = 0, probability p is not given as log(p); otherwise, p is given as log(p)
	System.out.printf("p(%.4f) = %.4f\n", x, jniRmath.phyper(x, r, b, n, lower_tail, log_p));
	System.out.println();

	//qhyper
	System.out.println("qhyper");
	p = 0.95; // p is quantiles
	r = 4; // r is the number of successes in a sequence
	b = 6; // b is the number of failures in a sequence
	n = 3; // n is the number of items in a sequence, n = r + b
	lower_tail = 1; // if lower_tail = 0, probability is P(X>x), otherwise, P(X<=x)
	log_p = 0; // if log_p = 0, probability p is not given as log(p); otherwise, p is given as log(p)
	System.out.printf("q(%.4f) = %.4f\n", p, jniRmath.qhyper(p, r, b, n, lower_tail, log_p));
	System.out.println();
	
	//rhyper
	System.out.println("rhyper");
	r = 4; // r is the number of successes in a sequence
	b = 6; // b is the number of failures in a sequence
	n = 3; // n is the number of items in a sequence, n = r + b
	System.out.printf("rhyper = %.4f\n",jniRmath.rhyper(r, b, n));
	System.out.println();

	/* Negative Binomial Distribution */
	double size;

	//dnbinom
	System.out.println("dnbinom");
	x = 1; // x is quantiles representing the number of successes in our sample
	size = 2; // size is target for number of successful trials, or dispersion parameter (the shape parameter of the gamma mixing distribution)
	pro = 0.3; // pro is the probability of success in each trail  
	give_log = 0; // if give_log = 0, probability p is not given as log(p); otherwise, p is given as log(p)
	System.out.printf("F(%.4f) = %.4f\n", x, jniRmath.dnbinom(x, size, pro, give_log));
	System.out.println();

	//pnbinom
	System.out.println("pnbinom");
	x = 1; // x is quantiles
	size = 2; // size is target for number of successful trials, or dispersion parameter (the shape parameter of the gamma mixing distribution)
	pro = 0.3; // pro is the probability of success in each trail  
	lower_tail = 1; // if lower_tail = 0, probability is P(X>x), otherwise, P(X<=x)
	log_p = 0; // if log_p = 0, probability p is not given as log(p); otherwise, p is given as log(p)
	System.out.printf("p(%.4f) = %.4f\n", x, jniRmath.pnbinom(x, size, pro, lower_tail, log_p));
	System.out.println();

	//qnbinom
	System.out.println("qnbinom");
	p = 0.95; // p is quantiles
	size = 2; // size is target for number of successful trials, or dispersion parameter (the shape parameter of the gamma mixing distribution)
	pro = 0.3; // pro is the probability of success in each trail  
	lower_tail = 1; // if lower_tail = 0, probability is P(X>x), otherwise, P(X<=x)
	log_p = 0; // if log_p = 0, probability p is not given as log(p); otherwise, p is given as log(p)
	System.out.printf("q(%.4f) = %.4f\n", p, jniRmath.qnbinom(p, size, pro, lower_tail, log_p));
	System.out.println();
	
	//rbinom
	System.out.println("rbinom");
	size = 2; // size is target for number of successful trials, or dispersion parameter (the shape parameter of the gamma mixing distribution)
	pro = 0.3; // pro is the probability of success in each trail  
	System.out.printf("rbinom = %.4f\n",jniRmath.rbinom(size, pro));
	System.out.println();


	/* Poisson Distribution */
	double lambda;

	//dpois
	System.out.println("dpois");
	x = 1; // x is quantiles representing the number of successes in our sample
	lambda = 1; // vector of means
	give_log = 0; // if give_log = 0, probability p is not given as log(p); otherwise, p is given as log(p)
	System.out.printf("F(%.4f) = %.4f\n", x, jniRmath.dpois(x, lambda, give_log));
	System.out.println();

	//ppois
	System.out.println("ppois");
	x = 1; // x is quantiles
	lambda = 1; // vector of means
	lower_tail = 1; // if lower_tail = 0, probability is P(X>x), otherwise, P(X<=x)
	log_p = 0; // if log_p = 0, probability p is not given as log(p); otherwise, p is given as log(p)
	System.out.printf("p(%.4f) = %.4f\n", x, jniRmath.ppois(x, lambda, lower_tail, log_p));
	System.out.println();

	//qpois
	System.out.println("qpois");
	p = 0.95; // p is quantiles
	lambda = 1; // vector of means
	lower_tail = 1; // if lower_tail = 0, probability is P(X>x), otherwise, P(X<=x)
	log_p = 0; // if log_p = 0, probability p is not given as log(p); otherwise, p is given as log(p)
	System.out.printf("q(%.4f) = %.4f\n", p, jniRmath.qpois(p, lambda, lower_tail, log_p));
	System.out.println();
	
	//rpois
	System.out.println("rpois");
	lambda = 1; // vector of means
	System.out.printf("rpois = %.4f\n",jniRmath.rpois(lambda));
	System.out.println();


	/* Weibull Distribution */

	//dweibull
	System.out.println("dweibull");
	x = 1; // x is quantiles representing the number of successes in our sample
	shape = 1; // shape is shape parameter
	scale = 1; // scale is scale parameter	
	give_log = 0; // if give_log = 0, probability p is not given as log(p); otherwise, p is given as log(p)
	System.out.printf("F(%.4f) = %.4f\n", x, jniRmath.dweibull(x, shape, scale, give_log));
	System.out.println();

	//pweibull
	System.out.println("pweibull");
	x = 1; // x is quantiles
	shape = 1; // shape is shape parameter
	scale = 1; // scale is scale parameter	
	lower_tail = 1; // if lower_tail = 0, probability is P(X>x), otherwise, P(X<=x)
	log_p = 0; // if log_p = 0, probability p is not given as log(p); otherwise, p is given as log(p)
	System.out.printf("p(%.4f) = %.4f\n", x, jniRmath.pweibull(x, shape, scale, lower_tail, log_p));
	System.out.println();

	//qweibull
	System.out.println("qweibull");
	p = 0.95; // p is quantiles
	shape = 1; // shape is shape parameter
	scale = 1; // scale is scale parameter	
	lower_tail = 1; // if lower_tail = 0, probability is P(X>x), otherwise, P(X<=x)
	log_p = 0; // if log_p = 0, probability p is not given as log(p); otherwise, p is given as log(p)
	System.out.printf("q(%.4f) = %.4f\n", p, jniRmath.qweibull(p, shape, scale, lower_tail, log_p));
	System.out.println();
	
	//rweibull
	System.out.println("rweibull");
	lambda = 1; // vector of means
	System.out.printf("rweibull = %.4f\n",jniRmath.rweibull(shape, scale));
	System.out.println();


	/* Logistic Distribution */

	//dlogis
	System.out.println("dlogis");
	x = 1; // x is quantiles representing the number of successes in our sample
	location = 0; // location is location parameter
	scale = 1; // scale is scale parameter	
	give_log = 0; // if give_log = 0, probability p is not given as log(p); otherwise, p is given as log(p)
	System.out.printf("F(%.4f) = %.4f\n", x, jniRmath.dlogis(x, location, scale, give_log));
	System.out.println();

	//plogis
	System.out.println("plogis");
	x = 1; // x is quantiles
	location = 0; // location is location parameter
	scale = 1; // scale is scale parameter	
	lower_tail = 1; // if lower_tail = 0, probability is P(X>x), otherwise, P(X<=x)
	log_p = 0; // if log_p = 0, probability p is not given as log(p); otherwise, p is given as log(p)
	System.out.printf("p(%.4f) = %.4f\n", x, jniRmath.plogis(x, location, scale, lower_tail, log_p));
	System.out.println();

	//qlogis
	System.out.println("qlogis");
	p = 0.95; // p is quantiles
	location = 0; // location is location parameter
	scale = 1; // scale is scale parameter	
	lower_tail = 1; // if lower_tail = 0, probability is P(X>x), otherwise, P(X<=x)
	log_p = 0; // if log_p = 0, probability p is not given as log(p); otherwise, p is given as log(p)
	System.out.printf("q(%.4f) = %.4f\n", p, jniRmath.qlogis(p, location, scale, lower_tail, log_p));
	System.out.println();
	
	//rlogis
	System.out.println("rlogis");
	location = 0; // location is location parameter
	scale = 1; // scale is scale parameter	
	System.out.printf("rlogis = %.4f\n",jniRmath.rlogis(location, scale));
	System.out.println();


	/* Non-central Beta Distribution */

	//dnbeta
	System.out.println("dnbeta");
	x = 1; // x is quantiles representing the number of successes in our sample
	a = 1; // a is non-negative parameter of the beta distribution
	b = 1; // b is non-negative parameter of the beta distribution
	ncp = 1; // ncp is non-centrality parameter	
	give_log = 0; // if give_log = 0, probability p is not given as log(p); otherwise, p is given as log(p)
	System.out.printf("F(%.4f) = %.4f\n", x, jniRmath.dnbeta(x, a, b, ncp, give_log));
	System.out.println();

	//pnbeta
	System.out.println("pnbeta");
	x = 1; // x is quantiles
	a = 1; // a is non-negative parameter of the beta distribution
	b = 1; // b is non-negative parameter of the beta distribution
	ncp = 1; // ncp is non-centrality parameter	
	lower_tail = 1; // if lower_tail = 0, probability is P(X>x), otherwise, P(X<=x)
	log_p = 0; // if log_p = 0, probability p is not given as log(p); otherwise, p is given as log(p)
	System.out.printf("p(%.4f) = %.4f\n", x, jniRmath.pnbeta(x, a, b, ncp, lower_tail, log_p));
	System.out.println();

	//qnbeta
	System.out.println("qnbeta");
	p = 0.95; // p is quantiles
	a = 1; // a is non-negative parameter of the beta distribution
	b = 1; // b is non-negative parameter of the beta distribution
	ncp = 1; // ncp is non-centrality parameter	
	lower_tail = 1; // if lower_tail = 0, probability is P(X>x), otherwise, P(X<=x)
	log_p = 0; // if log_p = 0, probability p is not given as log(p); otherwise, p is given as log(p)
	System.out.printf("q(%.4f) = %.4f\n", p, jniRmath.qnbeta(p, a, b, ncp, lower_tail, log_p));
	System.out.println();
	
	//rnbeta
	System.out.println("rnbeta");
	a = 1; // a is non-negative parameter of the beta distribution
	b = 1; // b is non-negative parameter of the beta distribution
	ncp = 1; // ncp is non-centrality parameter	
	//System.out.printf("rnbeta = %.4f\n", jniRmath.rnbeta(a, b, ncp));
	System.out.println();


	/* Non-central F Distribution */

	//dnf
	System.out.println("dnf");
	x = 1; // x is quantiles representing the number of successes in our sample
	df1 = 3; // df1 is the first degree of freedom
	df2 = 2; // df2 is the second degree of freedom
	ncp = 1; // ncp is non-centrality parameter	
	give_log = 0; // if give_log = 0, probability p is not given as log(p); otherwise, p is given as log(p)
	System.out.printf("F(%.4f) = %.4f\n", x, jniRmath.dnf(x, df1, df2, ncp, give_log));
	System.out.println();

	//pnf
	System.out.println("pnf");
	x = 1; // x is quantiles
	df1 = 3; // df1 is the first degree of freedom
	df2 = 2; // df2 is the second degree of freedom
	ncp = 1; // ncp is non-centrality parameter	
	lower_tail = 1; // if lower_tail = 0, probability is P(X>x), otherwise, P(X<=x)
	log_p = 0; // if log_p = 0, probability p is not given as log(p); otherwise, p is given as log(p)
	System.out.printf("p(%.4f) = %.4f\n", x, jniRmath.pnf(x, df1, df2, ncp, lower_tail, log_p));
	System.out.println();

	//qnf
	System.out.println("qnf");
	p = 0.95; // p is quantiles
	df1 = 3; // df1 is the first degree of freedom
	df2 = 2; // df2 is the second degree of freedom
	ncp = 1; // ncp is non-centrality parameter	
	lower_tail = 1; // if lower_tail = 0, probability is P(X>x), otherwise, P(X<=x)
	log_p = 0; // if log_p = 0, probability p is not given as log(p); otherwise, p is given as log(p)
	System.out.printf("q(%.4f) = %.4f\n", p, jniRmath.qnf(p, df1, df2, ncp, lower_tail, log_p));
	System.out.println();


	/* Non-central Student t Distribution */

	//dnt
	System.out.println("dnt");
	x = 1; // x is quantiles representing the number of successes in our sample
	df = 3; // df is the first degree of freedom
	ncp = 1; // ncp is non-centrality parameter	
	give_log = 0; // if give_log = 0, probability p is not given as log(p); otherwise, p is given as log(p)
	System.out.printf("F(%.4f) = %.4f\n", x, jniRmath.dnt(x, df, ncp, give_log));
	System.out.println();

	//pnt
	System.out.println("pnt");
	x = 1; // x is quantiles
	df = 3; // df is the first degree of freedom
	ncp = 1; // ncp is non-centrality parameter	
	lower_tail = 1; // if lower_tail = 0, probability is P(X>x), otherwise, P(X<=x)
	log_p = 0; // if log_p = 0, probability p is not given as log(p); otherwise, p is given as log(p)
	System.out.printf("p(%.4f) = %.4f\n", x, jniRmath.pnt(x, df, ncp, lower_tail, log_p));
	System.out.println();

	//qnt
	System.out.println("qnt");
	p = 0.95; // p is quantiles
	df = 3; // df is the first degree of freedom
	ncp = 1; // ncp is non-centrality parameter	
	lower_tail = 1; // if lower_tail = 0, probability is P(X>x), otherwise, P(X<=x)
	log_p = 0; // if log_p = 0, probability p is not given as log(p); otherwise, p is given as log(p)
	System.out.printf("q(%.4f) = %.4f\n", p, jniRmath.qnt(p, df, ncp, lower_tail, log_p));
	System.out.println();


	/* Studentized Range Distribution */
	double rr, cc;

	//ptukey
	System.out.println("ptukey");
	x = 1; // x is quantiles
	rr = 1; // rr is th maximum of studentized ranges
	cc = 5; // cc is means
	df = 10; // df is degrees of freedom
	lower_tail = 1; // if lower_tail = 0, probability is P(X>x), otherwise, P(X<=x)
	log_p = 0; // if log_p = 0, probability p is not given as log(p); otherwise, p is given as log(p)
	System.out.printf("p(%.4f) = %.4f\n", x, jniRmath.ptukey(x, rr, cc, df, lower_tail, log_p));
	System.out.println();

	//qtukey
	System.out.println("qtukey");
	p = 0.95; // p is quantiles
	rr = 1; // rr is th maximum of studentized ranges
	cc = 5; // cc is means
	df = 10; // df is degrees of freedom
	lower_tail = 1; // if lower_tail = 0, probability is P(X>x), otherwise, P(X<=x)
	log_p = 0; // if log_p = 0, probability p is not given as log(p); otherwise, p is given as log(p)
	System.out.printf("q(%.4f) = %.4f\n", p, jniRmath.qtukey(p, rr, cc, df, lower_tail, log_p));
	System.out.println();


	/* Wilcoxon Rank Sum Distribution */
	double m;

	//dwilcox
	System.out.println("dwilcox");
	x = 1; // x is quantiles representing the number of successes in our sample
	m = 5; // m is number of observations in the first sample	
	n = 4; // n is number of observations in the second sample
	give_log = 0; // if give_log = 0, probability p is not given as log(p); otherwise, p is given as log(p)
	System.out.printf("F(%.4f) = %.4f\n", x, jniRmath.dwilcox(x, m, n, give_log));
	System.out.println();

	//pwilcox
	System.out.println("pwilcox");
	x = 1; // x is quantiles
	m = 5; // m is number of observations in the first sample	
	n = 4; // n is number of observations in the second sample
	lower_tail = 1; // if lower_tail = 0, probability is P(X>x), otherwise, P(X<=x)
	log_p = 0; // if log_p = 0, probability p is not given as log(p); otherwise, p is given as log(p)
	System.out.printf("p(%.4f) = %.4f\n", x, jniRmath.pwilcox(x, m, n, lower_tail, log_p));
	System.out.println();

	//qwilcox
	System.out.println("qwilcox");
	p = 0.95; // p is quantiles
	m = 5; // m is number of observations in the first sample	
	n = 4; // n is number of observations in the second sample
	lower_tail = 1; // if lower_tail = 0, probability is P(X>x), otherwise, P(X<=x)
	log_p = 0; // if log_p = 0, probability p is not given as log(p); otherwise, p is given as log(p)
	System.out.printf("q(%.4f) = %.4f\n", p, jniRmath.qwilcox(p, m, n, lower_tail, log_p));
	System.out.println();
	
	//rwilcox
	System.out.println("rwilcox");
	m = 5; // m is number of observations in the first sample	
	n = 4; // n is number of observations in the second sample
	System.out.printf("rwilcox = %.4f\n", jniRmath.rwilcox(m, n));
	System.out.println();


	/* Wilcoxon Signed Rank Distribution */

	//dsignrank
	System.out.println("dsignrank");
	x = 1; // x is quantiles representing the number of successes in our sample
	n = 4; // n is number of observations in the sample
	give_log = 0; // if give_log = 0, probability p is not given as log(p); otherwise, p is given as log(p)
	System.out.printf("F(%.4f) = %.4f\n", x, jniRmath.dsignrank(x, n, give_log));
	System.out.println();

	//psignrank
	System.out.println("psignrank");
	x = 1; // x is quantiles
	n = 4; // n is number of observations in the sample
	lower_tail = 1; // if lower_tail = 0, probability is P(X>x), otherwise, P(X<=x)
	log_p = 0; // if log_p = 0, probability p is not given as log(p); otherwise, p is given as log(p)
	System.out.printf("p(%.4f) = %.4f\n", x, jniRmath.psignrank(x, n, lower_tail, log_p));
	System.out.println();

	//qsignrank
	System.out.println("qsignrank");
	p = 0.95; // p is quantiles
	n = 4; // n is number of observations in the sample
	lower_tail = 1; // if lower_tail = 0, probability is P(X>x), otherwise, P(X<=x)
	log_p = 0; // if log_p = 0, probability p is not given as log(p); otherwise, p is given as log(p)
	System.out.printf("q(%.4f) = %.4f\n", p, jniRmath.qsignrank(p, n, lower_tail, log_p));
	System.out.println();
	
	//rsignrank
	System.out.println("rsignrank");
	n = 4; // n is number of observations in the sample
	System.out.printf("rsignrank = %.4f\n", jniRmath.rsignrank(n));
	System.out.println();


	/* Gamma and Related Functions */
	int[] sgn;
	int kode;
	double[] ans;
	int[] nz;
	int[] ierr;
	double deriv;

	//gammafn
	System.out.println("gammafn");
	x = 3; // x is quantiles
	System.out.printf("gammafn(%.4f) = %.4f\n", x, jniRmath.gammafn(x));
	System.out.println();

	//lgammafn
	System.out.println("lgammafn");
	x = 3; // x is quantiles
	System.out.printf("lgammafn(%.4f) = %.4f\n", x, jniRmath.lgammafn(x));
	System.out.println();

	//lgammafn_sign
	System.out.println("lgammafn_sign");
	x = 3; // x is quantiles
	sgn = new int[]{0}; // address to store the sign of the gamma function when output
	System.out.printf("lgammafn_sign(%.4f) = %.4f\t, sgn = %d\n", x, jniRmath.lgammafn_sign(x, sgn), sgn[0]);
	System.out.println();

	//dpsifn
	System.out.println("dpsifn");
	x = 3; // x is quantiles
	n = 0; // n is first member of the sequence, 0 <= n <= 100
	kode = 1; // k is selection parameter, takes value over {1,2}, for kode = 1, dpsifn returns the scaled derivaties of the psi function; kode = 2 returns the scaled derivatives of the psi function except when n = 0. In this case, ans(1) = -psi(x) + ln(x) is returned
	m = 3; // m is number of members of the sequence, m >= 1
	ans = new double[]{0,0,0}; // ans is output, a vector of length at least m whose first m components contain the sequence of derivatives scaled according to kode
	nz = new int[] {0}; // nz is underflow flag, if nz = 0, a normal return; if nz != 0, underflow, last nz components of ans are set to zero
	ierr = new int[] {0}; // ierr is error flag, if ierr = 0, a normal return; if ierr = 1, input error, no computation; if ierr = 2, overflow, x too small or n+m-1 too large or both; if ierr = 3, error, n too large
	jniRmath.dpsifn(x, (int)n, kode, (int)m, ans, nz, ierr);
	System.out.printf("dpsifn(%.4f, %.4f, %d, %.4f): ans = %.4f %.4f %.4f\n", x, n, kode, m, ans[0], ans[1], ans[2]);
	System.out.println();

	//psigamma
	System.out.println("psigamma");
	x = 3; // x is quantiles
	deriv = 2; // deriv specifies the order of derivative
	System.out.printf("psigamma(%.4f, %.4f) = %.4f\n", x, deriv, jniRmath.psigamma(x, deriv));
	System.out.println();

	//digamma
	System.out.println("digamma");
	x = 3; // x is quantiles
	System.out.printf("digamma(%.4f) = %.4f\n", x, jniRmath.digamma(x));
	System.out.println();

	//trigamma
	System.out.println("trigamma");
	x = 3; // x is quantiles
	System.out.printf("trigamma(%.4f) = %.4f\n", x, jniRmath.trigamma(x));
	System.out.println();

	//tetragamma
	System.out.println("tetragamma");
	x = 3; // x is quantiles
	System.out.printf("tetragamma(%.4f) = %.4f\n", x, jniRmath.tetragamma(x));
	System.out.println();

	//pentagamma
	System.out.println("pentagamma");
	x = 3; // x is quantiles
	System.out.printf("pentagamma(%.4f) = %.4f\n", x, jniRmath.pentagamma(x));
	System.out.println();

	//beta
	System.out.println("beta");
	a = 2; // a is the first shape parameter
	b = 5; // b is the second shape parameter
	System.out.printf("beta(%.4f, %.4f) = %.4f\n", a, b, jniRmath.beta(a, b));
	System.out.println();

	//lbeta
	System.out.println("lbeta");
	a = 2; // a is the first shape parameter
	b = 5; // b is the second shape parameter
	System.out.printf("lbeta(%.4f, %.4f) = %.4f\n", a, b, jniRmath.lbeta(a, b));
	System.out.println();

	//choose
	System.out.println("choose");
	n = 5.0; // n is the n in binomial coefficient C(n,k)
	k = 2; // k is the k in binomial coefficient C(n,k)
	System.out.printf("choose(%.4f, %.4f) = %.4f\n", n, (double)k, jniRmath.choose(n, k));
	System.out.println();

	//lchoose
	System.out.println("lchoose");
	n = 5.0; // n is the n in binomial coefficient C(n,k)
	k = 2; // k is the k in binomial coefficient C(n,k)
	System.out.printf("lchoose(%.4f, %.4f) = %.4f\n", n, (double)k, jniRmath.lchoose(n, k));
	System.out.println();


	/* Bessel Functions */
	double alpha, expo;

	//bessel_i
	System.out.println("bessel_i");
	x = 0.3; // x is quantiles
	alpha = 1; // alpha is order of the corresponding bessel function
	expo = 2; //takes value over {1,2}, if expo = 2, the results are exponentially scaled in order to avoid overflow or underflow, respectively	
	System.out.printf("bessel_i(%.4f, %.4f, %.4f) = %.4f\n", x, alpha, expo, jniRmath.bessel_i(x, alpha, expo));
	System.out.println();

	//bessel_j
	System.out.println("bessel_j");
	x = 0.3; // x is quantiles
	alpha = 1; // alpha is order of the corresponding bessel function
	System.out.printf("bessel_j(%.4f, %.4f) = %.4f\n", x, alpha, jniRmath.bessel_j(x, alpha));
	System.out.println();

	//bessel_k
	System.out.println("bessel_k");
	x = 0.3; // x is quantiles
	alpha = 1; // alpha is order of the corresponding bessel function
	expo = 2; //takes value over {1,2}, if expo = 2, the results are exponentially scaled in order to avoid overflow or underflow, respectively	
	System.out.printf("bessel_k(%.4f, %.4f, %.4f) = %.4f\n", x, alpha, expo, jniRmath.bessel_k(x, alpha, expo));
	System.out.println();

	//bessel_y
	System.out.println("bessel_y");
	x = 0.3; // x is quantiles
	alpha = 1; // alpha is order of the corresponding bessel function
	System.out.printf("bessel_y(%.4f, %.4f) = %.4f\n", x, alpha, jniRmath.bessel_y(x, alpha));
	System.out.println();


    }
    
}







