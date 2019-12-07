context("Generalized Gamma distribution");

require(flexsurv);
require(ggamma);

test_that("dggamma integrates to 1", {
	result = integrate(function(x) dggamma(x, 1, 1, 1), 0, Inf)$value;
	expect_equal(result, 1, tolerance=0.001);

	result = integrate(function(x) dggamma(x, 0.8, 0.2, 1), 0, Inf)$value;
	expect_equal(result, 1, tolerance=0.001);

	result = integrate(function(x) dggamma(x, 0.5, 0.5, 0.2), 0, Inf)$value;
	expect_equal(result, 1, tolerance=0.001);

	result = integrate(function(x) dggamma(x, 0.5, 5, 2), 0, Inf)$value;
	expect_equal(result, 1, tolerance=0.001);

	result = integrate(function(x) dggamma(x, 0.5, 1, 2), 0, Inf)$value;
	expect_equal(result, 1, tolerance=0.001);

	result = integrate(function(x) dggamma(x, 5, 1, 0.5), 0, Inf)$value;
	expect_equal(result, 1, tolerance=0.001);
});

test_that("pggamma correctly reflects dggamma", {
	expected = integrate(function(x) dggamma(x, 1, 1, 1), 0, 10)$value;
	result   = pggamma(10, 1, 1, 1);
	expect_equal(result, expected);

	result   = pggamma(1e-10, 1, 1, 1);
	expect_equal(result, 0);

	result   = pggamma(20, 1, 1, 1);
	expect_equal(result, 1);

	expected = integrate(function(x) dggamma(x, 0.8, 0.5, 1), 0, 10)$value;
	result   = pggamma(10, 0.8, 0.5, 1);
	expect_equal(result, expected);

	result   = pggamma(0, 0.8, 0.5, 1);
	expect_equal(result, 0);

	result   = pggamma(1000, 0.8, 0.5, 1);
	expect_equal(result, 1);

	expected = integrate(function(x) dggamma(x, 1.2, 0.5, 1.5), 0, 10)$value;
	result   = pggamma(10, 1.2, 0.5, 1.5);
	expect_equal(result, expected, tolerance=1e-7);

	result   = pggamma(0, 1.2, 0.5, 1.5);
	expect_equal(result, 0);

	result   = pggamma(1000, 1.2, 0.5, 1.5);
	expect_equal(result, 1);
});

test_that("qggamma correctly inverts pggamma", {
	a = 1; b = 1; k = 1;
	expected = c(1, 10);
	expect_equal(qggamma(pggamma(expected, a, b, k), a, b, k), expected);

	a = 0.5; b = 1.5; k = 1.5;
	expected = c(0.5, 1, 3);
	expect_equal(qggamma(pggamma(expected, a, b, k), a, b, k), expected);

	a = 1.5; b = 0.5; k = 3;
	expected = c(0.5, 1, 3);
	expect_equal(qggamma(pggamma(expected, a, b, k), a, b, k), expected);
});

test_that("random number generation", {
	set.seed(72);

	a = 1; b = 1; k = 1;
	r = rggamma(100000, a, b, k);
	e = ecdf(r);
	x = seq(0.001, 10, length=10000);
	maxError = max(abs(e(x) - pggamma(x, a, b, k)));
	expect_lte(maxError, 1e-2);

	a = 0.3; b = 1.5; k = 1.5;
	r = rggamma(100000, a, b, k);
	e = ecdf(r);
	x = seq(0.001, 10, length=10000);
	maxError = max(abs(e(x) - pggamma(x, a, b, k)));
	expect_lte(maxError, 1e-2);

	a = 1.5; b = 5; k = 0.3;
	r = rggamma(100000, a, b, k);
	e = ecdf(r);
	x = seq(0.001, 10, length=10000);
	maxError = max(abs(e(x) - pggamma(x, a, b, k)));
	expect_lte(maxError, 1e-2);

	a = 3; b = 1.8; k = 0.5;
	r = rggamma(100000, a, b, k);
	e = ecdf(r);
	x = seq(0.001, 10, length=10000);
	maxError = max(abs(e(x) - pggamma(x, a, b, k)));
	expect_lte(maxError, 1e-2);

	set.seed(NULL);
});

test_that("ggamma equals gengamma.orig from flexsurv", {
	x = seq(0.001, 10, length=100);
	q = seq(0.001, 0.999, length=10);
	
	a = 1; b = 1; k = 1;
	expect_equal(dgengamma.orig(x, b, a, k), dggamma(x, a, b, k));
	expect_equal(pgengamma.orig(x, b, a, k), pggamma(x, a, b, k));
	expect_equal(qgengamma.orig(q, b, a, k), qggamma(q, a, b, k));

	a = 0.3; b = 1.5; k = 1.5;
	expect_equal(dgengamma.orig(x, b, a, k), dggamma(x, a, b, k));
	expect_equal(pgengamma.orig(x, b, a, k), pggamma(x, a, b, k));
	expect_equal(qgengamma.orig(q, b, a, k), qggamma(q, a, b, k));

	a = 1.5; b = 5; k = 0.3;
	expect_equal(dgengamma.orig(x, b, a, k), dggamma(x, a, b, k));
	expect_equal(pgengamma.orig(x, b, a, k), pggamma(x, a, b, k));
	expect_equal(qgengamma.orig(q, b, a, k), qggamma(q, a, b, k));

	a = 3; b = 1.8; k = 0.5;
	expect_equal(dgengamma.orig(x, b, a, k), dggamma(x, a, b, k));
	expect_equal(pgengamma.orig(x, b, a, k), pggamma(x, a, b, k));
	expect_equal(qgengamma.orig(q, b, a, k), qggamma(q, a, b, k));
});
