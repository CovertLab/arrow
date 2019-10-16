.PHONY: clean, compile, dist
.DEFAULT_GOAL := compile

clean:
	rm -rf arrow/arrowhead*.so arrow/arrowhead.c arrow/arrowhead.html build/ dist/ MANIFEST .pytest_cache/ stochastic_arrow.egg-info/
	find . -name "*.pyc" -delete
	find . -name "__pycache__" -delete

compile:
	USE_CYTHON=1 python setup.py build_ext --inplace

dist:
	### bdist_wheel is disabled on linux since the distribution process doesn't
	### yet have a way to specify compatible linux distros.
	# USE_CYTHON=1 python setup.py sdist bdist_wheel
	USE_CYTHON=1 python setup.py sdist
