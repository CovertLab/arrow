.PHONY: clean, compile, dist
.DEFAULT_GOAL := compile

clean:
	### Files for older versions of stochastic_arrow are in arrow folder
	rm -rf arrow/arrowhead*.so arrow/arrowhead.c arrow/arrowhead.html
	rm -rf stochastic_arrow/arrowhead*.so stochastic_arrow/arrowhead.c stochastic_arrow/arrowhead.html build/ dist/ MANIFEST .pytest_cache/ stochastic_arrow.egg-info/
	find . -name "*.pyc" -delete
	find . -name "__pycache__" -delete

compile:
	USE_CYTHON=1 python setup.py build_ext --inplace

dist:
	### bdist_wheel is disabled on linux since the distribution machinery doesn't
	### yet have a way to specify compatible linux distros.
	USE_CYTHON=1 python setup.py sdist # bdist_wheel
