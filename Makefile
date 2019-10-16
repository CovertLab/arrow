.PHONY: clean, compile, dist

clean:
	rm -rf arrow/arrowhead*.so arrow/arrowhead.c arrow/arrowhead.html build/ dist/ MANIFEST .pytest_cache/ stochastic_arrow.egg-info/
	find . -name "*.pyc" -delete
	find . -name "__pycache__" -delete

compile:
	USE_CYTHON=1 python setup.py build_ext --inplace

dist:
	USE_CYTHON=1 python setup.py sdist bdist_wheel
