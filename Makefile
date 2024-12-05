.PHONY: clean, compile, dist
.DEFAULT_GOAL := compile

clean:
	### Files for older versions of stochastic_arrow are in arrow or stochastic_arrow folder
	rm -rf arrow/arrowhead*.so arrow/arrowhead.c arrow/arrowhead.html
	rm -rf stochastic_arrow/arrowhead*.so stochastic_arrow/arrowhead.c stochastic_arrow/arrowhead.html build/ dist/ MANIFEST .pytest_cache/ stochastic_arrow.egg-info/
	### Newer versions of stochastic_arrow are in src/stochastic_arrow folder
	rm -rf src/stochastic_arrow/arrowhead*.so src/stochastic_arrow/arrowhead.c src/stochastic_arrow/arrowhead.html src/stochastic_arrow.egg-info/
	find . -name "*.pyc" -delete
	find . -name "__pycache__" -delete

compile:
	USE_CYTHON=1 python -m pip install -e .

dist:
	USE_CYTHON=1 python -m build --sdist
