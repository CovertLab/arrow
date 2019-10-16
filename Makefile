.PHONY: compile, clean

compile:
	python setup.py build_ext --inplace

clean:
	rm -rf arrow/arrowhead*.so arrow/arrowhead.c build/ .pytest_cache/
	find . -name "*.pyc" -delete
	find . -name "__pycache__" -delete
