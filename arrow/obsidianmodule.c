#include <Python.h>

int _print_array(double pr[], int length) {
  for (int index = 0; index < length; index++)
    printf("pr[%d] = %f\n", index, pr[index]);
  return 0;
}

static PyObject * print_array(PyObject * self, PyObject * args) {
  PyObject *float_list;
  int pr_length;
  double *pr;

  if (!PyArg_ParseTuple(args, "O", &float_list))
    return NULL;

  pr_length = PyObject_Length(float_list);
  if (pr_length < 0)
    return NULL;

  pr = (double *) malloc(sizeof(double *) * pr_length);
  if (pr == NULL)
    return NULL;

  for (int index = 0; index < pr_length; index++) {
    PyObject *item;
    item = PyList_GetItem(float_list, index);
    if (!PyFloat_Check(item))
      pr[index] = 0.0;

    pr[index] = PyFloat_AsDouble(item);
  }

  return Py_BuildValue("i", _print_array(pr, pr_length));
}

static PyMethodDef ObsidianMethods[] = {
  {"print_array", print_array, METH_VARARGS, "print array of floats"},
  {NULL, NULL, 0, NULL}
};

PyMODINIT_FUNC initobsidian(void) {
  (void) Py_InitModule("obsidian", ObsidianMethods);
}
