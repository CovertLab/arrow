#include <Python.h>
#include <numpy/arrayobject.h>
#include "obsidian.h"

typedef struct {
  PyObject_HEAD
  PyObject * x_attr;
  int reactions_length;
  int substrates_length;
  double * stoichiometry;
  double * rates;

  double * reactants_lengths;
  double * reactants_indexes;
  double * reactants;
  double * reactions;

  double * dependencies_lengths;
  double * dependencies_indexes;
  double * dependencies;
} ObsidianObject;

static PyTypeObject Obsidian_Type;

#define ObsidianObject_Check(v) (Py_TYPE(v) == &Obsidian_Type)

static ObsidianObject *
newObsidianObject(int reactions_length,
                  int substrates_length,
                  double * stoichiometry,
                  double * rates,

                  double * reactants_lengths,
                  double * reactants_indexes,
                  double * reactants,
                  double * reactions,

                  double * dependencies_lengths,
                  double * dependencies_indexes,
                  double * dependencies)
{
  ObsidianObject * self;
  self = PyObject_New(ObsidianObject, &Obsidian_Type);

  if (self == NULL)
    return NULL;

  self->x_attr = NULL;
  self->reactions_length = reactions_length;
  self->substrates_length = substrates_length;
  self->stoichiometry = stoichiometry;
  self->rates = rates;

  self->reactants_lengths = reactants_lengths;
  self->reactants_indexes = reactants_indexes;
  self->reactants = reactants;
  self->reactions = reactions;
  self->dependencies_lengths = dependencies_lengths;
  self->dependencies_indexes = dependencies_indexes;
  self->dependencies = dependencies;

  return self;
}

static void
Obsidian_dealloc(ObsidianObject *self)
{
  Py_XDECREF(self->x_attr);
  PyObject_Del(self);
}

static PyObject *
Obsidian_demo(ObsidianObject *self, PyObject *args)
{
  if (!PyArg_ParseTuple(args, ":demo"))
    return NULL;

  print_array(self->rates, self->reactions_length);

  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *
Obsidian_reactions_length(ObsidianObject *self, PyObject *args)
{
  if (!PyArg_ParseTuple(args, ":reactions_length"))
    return NULL;

  return Py_BuildValue("i", self->reactions_length);
}

static PyObject *
Obsidian_substrates_length(ObsidianObject *self, PyObject *args)
{
  if (!PyArg_ParseTuple(args, ":substrates_length"))
    return NULL;

  return Py_BuildValue("i", self->substrates_length);
}

static PyMethodDef Obsidian_methods[] = {
  {"demo", (PyCFunction) Obsidian_demo,  METH_VARARGS, PyDoc_STR("demo() -> None")},
  {"reactions_length", (PyCFunction) Obsidian_reactions_length,  METH_VARARGS, PyDoc_STR("number of reactions this system is capable of.")},
  {"substrates_length", (PyCFunction) Obsidian_substrates_length,  METH_VARARGS, PyDoc_STR("the number of substrates the reactions in this system operate upon.")},
  {NULL, NULL} // sentinel
};

static PyObject *
Obsidian_getattro(ObsidianObject *self, PyObject *name)
{
  if (self->x_attr != NULL) {
    PyObject *v = PyDict_GetItem(self->x_attr, name);

    if (v != NULL) {
      Py_INCREF(v);
      return v;
    }
  }

  return PyObject_GenericGetAttr((PyObject *)self, name);
}

static int
Obsidian_setattr(ObsidianObject *self, const char *name, PyObject *v)
{
  if (self->x_attr == NULL) {
    self->x_attr = PyDict_New();

    if (self->x_attr == NULL)
      return -1;
  }
  if (v == NULL) {
    int rv = PyDict_DelItemString(self->x_attr, name);

    if (rv < 0)
      PyErr_SetString(PyExc_AttributeError,
                      "delete non-existing Obsidian attribute");
    return rv;
  }
  else
    return PyDict_SetItemString(self->x_attr, name, v);
}

static PyTypeObject Obsidian_Type = {
  /* The ob_type field must be initialized in the module init function
   * to be portable to Windows without using C++. */
  PyVarObject_HEAD_INIT(NULL, 0)
  "obsidianmodule.Obsidian",        /*tp_name*/
  sizeof(ObsidianObject),           /*tp_basicsize*/
  0,                                /*tp_itemsize*/
  /* methods */
  (destructor) Obsidian_dealloc,    /*tp_dealloc*/
  0,                                /*tp_print*/
  (getattrfunc) 0,                  /*tp_getattr*/
  (setattrfunc) Obsidian_setattr,   /*tp_setattr*/
  0,                                /*tp_reserved*/
  0,                                /*tp_repr*/
  0,                                /*tp_as_number*/
  0,                                /*tp_as_sequence*/
  0,                                /*tp_as_mapping*/
  0,                                /*tp_hash*/
  0,                                /*tp_call*/
  0,                                /*tp_str*/
  (getattrofunc) Obsidian_getattro, /*tp_getattro*/
  0,                                /*tp_setattro*/
  0,                                /*tp_as_buffer*/
  Py_TPFLAGS_DEFAULT,               /*tp_flags*/
  0,                                /*tp_doc*/
  0,                                /*tp_traverse*/
  0,                                /*tp_clear*/
  0,                                /*tp_richcompare*/
  0,                                /*tp_weaklistoffset*/
  0,                                /*tp_iter*/
  0,                                /*tp_iternext*/
  Obsidian_methods,                 /*tp_methods*/
  0,                                /*tp_members*/
  0,                                /*tp_getset*/
  0,                                /*tp_base*/
  0,                                /*tp_dict*/
  0,                                /*tp_descr_get*/
  0,                                /*tp_descr_set*/
  0,                                /*tp_dictoffset*/
  0,                                /*tp_init*/
  0,                                /*tp_alloc*/
  0,                                /*tp_new*/
  0,                                /*tp_free*/
  0,                                /*tp_is_gc*/
};


// https://dfm.io/posts/python-c-extensions/

static PyObject *
_print_array(PyObject * self, PyObject * args) {
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

  int value = print_array(pr, pr_length);

  free(pr);

  return Py_BuildValue("i", value);
}

static PyObject *
array_for(PyObject * array_obj) {
  PyObject * array = PyArray_FROM_OTF(array_obj, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY);

  if (array == NULL) {
    Py_XDECREF(array);
    return NULL;
  }

  return array;
}

static double *
data_for(PyObject * array_obj) {
  PyObject * array = array_for(array_obj);
  return (double *) PyArray_DATA(array);
}

static PyObject *
_invoke_obsidian(PyObject * self, PyObject * args) {
  ObsidianObject * obsidian;
  PyObject * stoichiometry_obj,
    * rates_obj,

    * reactants_lengths_obj,
    * reactants_indexes_obj,
    * reactants_obj,
    * reactions_obj,

    * dependencies_lengths_obj,
    * dependencies_indexes_obj,
    * dependencies_obj;

  if (!PyArg_ParseTuple(args,
                        "OOOOOOOOO",
                        &stoichiometry_obj,
                        &rates_obj,

                        &reactants_lengths_obj,
                        &reactants_indexes_obj,
                        &reactants_obj,
                        &reactions_obj,

                        &dependencies_lengths_obj,
                        &dependencies_indexes_obj,
                        &dependencies_obj))
    return NULL;

  // import the stoichiometric_matrix as a 2d numpy array
  PyObject * stoichiometry_array = array_for(stoichiometry_obj);
  int reactions_length = (int) PyArray_DIM(stoichiometry_array, 0);
  int substrates_length = (int) PyArray_DIM(stoichiometry_array, 1);
  double * stoichiometry = (double *) PyArray_DATA(stoichiometry_array);

  printf("dims - %d x %d\n", reactions_length, substrates_length);
  print_array(stoichiometry, reactions_length * substrates_length);

  // import the rates as a 1d numpy array
  PyObject * rates_array = array_for(rates_obj);
  double * rates = (double *) PyArray_DATA(rates_array);

  print_array(rates, reactions_length);

  double * reactants_lengths = data_for(reactants_lengths_obj);
  double * reactants_indexes = data_for(reactants_indexes_obj);
  double * reactants = data_for(reactants_obj);
  double * reactions = data_for(reactions_obj);

  double * dependencies_lengths = data_for(dependencies_lengths_obj);
  double * dependencies_indexes = data_for(dependencies_indexes_obj);
  double * dependencies = data_for(dependencies_obj);

  // create the obsidian object
  obsidian = newObsidianObject(reactions_length,
                               substrates_length,
                               stoichiometry,
                               rates,

                               reactants_lengths,
                               reactants_indexes,
                               reactants,
                               reactions,

                               dependencies_lengths,
                               dependencies_indexes,
                               dependencies);

  if (obsidian == NULL) {
    return NULL;
  }

  return (PyObject *) obsidian;
}

static PyObject *
_evolve(PyObject * self, PyObject * args) {
  return self;
}

static PyMethodDef ObsidianMethods[] = {
  {"print_array", _print_array, METH_VARARGS, "print array of floats"},
  {"obsidian", _invoke_obsidian, METH_VARARGS, "create new obsidian object"},
  {NULL, NULL, 0, NULL}
};

PyMODINIT_FUNC initobsidian(void) {
  (void) Py_InitModule("obsidian", ObsidianMethods);
  import_array();
}

int main(int argc, char** argv) {
  Py_SetProgramName(argv[0]);
  Py_Initialize();

  initobsidian();
}
