#include <Python.h>
#include <numpy/arrayobject.h>
#include "obsidian.h"

/* static PyObject *ErrorObject; */

/* typedef struct { */
/*   PyObject_HEAD */
/*   PyObject *x_attr; */
/* } ObsidianObject; */

/* static PyTypeObject Obsidian_Type; */

/* #define ObsidianObject_Check(v)      (Py_TYPE(v) == &Obsidian_Type) */

/* static ObsidianObject * */
/* newObsidianObject(PyObject *arg) */
/* { */
/*   ObsidianObject *self; */
/*   self = PyObject_New(ObsidianObject, &Obsidian_Type); */

/*   if (self == NULL) */
/*     return NULL; */

/*   self->x_attr = NULL; */
/*   return self; */
/* } */

/* static void */
/* Obsidian_dealloc(ObsidianObject *self) */
/* { */
/*   Py_XDECREF(self->x_attr); */
/*   PyObject_Del(self); */
/* } */

/* static PyObject * */
/* Obsidian_demo(ObsidianObject *self, PyObject *args) */
/* { */
/*   if (!PyArg_ParseTuple(args, ":demo")) */
/*     return NULL; */

/*   Py_INCREF(Py_None); */
/*   return Py_None; */
/* } */

/* static PyMethodDef Obsidian_methods[] = { */
/*   {"demo", (PyCFunction) Obsidian_demo,  METH_VARARGS, PyDoc_STR("demo() -> None")}, */
/*   {NULL, NULL}           /\* sentinel *\/ */
/* }; */

/* static PyObject * */
/* Obsidian_getattro(ObsidianObject *self, PyObject *name) */
/* { */
/*   if (self->x_attr != NULL) { */
/*     PyObject *v = PyDict_GetItem(self->x_attr, name); */

/*     if (v != NULL) { */
/*       Py_INCREF(v); */
/*       return v; */
/*     } */
/*   } */

/*   return PyObject_GenericGetAttr((PyObject *)self, name); */
/* } */

/* static int */
/* Obsidian_setattr(ObsidianObject *self, const char *name, PyObject *v) */
/* { */
/*   if (self->x_attr == NULL) { */
/*     self->x_attr = PyDict_New(); */

/*     if (self->x_attr == NULL) */
/*       return -1; */
/*   } */
/*   if (v == NULL) { */
/*     int rv = PyDict_DelItemString(self->x_attr, name); */

/*     if (rv < 0) */
/*       PyErr_SetString(PyExc_AttributeError, */
/*                       "delete non-existing Obsidian attribute"); */
/*     return rv; */
/*   } */
/*   else */
/*     return PyDict_SetItemString(self->x_attr, name, v); */
/* } */

/* static PyTypeObject Obsidian_Type = { */
/*   /\* The ob_type field must be initialized in the module init function */
/*    * to be portable to Windows without using C++. *\/ */
/*   PyVarObject_HEAD_INIT(NULL, 0) */
/*   "xxmodule.Obsidian",              /\*tp_name*\/ */
/*   sizeof(ObsidianObject),           /\*tp_basicsize*\/ */
/*   0,                                /\*tp_itemsize*\/ */
/*   /\* methods *\/ */
/*   (destructor) Obsidian_dealloc,    /\*tp_dealloc*\/ */
/*   0,                                /\*tp_print*\/ */
/*   (getattrfunc) 0,                  /\*tp_getattr*\/ */
/*   (setattrfunc) Obsidian_setattr,   /\*tp_setattr*\/ */
/*   0,                                /\*tp_reserved*\/ */
/*   0,                                /\*tp_repr*\/ */
/*   0,                                /\*tp_as_number*\/ */
/*   0,                                /\*tp_as_sequence*\/ */
/*   0,                                /\*tp_as_mapping*\/ */
/*   0,                                /\*tp_hash*\/ */
/*   0,                                /\*tp_call*\/ */
/*   0,                                /\*tp_str*\/ */
/*   (getattrofunc) Obsidian_getattro, /\*tp_getattro*\/ */
/*   0,                                /\*tp_setattro*\/ */
/*   0,                                /\*tp_as_buffer*\/ */
/*   Py_TPFLAGS_DEFAULT,               /\*tp_flags*\/ */
/*   0,                                /\*tp_doc*\/ */
/*   0,                                /\*tp_traverse*\/ */
/*   0,                                /\*tp_clear*\/ */
/*   0,                                /\*tp_richcompare*\/ */
/*   0,                                /\*tp_weaklistoffset*\/ */
/*   0,                                /\*tp_iter*\/ */
/*   0,                                /\*tp_iternext*\/ */
/*   Obsidian_methods,                 /\*tp_methods*\/ */
/*   0,                                /\*tp_members*\/ */
/*   0,                                /\*tp_getset*\/ */
/*   0,                                /\*tp_base*\/ */
/*   0,                                /\*tp_dict*\/ */
/*   0,                                /\*tp_descr_get*\/ */
/*   0,                                /\*tp_descr_set*\/ */
/*   0,                                /\*tp_dictoffset*\/ */
/*   0,                                /\*tp_init*\/ */
/*   0,                                /\*tp_alloc*\/ */
/*   0,                                /\*tp_new*\/ */
/*   0,                                /\*tp_free*\/ */
/*   0,                                /\*tp_is_gc*\/ */
/* }; */




static PyObject * _print_array(PyObject * self, PyObject * args) {
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

  return Py_BuildValue("i", print_array(pr, pr_length));
}

static PyMethodDef ObsidianMethods[] = {
  {"print_array", _print_array, METH_VARARGS, "print array of floats"},
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
