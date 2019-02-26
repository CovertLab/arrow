#include <Python.h>
#include <numpy/arrayobject.h>
#include "mersenne.h"
#include "obsidian.h"

// The code in this file acts as a bridge between Python (the caller) and Obsidian
// (the called). Most of this is to reify an object that will hold the references to
// the C values that will be used to perform the actual Gillespie algorithm in pure C.
// The rest is translation of the incoming numpy arrays into C double * and int64_t *. 

// The actual Gillespie algorithm is performed in pure C in obsidian.c, with a header
// in obsidian.h, and some effort has been made to keep it clean of any reference or
// knowledge of python, which is restricted to this file.


// Convert an arbitrary PyObject * into another PyObject * that represents a numpy
// array. The reference to this array is transferred to the caller, who is responsible
// for decrementing the reference once they are finished with it.
static PyObject *
array_for(PyObject *array_obj, int npy_type) {
  PyObject *array = PyArray_FROM_OTF(array_obj, npy_type, NPY_ARRAY_IN_ARRAY);

  if (array == NULL) {
    Py_XDECREF(array);
    return NULL;
  }

  return array;
}

// The definition of the data fields for an Obsidian object as a C struct
typedef struct {
  PyObject_HEAD
  PyObject *x_attr;

  int random_seed;
  MTState *random_state;

  int reactions_count;
  int substrates_count;
  int64_t *stoichiometry;
  double *rates;

  int64_t *reactants_lengths;
  int64_t *reactants_indexes;
  int64_t *reactants;
  int64_t *reactions;

  int64_t *dependencies_lengths;
  int64_t *dependencies_indexes;
  int64_t *dependencies;

  int64_t *substrates_lengths;
  int64_t *substrates_indexes;
  int64_t *substrates;
} ObsidianObject;

// Declaring a new python type for Obsidian
static PyTypeObject Obsidian_Type;
#define ObsidianObject_Check(v) (Py_TYPE(v) == &Obsidian_Type)

// Accept all the information needed to construct a new Obsidian object
// and return a reference to it
static ObsidianObject *
newObsidianObject(int random_seed,
                  int reactions_count,
                  int substrates_count,
                  int64_t *stoichiometry,
                  double *rates,

                  int64_t *reactants_lengths,
                  int64_t *reactants_indexes,
                  int64_t *reactants,
                  int64_t *reactions,

                  int64_t *dependencies_lengths,
                  int64_t *dependencies_indexes,
                  int64_t *dependencies,

                  int64_t *substrates_lengths,
                  int64_t *substrates_indexes,
                  int64_t *substrates)
{
  ObsidianObject *self;
  self = PyObject_New(ObsidianObject, &Obsidian_Type);

  if (self == NULL)
    return NULL;

  self->x_attr = NULL;

  // set up mersenne twister state from the provided seed
  self->random_seed = random_seed;

  MTState *random_state = malloc(sizeof (MTState));
  if (random_state == NULL) {
    return NULL;
  }

  seed(random_state, random_seed);
  self->random_state = random_state;

  self->reactions_count = reactions_count;
  self->substrates_count = substrates_count;
  self->stoichiometry = stoichiometry;
  self->rates = rates;

  self->reactants_lengths = reactants_lengths;
  self->reactants_indexes = reactants_indexes;
  self->reactants = reactants;
  self->reactions = reactions;

  self->dependencies_lengths = dependencies_lengths;
  self->dependencies_indexes = dependencies_indexes;
  self->dependencies = dependencies;

  self->substrates_lengths = substrates_lengths;
  self->substrates_indexes = substrates_indexes;
  self->substrates = substrates;

  return self;
}

// Provide a means for deallocating the Obsidian object
static void
Obsidian_dealloc(ObsidianObject *self)
{
  free(self->random_state);

  Py_XDECREF(self->x_attr);
  PyObject_Del(self);
}

// A simple demo function that prints the reaction rates to stdout
static PyObject *
Obsidian_demo(ObsidianObject *self, PyObject *args)
{
  if (!PyArg_ParseTuple(args, ":demo"))
    return NULL;

  print_array(self->rates, self->reactions_count);

  Py_INCREF(Py_None);
  return Py_None;
}

// Obtain the number of reactions for this system
static PyObject *
Obsidian_reactions_count(ObsidianObject *self, PyObject *args)
{
  if (!PyArg_ParseTuple(args, ":reactions_count"))
    return NULL;

  return Py_BuildValue("i", self->reactions_count);
}

// Find the number of substrates this system operates upon
static PyObject *
Obsidian_substrates_count(ObsidianObject *self, PyObject *args)
{
  if (!PyArg_ParseTuple(args, ":substrates_count"))
    return NULL;

  return Py_BuildValue("i", self->substrates_count);
}

// The main entry point into the system, this function accepts all of the arguments
// from python for running the Gillespie algorithm with the stoichiometry and rates
// initialized earlier on the provided initial state for the given duration.
static PyObject *
Obsidian_evolve(ObsidianObject *self, PyObject *args)
{
  // The variables that will be extracted from python
  double duration;
  PyObject *state_obj;

  // Obtain the arguments as the references declared above
  if (!PyArg_ParseTuple(args, "dO", &duration, &state_obj))
    return NULL;

  // Pull the int64_t * data out of the state numpy array
  PyObject *state_array = array_for(state_obj, NPY_INT64);
  int64_t *state = (int64_t *) PyArray_DATA(state_array);

  // Invoke the actual algorithm with all of the required information
  evolve_result result = evolve(self->random_state,
                                self->reactions_count,
                                self->substrates_count,
                                self->stoichiometry,
                                self->rates,
                                self->reactants_lengths,
                                self->reactants_indexes,
                                self->reactants,
                                self->reactions,
                                self->dependencies_lengths,
                                self->dependencies_indexes,
                                self->dependencies,
                                self->substrates_lengths,
                                self->substrates_indexes,
                                self->substrates,
                                duration,
                                state);
  
  if (result.steps == -1) {
    Py_XDECREF(state_array);
    return NULL;
  }

  // Declare containers for the results
  int64_t steps[1];
  steps[0] = result.steps;

  int64_t substrates[1];
  substrates[0] = self->substrates_count;

  // Create new python numpy arrays from the raw C results
  PyArrayObject *time_obj = (PyArrayObject *) PyArray_SimpleNew(1, steps, NPY_DOUBLE);
  memcpy(PyArray_DATA(time_obj), result.time, (sizeof (result.time[0])) * result.steps);
  PyArrayObject *events_obj = (PyArrayObject *) PyArray_SimpleNew(1, steps, NPY_INT64);
  memcpy(PyArray_DATA(events_obj), result.events, (sizeof (result.events[0])) * result.steps);
  PyArrayObject *outcome_obj = (PyArrayObject *) PyArray_SimpleNew(1, substrates, NPY_INT64);
  memcpy(PyArray_DATA(outcome_obj), result.outcome, (sizeof (result.outcome[0])) * self->substrates_count);

  free(result.time);
  free(result.events);
  free(result.outcome);

  // Decrement the reference to the state array now that we are done with it
  Py_XDECREF(state_array);

  // Construct the return value that will be ultimately visible to python
  return Py_BuildValue("iOOO",
                       result.steps,
                       time_obj,
                       events_obj,
                       outcome_obj);
}

// Declare the various methods that an Obsidian object will have and align them with
// their respective C definitions
static PyMethodDef Obsidian_methods[] = {
  {"demo", (PyCFunction) Obsidian_demo,  METH_VARARGS, PyDoc_STR("demo() -> None")},
  {"reactions_count", (PyCFunction) Obsidian_reactions_count,  METH_VARARGS, PyDoc_STR("number of reactions this system is capable of.")},
  {"substrates_count", (PyCFunction) Obsidian_substrates_count,  METH_VARARGS, PyDoc_STR("the number of substrates the reactions in this system operate upon.")},
  {"evolve", (PyCFunction) Obsidian_evolve,  METH_VARARGS, PyDoc_STR("evolve the system for the given duration and given initial state.")},  
  {NULL, NULL} // sentinel
};

// Provide a means for reading attributes from an Obsidian object
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

// Provide a means for setting attributes in an Obsidian object
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

// Create the object table for an Obsidian object as defined by python internals
static PyTypeObject Obsidian_Type = {
  PyVarObject_HEAD_INIT(NULL, 0)
  "obsidianmodule.Obsidian",        // tp_name
  sizeof(ObsidianObject),           // tp_basicsize
  0,                                // tp_itemsize
  //  methods 
  (destructor) Obsidian_dealloc,    // tp_dealloc
  0,                                // tp_print
  (getattrfunc) 0,                  // tp_getattr
  (setattrfunc) Obsidian_setattr,   // tp_setattr
  0,                                // tp_reserved
  0,                                // tp_repr
  0,                                // tp_as_number
  0,                                // tp_as_sequence
  0,                                // tp_as_mapping
  0,                                // tp_hash
  0,                                // tp_call
  0,                                // tp_str
  (getattrofunc) Obsidian_getattro, // tp_getattro
  0,                                // tp_setattro
  0,                                // tp_as_buffer
  Py_TPFLAGS_DEFAULT,               // tp_flags
  0,                                // tp_doc
  0,                                // tp_traverse
  0,                                // tp_clear
  0,                                // tp_richcompare
  0,                                // tp_weaklistoffset
  0,                                // tp_iter
  0,                                // tp_iternext
  Obsidian_methods,                 // tp_methods
  0,                                // tp_members
  0,                                // tp_getset
  0,                                // tp_base
  0,                                // tp_dict
  0,                                // tp_descr_get
  0,                                // tp_descr_set
  0,                                // tp_dictoffset
  0,                                // tp_init
  0,                                // tp_alloc
  0,                                // tp_new
  0,                                // tp_free
  0,                                // tp_is_gc
};

// Print a python array of floats
static PyObject *
_print_array(PyObject *self, PyObject *args) {
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

  int index;
  for (index = 0; index < pr_length; index++) {
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

// Accept the necessary information to construct the Obsidian object.
static PyObject *
_invoke_obsidian(PyObject *self, PyObject *args) {
  ObsidianObject *obsidian;

  int random_seed;
  PyObject *stoichiometry_obj,
    *rates_obj,

    *reactants_lengths_obj,
    *reactants_indexes_obj,
    *reactants_obj,
    *reactions_obj,

    *dependencies_lengths_obj,
    *dependencies_indexes_obj,
    *dependencies_obj,

    *substrates_lengths_obj,
    *substrates_indexes_obj,
    *substrates_obj;

  if (!PyArg_ParseTuple(args,
                        "iOOOOOOOOOOOO",
                        &random_seed,
                        &stoichiometry_obj,
                        &rates_obj,

                        &reactants_lengths_obj,
                        &reactants_indexes_obj,
                        &reactants_obj,
                        &reactions_obj,

                        &dependencies_lengths_obj,
                        &dependencies_indexes_obj,
                        &dependencies_obj,

                        &substrates_lengths_obj,
                        &substrates_indexes_obj,
                        &substrates_obj))
    return NULL;

  // Import the 2d stoichiometric_matrix as a 1d numpy array
  PyObject *stoichiometry_array = array_for(stoichiometry_obj, NPY_INT64);
  int reactions_count = (int) PyArray_DIM(stoichiometry_array, 0);
  int substrates_count = (int) PyArray_DIM(stoichiometry_array, 1);
  int64_t *stoichiometry = (int64_t *) PyArray_DATA(stoichiometry_array);

  // Import the rates for each reaction
  PyObject *rates_array = array_for(rates_obj, NPY_DOUBLE);
  double *rates = (double *) PyArray_DATA(rates_array);

  // Pull all the precalculated nested arrays
  PyObject *reactants_lengths_array = array_for(reactants_lengths_obj, NPY_INT64);
  int64_t *reactants_lengths = (int64_t *) PyArray_DATA(reactants_lengths_array);
  PyObject *reactants_indexes_array = array_for(reactants_indexes_obj, NPY_INT64);
  int64_t *reactants_indexes = (int64_t *) PyArray_DATA(reactants_indexes_array);
  PyObject *reactants_array = array_for(reactants_obj, NPY_INT64);
  int64_t *reactants = (int64_t *) PyArray_DATA(reactants_array);
  PyObject *reactions_array = array_for(reactions_obj, NPY_INT64);
  int64_t *reactions = (int64_t *) PyArray_DATA(reactions_array);

  PyObject *dependencies_lengths_array = array_for(dependencies_lengths_obj, NPY_INT64);
  int64_t *dependencies_lengths = (int64_t *) PyArray_DATA(dependencies_lengths_array);
  PyObject *dependencies_indexes_array = array_for(dependencies_indexes_obj, NPY_INT64);
  int64_t *dependencies_indexes = (int64_t *) PyArray_DATA(dependencies_indexes_array);
  PyObject *dependencies_array = array_for(dependencies_obj, NPY_INT64);
  int64_t *dependencies = (int64_t *) PyArray_DATA(dependencies_array);

  PyObject *substrates_lengths_array = array_for(substrates_lengths_obj, NPY_INT64);
  int64_t *substrates_lengths = (int64_t *) PyArray_DATA(substrates_lengths_array);
  PyObject *substrates_indexes_array = array_for(substrates_indexes_obj, NPY_INT64);
  int64_t *substrates_indexes = (int64_t *) PyArray_DATA(substrates_indexes_array);
  PyObject *substrates_array = array_for(substrates_obj, NPY_INT64);
  int64_t *substrates = (int64_t *) PyArray_DATA(substrates_array);

  // Create the obsidian object
  obsidian = newObsidianObject(random_seed,
                               reactions_count,
                               substrates_count,
                               stoichiometry,
                               rates,

                               reactants_lengths,
                               reactants_indexes,
                               reactants,
                               reactions,

                               dependencies_lengths,
                               dependencies_indexes,
                               dependencies,

                               substrates_lengths,
                               substrates_indexes,
                               substrates);

  // Clean up all the PyObject * references
  Py_XDECREF(stoichiometry_array);
  Py_XDECREF(rates_array);

  Py_XDECREF(reactants_lengths_array);
  Py_XDECREF(reactants_indexes_array);
  Py_XDECREF(reactants_array);
  Py_XDECREF(reactions_array);

  Py_XDECREF(dependencies_lengths_array);
  Py_XDECREF(dependencies_indexes_array);
  Py_XDECREF(dependencies_array);

  Py_XDECREF(substrates_lengths_array);
  Py_XDECREF(substrates_indexes_array);
  Py_XDECREF(substrates_array);

  // if something went wrong throw an error
  if (obsidian == NULL) {
    return NULL;
  }

  // Return the Obsidian object as a PyObject *
  return (PyObject *) obsidian;
}

// Declare the methods that the Obsidian module will have (as opposed to the Obsidian
// object)
static PyMethodDef ObsidianMethods[] = {
  {"print_array", _print_array, METH_VARARGS, "print array of floats"},
  {"obsidian", _invoke_obsidian, METH_VARARGS, "create new obsidian object"},
  {NULL, NULL, 0, NULL}
};

// Initialize the Obsidian module within the python runtime
PyMODINIT_FUNC initobsidian(void) {
  (void) Py_InitModule("obsidian", ObsidianMethods);
  import_array();
}

// Perform fundamental initialization
int main(int argc, char** argv) {
  Py_SetProgramName(argv[0]);
  Py_Initialize();

  initobsidian();
}
