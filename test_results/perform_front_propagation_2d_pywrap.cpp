/*=================================================================
% perform_front_propagation_2d - perform a Fast Marching front propagation.
%
%   Copyright (c) 2004 Gabriel Peyré
*=================================================================*/
/* Must include python.h before any standard headers */
# include "Python.h"
# include "perform_front_propagation_2d.h"

/* Define the wrapper functions exposed to python (must be static) */
static PyObject* wrap_fmm(PyObject *self, PyObject *args) {
	int n, result;
	/* Python -> C conversion*/
	if (!PyArg_ParseTuple(args, "i",&n))
	return NULL;
	/* Call your function*/
	result = perform_front_propagation_2d(n);
	
	/* C -> Python conversion*/
	return Py_BuildValue("i",result);
	
}

/* Methods table declaring the names of functions exposed to python */
static PyMethodDef ExampleMethods[] = {
	{"perform_front_propagation_2d", wrap_fmm,METH_VARARGS, "Compute fmm"},
	{NULL,NULL,0,NULL}  /* sentinel  */
};

/* Module initialization function called at "import example" */
PyMODINIT_FUNC initexample(void){
	(void) Py_InitModule ("example", ExampleMethods)
}