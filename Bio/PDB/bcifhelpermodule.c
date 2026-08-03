#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <math.h>
#include <stdint.h>
#include <string.h>

static int
integer_unpack_u8(Py_buffer *in_view, Py_buffer *out_view)
{
    Py_ssize_t in_size = in_view->shape[0];
    Py_ssize_t out_size = out_view->shape[0];
    Py_ssize_t in_index = 0;
    Py_ssize_t out_index = 0;

    uint8_t *in_data = in_view->buf;
    uint32_t *out_data = out_view->buf;

    while (in_index < in_size) {
        uint32_t sum = in_data[in_index];

        if (sum == UINT8_MAX) {
            while (in_index + 1 < in_size) {
                in_index += 1;
                sum += in_data[in_index];

                if (in_data[in_index] != UINT8_MAX) {
                    break;
                }
            }
            if (in_data[in_index] == UINT8_MAX) {
                return -2;
            }
        }

        if (out_index >= out_size) {
            return -1;
        }
        out_data[out_index] = sum;
        in_index += 1;
        out_index += 1;
    }
    return out_index == out_size ? 0 : -3;
}

static int
integer_unpack_u16(Py_buffer *in_view, Py_buffer *out_view)
{
    Py_ssize_t in_size = in_view->shape[0];
    Py_ssize_t out_size = out_view->shape[0];
    Py_ssize_t in_index = 0;
    Py_ssize_t out_index = 0;

    uint16_t *in_data = in_view->buf;
    uint32_t *out_data = out_view->buf;

    while (in_index < in_size) {
        uint32_t sum = in_data[in_index];

        if (sum == UINT16_MAX) {
            while (in_index + 1 < in_size) {
                in_index += 1;
                sum += in_data[in_index];

                if (in_data[in_index] != UINT16_MAX) {
                    break;
                }
            }
            if (in_data[in_index] == UINT16_MAX) {
                return -2;
            }
        }

        if (out_index >= out_size) {
            return -1;
        }
        out_data[out_index] = sum;
        in_index += 1;
        out_index += 1;
    }
    return out_index == out_size ? 0 : -3;
}

static int
integer_unpack_i8(Py_buffer *in_view, Py_buffer *out_view)
{
    Py_ssize_t in_size = in_view->shape[0];
    Py_ssize_t out_size = out_view->shape[0];
    Py_ssize_t in_index = 0;
    Py_ssize_t out_index = 0;

    int8_t *in_data = in_view->buf;
    int32_t *out_data = out_view->buf;

    while (in_index < in_size) {
        int32_t sum = in_data[in_index];

        if (sum == INT8_MAX || sum == INT8_MIN) {
            while (in_index + 1 < in_size) {
                in_index += 1;
                sum += in_data[in_index];

                if (in_data[in_index] != INT8_MAX && in_data[in_index] != INT8_MIN) {
                    break;
                }
            }
            if (in_data[in_index] == INT8_MAX ||
                in_data[in_index] == INT8_MIN) {
                return -2;
            }
        }

        if (out_index >= out_size) {
            return -1;
        }
        out_data[out_index] = sum;
        in_index += 1;
        out_index += 1;
    }
    return out_index == out_size ? 0 : -3;
}

static int
integer_unpack_i16(Py_buffer *in_view, Py_buffer *out_view)
{
    Py_ssize_t in_size = in_view->shape[0];
    Py_ssize_t out_size = out_view->shape[0];
    Py_ssize_t in_index = 0;
    Py_ssize_t out_index = 0;

    int16_t *in_data = in_view->buf;
    int32_t *out_data = out_view->buf;

    while (in_index < in_size) {
        int32_t sum = in_data[in_index];

        if (sum == INT16_MAX || sum == INT16_MIN) {
            while (in_index + 1 < in_size) {
                in_index += 1;
                sum += in_data[in_index];

                if (in_data[in_index] != INT16_MAX && in_data[in_index] != INT16_MIN) {
                    break;
                }
            }
            if (in_data[in_index] == INT16_MAX ||
                in_data[in_index] == INT16_MIN) {
                return -2;
            }
        }

        if (out_index >= out_size) {
            return -1;
        }
        out_data[out_index] = sum;
        in_index += 1;
        out_index += 1;
    }
    return out_index == out_size ? 0 : -3;
}

/* Return the type code of a PEP 3118 format string, or '\0' if it is not a
   single type code in native byte order.

   The unpack loops above read and write through native-width pointers and
   ignore the declared byte order entirely, so a buffer in the opposite order
   would be decoded to the wrong values rather than rejected. A format string
   may carry a leading byte-order character, and NumPy supplies one whenever
   the array's order is not the platform's: on a big-endian machine the
   little-endian dtypes in Bio/PDB/binary_cif.py arrive here as ">H" and the
   like. Accept the orders that mean native and refuse the rest, so that the
   caller gets an error instead of silently wrong coordinates. */
static char
native_type_code(const char *format)
{
    if (format == NULL) {
        return '\0';
    }
    switch (format[0]) {
        case '@':
        case '=':
#if PY_LITTLE_ENDIAN
        case '<':
#else
        case '>':
        case '!':
#endif
            format++;
            break;
        default:
            break;
    }
    if (format[0] == '\0' || format[1] != '\0') {
        return '\0';
    }
    return format[0];
}

static PyObject *
integer_unpack(PyObject *self, PyObject *args)
{
    PyObject *in = NULL;
    PyObject *out = NULL;

    if (!PyArg_ParseTuple(args, "OO", &in, &out)) {
        return NULL;
    }

    Py_buffer in_view, out_view;
    const int flags = PyBUF_ND | PyBUF_FORMAT;

    if (PyObject_GetBuffer(in, &in_view, flags) != 0) {
        return NULL;
    }
    if (PyObject_GetBuffer(out, &out_view, flags | PyBUF_WRITABLE) != 0) {
        PyBuffer_Release(&in_view);
        return NULL;
    }

    if (in_view.ndim != 1) {
        PyErr_SetString(PyExc_ValueError, "First argument should be one-dimensional.");
        goto exit;
    }
    if (out_view.ndim != 1) {
        PyErr_SetString(PyExc_ValueError, "Second argument should be one-dimensional.");
        goto exit;
    }

    char format;
    char out_format;
    int status = 0;

    if (in_view.format == NULL || out_view.format == NULL) {
        PyErr_SetString(PyExc_ValueError, "Buffer format is not available.");
        goto exit;
    }

    format = native_type_code(in_view.format);
    if (format == '\0') {
        PyErr_Format(PyExc_ValueError,
            "Unexpected buffer format: %s (expected one type code in native "
            "byte order)",
            in_view.format);
        goto exit;
    }
    out_format = native_type_code(out_view.format);
    if (out_format == '\0') {
        PyErr_Format(PyExc_ValueError,
            "Unexpected output buffer format: %s (expected one type code in "
            "native byte order)",
            out_view.format);
        goto exit;
    }

    if (format == 'B' && in_view.itemsize == sizeof(uint8_t)) {
        if ((out_format != 'I' && out_format != 'L') ||
            out_view.itemsize != sizeof(uint32_t)) {
            PyErr_SetString(PyExc_ValueError,
                "Output buffer should contain 32-bit unsigned integers.");
            goto exit;
        }
        status = integer_unpack_u8(&in_view, &out_view);
    }
    else if (format == 'H' && in_view.itemsize == sizeof(uint16_t)) {
        if ((out_format != 'I' && out_format != 'L') ||
            out_view.itemsize != sizeof(uint32_t)) {
            PyErr_SetString(PyExc_ValueError,
                "Output buffer should contain 32-bit unsigned integers.");
            goto exit;
        }
        status = integer_unpack_u16(&in_view, &out_view);
    }
    else if (format == 'b' && in_view.itemsize == sizeof(int8_t)) {
        if ((out_format != 'i' && out_format != 'l') ||
            out_view.itemsize != sizeof(int32_t)) {
            PyErr_SetString(PyExc_ValueError,
                "Output buffer should contain 32-bit signed integers.");
            goto exit;
        }
        status = integer_unpack_i8(&in_view, &out_view);
    }
    else if (format == 'h' && in_view.itemsize == sizeof(int16_t)) {
        if ((out_format != 'i' && out_format != 'l') ||
            out_view.itemsize != sizeof(int32_t)) {
            PyErr_SetString(PyExc_ValueError,
                "Output buffer should contain 32-bit signed integers.");
            goto exit;
        }
        status = integer_unpack_i16(&in_view, &out_view);
    }
    else {
        PyErr_Format(PyExc_ValueError,
            "Unexpected buffer format: %s",
            in_view.format);
        goto exit;
    }

    if (status != 0) {
        if (status == -1) {
            PyErr_SetString(PyExc_ValueError, "Output buffer is too small.");
        }
        else if (status == -2) {
            PyErr_SetString(PyExc_ValueError, "Packed integer is truncated.");
        }
        else {
            PyErr_SetString(PyExc_ValueError, "Output buffer is too large.");
        }
    }

exit:
    PyBuffer_Release(&in_view);
    PyBuffer_Release(&out_view);
    if (PyErr_Occurred()) {
        return NULL;
    }
    Py_RETURN_NONE;
}

static PyMethodDef IntegerUnpackMethods[] = {
    {"integer_unpack", integer_unpack, METH_VARARGS, NULL},
    {NULL, NULL, 0, NULL}
};

static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    "_bcif_helper",
    NULL,
    -1,
    IntegerUnpackMethods
};

PyMODINIT_FUNC
PyInit__bcif_helper(void)
{
    PyObject *m;

    m = PyModule_Create(&moduledef);
    if (!m) {
        return NULL;
    }

    return m;
}
