/* Derived from stringobject.c in CPython 2.7

Only supports bytes for Python 3.  Don't compile with Python 2.7

This can hopefully go away when % is supported for PyBytes in Python 3.5
in September 2015 
*/

#define PY_SSIZE_T_CLEAN

#include "Python.h"
#include <ctype.h>
#include <stddef.h>

Py_LOCAL_INLINE(PyObject *)
getnextarg(PyObject *args, Py_ssize_t arglen, Py_ssize_t *p_argidx)
{
    Py_ssize_t argidx = *p_argidx;
    if (argidx < arglen) {
        (*p_argidx)++;
        if (arglen < 0)
            return args;
        else
            return PyTuple_GetItem(args, argidx);
    }
    PyErr_SetString(PyExc_TypeError,
                    "not enough arguments for format string");
    return NULL;
}

/* Format codes
 * F_LJUST      '-'
 * F_SIGN       '+'
 * F_BLANK      ' '
 * F_ALT        '#'
 * F_ZERO       '0'
 */
#define F_LJUST (1<<0)
#define F_SIGN  (1<<1)
#define F_BLANK (1<<2)
#define F_ALT   (1<<3)
#define F_ZERO  (1<<4)

/* Returns a new reference to a PyString object, or NULL on failure. */

static PyObject *
formatfloat(PyObject *v, int flags, int prec, int type)
{
    char *p;
    PyObject *result;
    double x;

    x = PyFloat_AsDouble(v);
    if (x == -1.0 && PyErr_Occurred()) {
        PyErr_Format(PyExc_TypeError, "float argument required, "
                     "not %.200s", Py_TYPE(v)->tp_name);
        return NULL;
    }

    if (prec < 0)
        prec = 6;

    p = PyOS_double_to_string(x, type, prec,
                              (flags & F_ALT) ? Py_DTSF_ALT : 0, NULL);

    if (p == NULL)
        return NULL;
    // result = PyString_FromStringAndSize(p, strlen(p));
    result = PyBytes_FromStringAndSize(p, strlen(p));
    PyMem_Free(p);
    return result;
}

/* _PyString_FormatLong emulates the format codes d, u, o, x and X, and
 * the F_ALT flag, for Python's long (unbounded) ints.  It's not used for
 * Python's regular ints.
 * Return value:  a new PyString*, or NULL if error.
 *  .  *pbuf is set to point into it,
 *     *plen set to the # of chars following that.
 *     Caller must decref it when done using pbuf.
 *     The string starting at *pbuf is of the form
 *         "-"? ("0x" | "0X")? digit+
 *     "0x"/"0X" are present only for x and X conversions, with F_ALT
 *         set in flags.  The case of hex digits will be correct,
 *     There will be at least prec digits, zero-filled on the left if
 *         necessary to get that many.
 * val          object to be converted
 * flags        bitmask of format flags; only F_ALT is looked at
 * prec         minimum number of digits; 0-fill on left if needed
 * type         a character in [duoxX]; u acts the same as d
 *
 * CAUTION:  o, x and X conversions on regular ints can never
 * produce a '-' sign, but can for Python's unbounded ints.
 */
PyObject*
_PyString_FormatLong(PyObject *val, int flags, int prec, int type,
                     char **pbuf, int *plen)
{
    PyObject *result = NULL;
    char *buf;
    Py_ssize_t i;
    int sign;           /* 1 if '-', else 0 */
    int len;            /* number of characters */
    Py_ssize_t llen;
    int numdigits;      /* len == numnondigits + numdigits */
    int numnondigits = 0;

    switch (type) {
    case 'd':
    case 'u':
        // result = Py_TYPE(val)->tp_str(val);
        result = PyNumber_ToBase(val, 10);
        break;
    case 'o':
        // result = Py_TYPE(val)->tp_as_number->nb_oct(val);
    result = PyNumber_ToBase(val, 8);
        break;
    case 'x':
    case 'X':
        numnondigits = 2;
        // result = Py_TYPE(val)->tp_as_number->nb_hex(val);
        result = PyNumber_ToBase(val, 16);
        break;
    default:
        assert(!"'type' not in [duoxX]");
    }
    if (!result)
        return NULL;

    // buf = PyString_AsString(result);
    buf = PyBytes_AsString(result);
    if (!buf) {
        Py_DECREF(result);
        return NULL;
    }

    /* To modify the string in-place, there can only be one reference. */
    if (Py_REFCNT(result) != 1) {
        PyErr_BadInternalCall();
        return NULL;
    }
    // llen = PyString_Size(result);
    llen = PyBytes_Size(result);
    if (llen > INT_MAX) {
        PyErr_SetString(PyExc_ValueError, "string too large in _PyString_FormatLong");
        return NULL;
    }
    len = (int)llen;
    if (buf[len-1] == 'L') {
        --len;
        buf[len] = '\0';
    }
    sign = buf[0] == '-';
    numnondigits += sign;
    numdigits = len - numnondigits;
    assert(numdigits > 0);

    /* Get rid of base marker unless F_ALT */
    if ((flags & F_ALT) == 0) {
        /* Need to skip 0x, 0X or 0. */
        int skipped = 0;
        switch (type) {
        case 'o':
            assert(buf[sign] == '0');
            /* If 0 is only digit, leave it alone. */
            if (numdigits > 1) {
                skipped = 1;
                --numdigits;
            }
            break;
        case 'x':
        case 'X':
            assert(buf[sign] == '0');
            assert(buf[sign + 1] == 'x');
            skipped = 2;
            numnondigits -= 2;
            break;
        }
        if (skipped) {
            buf += skipped;
            len -= skipped;
            if (sign)
                buf[0] = '-';
        }
        assert(len == numnondigits + numdigits);
        assert(numdigits > 0);
    }

    /* Fill with leading zeroes to meet minimum width. */
    if (prec > numdigits) {
        // PyObject *r1 = PyString_FromStringAndSize(NULL,
        //                         numnondigits + prec);
        PyObject *r1 = PyBytes_FromStringAndSize(NULL,
                                                    numnondigits + prec);
        char *b1;
        if (!r1) {
            Py_DECREF(result);
            return NULL;
        }
        // b1 = PyString_AS_STRING(r1);
        b1 = PyBytes_AS_STRING(r1);
        for (i = 0; i < numnondigits; ++i)
            *b1++ = *buf++;
        for (i = 0; i < prec - numdigits; i++)
            *b1++ = '0';
        for (i = 0; i < numdigits; i++)
            *b1++ = *buf++;
        *b1 = '\0';
        Py_DECREF(result);
        result = r1;
        // buf = PyString_AS_STRING(result);
        buf = PyBytes_AS_STRING(result);
        len = numnondigits + prec;
    }

    /* Fix up case for hex conversions. */
    if (type == 'X') {
        /* Need to convert all lower case letters to upper case.
           and need to convert 0x to 0X (and -0x to -0X). */
        for (i = 0; i < len; i++)
            if (buf[i] >= 'a' && buf[i] <= 'x')
                buf[i] -= 'a'-'A';
    }
    *pbuf = buf;
    *plen = len;
    return result;
}

Py_LOCAL_INLINE(int)
formatint(char *buf, size_t buflen, int flags,
          int prec, int type, PyObject *v)
{
    /* fmt = '%#.' + `prec` + 'l' + `type`
       worst case length = 3 + 19 (worst len of INT_MAX on 64-bit machine)
       + 1 + 1 = 24 */
    char fmt[64];       /* plenty big enough! */
    char *sign;
    long x;

    // x = PyInt_AsLong(v);
    x = PyLong_AsLong(v);
    if (x == -1 && PyErr_Occurred()) {
        PyErr_Format(PyExc_TypeError, "int argument required, not %.200s",
                     Py_TYPE(v)->tp_name);
        return -1;
    }
    if (x < 0 && type == 'u') {
        type = 'd';
    }
    if (x < 0 && (type == 'x' || type == 'X' || type == 'o'))
        sign = "-";
    else
        sign = "";
    if (prec < 0)
        prec = 1;

    if ((flags & F_ALT) &&
        (type == 'x' || type == 'X')) {
        /* When converting under %#x or %#X, there are a number
         * of issues that cause pain:
         * - when 0 is being converted, the C standard leaves off
         *   the '0x' or '0X', which is inconsistent with other
         *   %#x/%#X conversions and inconsistent with Python's
         *   hex() function
         * - there are platforms that violate the standard and
         *   convert 0 with the '0x' or '0X'
         *   (Metrowerks, Compaq Tru64)
         * - there are platforms that give '0x' when converting
         *   under %#X, but convert 0 in accordance with the
         *   standard (OS/2 EMX)
         *
         * We can achieve the desired consistency by inserting our
         * own '0x' or '0X' prefix, and substituting %x/%X in place
         * of %#x/%#X.
         *
         * Note that this is the same approach as used in
         * formatint() in unicodeobject.c
         */
        PyOS_snprintf(fmt, sizeof(fmt), "%s0%c%%.%dl%c",
                      sign, type, prec, type);
    }
    else {
        PyOS_snprintf(fmt, sizeof(fmt), "%s%%%s.%dl%c",
                      sign, (flags&F_ALT) ? "#" : "",
                      prec, type);
    }

    /* buf = '+'/'-'/'' + '0'/'0x'/'' + '[0-9]'*max(prec, len(x in octal))
     * worst case buf = '-0x' + [0-9]*prec, where prec >= 11
     */
    if (buflen <= 14 || buflen <= (size_t)3 + (size_t)prec) {
        PyErr_SetString(PyExc_OverflowError,
            "formatted integer is too long (precision too large?)");
        return -1;
    }
    if (sign[0])
        PyOS_snprintf(buf, buflen, fmt, -x);
    else
        PyOS_snprintf(buf, buflen, fmt, x);
    return (int)strlen(buf);
}

Py_LOCAL_INLINE(int)
formatchar(char *buf, size_t buflen, PyObject *v)
{
    /* presume that the buffer is at least 2 characters long */
    // if (PyString_Check(v)) {
    if (PyBytes_Check(v)) {
        if (!PyArg_Parse(v, "c;%c requires int or char", &buf[0]))
            return -1;
    }
    else {
        if (!PyArg_Parse(v, "b;%c requires int or char", &buf[0]))
            return -1;
    }
    buf[1] = '\0';
    return 1;
}

/* fmt%(v1,v2,...) is roughly equivalent to sprintf(fmt, v1, v2, ...)
   FORMATBUFLEN is the length of the buffer in which the ints &
   chars are formatted. XXX This is a magic number. Each formatting
   routine does bounds checking to ensure no overflow, but a better
   solution may be to malloc a buffer of appropriate size for each
   format. For now, the current solution is sufficient.
*/
#define FORMATBUFLEN (size_t)120    

PyObject *
// PyString_Format(PyObject *format, PyObject *args)
PyBytes_Format(PyObject *format, PyObject *args)
{
    char *fmt, *res;
    Py_ssize_t arglen, argidx;
    Py_ssize_t reslen, rescnt, fmtcnt;
    int args_owned = 0;
    PyObject *result, *orig_args;
// #ifdef Py_USING_UNICODE
//     PyObject *v, *w;
// #endif
    PyObject *dict = NULL;
    // if (format == NULL || !PyString_Check(format) || args == NULL) {
    if (format == NULL || !PyBytes_Check(format) || args == NULL) {
        PyErr_BadInternalCall();
        return NULL;
    }
    orig_args = args;
    // fmt = PyString_AS_STRING(format);
    fmt = PyBytes_AS_STRING(format);
    // fmtcnt = PyString_GET_SIZE(format);
    fmtcnt = PyBytes_GET_SIZE(format);
    reslen = rescnt = fmtcnt + 100;
    // result = PyString_FromStringAndSize((char *)NULL, reslen);
    result = PyBytes_FromStringAndSize((char *)NULL, reslen);
    if (result == NULL)
        return NULL;
    // res = PyString_AsString(result);
    res = PyBytes_AsString(result);
    if (PyTuple_Check(args)) {
        arglen = PyTuple_GET_SIZE(args);
        argidx = 0;
    }
    else {
        arglen = -1;
        argidx = -2;
    }
    // if (Py_TYPE(args)->tp_as_mapping && Py_TYPE(args)->tp_as_mapping->mp_subscript &&
    //     !PyTuple_Check(args) && !PyObject_TypeCheck(args, &PyBaseString_Type))
    if (Py_TYPE(args)->tp_as_mapping && Py_TYPE(args)->tp_as_mapping->mp_subscript &&
            !PyTuple_Check(args) && !PyObject_TypeCheck(args, &PyBytes_Type))
        dict = args;
    while (--fmtcnt >= 0) {
        if (*fmt != '%') {
            if (--rescnt < 0) {
                rescnt = fmtcnt + 100;
                reslen += rescnt;
                // if (_PyString_Resize(&result, reslen))
                if (_PyBytes_Resize(&result, reslen)) {
                    return NULL;
                }
                // res = PyString_AS_STRING(result)
                res = PyBytes_AS_STRING(result)
                            + reslen - rescnt;
                --rescnt;
            }
            *res++ = *fmt++;
        }
        else {
            /* Got a format specifier */
            int flags = 0;
            Py_ssize_t width = -1;
            int prec = -1;
            int c = '\0';
            int fill;
            int isnumok;
            PyObject *v = NULL;
            PyObject *temp = NULL;
            char *pbuf;
            int sign;
            Py_ssize_t len;
            char formatbuf[FORMATBUFLEN];
                 /* For format{int,char}() */
// #ifdef Py_USING_UNICODE
//             char *fmt_start = fmt;
//             Py_ssize_t argidx_start = argidx;
// #endif

            fmt++;
            if (*fmt == '(') {
                char *keystart;
                Py_ssize_t keylen;
                PyObject *key;
                int pcount = 1;

                if (dict == NULL) {
                    PyErr_SetString(PyExc_TypeError,
                             "format requires a mapping");
                    goto error;
                }
                ++fmt;
                --fmtcnt;
                keystart = fmt;
                /* Skip over balanced parentheses */
                while (pcount > 0 && --fmtcnt >= 0) {
                    if (*fmt == ')')
                        --pcount;
                    else if (*fmt == '(')
                        ++pcount;
                    fmt++;
                }
                keylen = fmt - keystart - 1;
                if (fmtcnt < 0 || pcount > 0) {
                    PyErr_SetString(PyExc_ValueError,
                               "incomplete format key");
                    goto error;
                }
                // key = PyString_FromStringAndSize(keystart,
                //                                  keylen);
                key = PyBytes_FromStringAndSize(keystart,
                                 keylen);
                if (key == NULL)
                    goto error;
                if (args_owned) {
                    Py_DECREF(args);
                    args_owned = 0;
                }
                args = PyObject_GetItem(dict, key);
                Py_DECREF(key);
                if (args == NULL) {
                    goto error;
                }
                args_owned = 1;
                arglen = -1;
                argidx = -2;
            }
            while (--fmtcnt >= 0) {
                switch (c = *fmt++) {
                case '-': flags |= F_LJUST; continue;
                case '+': flags |= F_SIGN; continue;
                case ' ': flags |= F_BLANK; continue;
                case '#': flags |= F_ALT; continue;
                case '0': flags |= F_ZERO; continue;
                }
                break;
            }
            if (c == '*') {
                v = getnextarg(args, arglen, &argidx);
                if (v == NULL)
                    goto error;
                // if (!PyInt_Check(v)) {
                if (!PyLong_Check(v)) {
                    PyErr_SetString(PyExc_TypeError,
                                    "* wants int");
                    goto error;
                }
                // width = PyInt_AsSsize_t(v);
                width = PyLong_AsSsize_t(v);
                if (width == -1 && PyErr_Occurred())
                    goto error;
                if (width < 0) {
                    flags |= F_LJUST;
                    width = -width;
                }
                if (--fmtcnt >= 0)
                    c = *fmt++;
            }
            else if (c >= 0 && isdigit(c)) {
                width = c - '0';
                while (--fmtcnt >= 0) {
                    c = Py_CHARMASK(*fmt++);
                    if (!isdigit(c))
                        break;
                    if (width > (PY_SSIZE_T_MAX - ((int)c - '0')) / 10) {
                        PyErr_SetString(
                            PyExc_ValueError,
                            "width too big");
                        goto error;
                    }
                    width = width*10 + (c - '0');
                }
            }
            if (c == '.') {
                prec = 0;
                if (--fmtcnt >= 0)
                    c = *fmt++;
                if (c == '*') {
                    v = getnextarg(args, arglen, &argidx);
                    if (v == NULL)
                        goto error;
                    // if (!PyInt_Check(v)) {
                    if (!PyLong_Check(v)) {
                        PyErr_SetString(
                            PyExc_TypeError,
                            "* wants int");
                        goto error;
                    }
                    // prec = _PyInt_AsInt(v);
                    prec = _PyLong_AsInt(v);
                    if (prec == -1 && PyErr_Occurred())
                        goto error;

                    if (prec < 0)
                        prec = 0;
                    if (--fmtcnt >= 0)
                        c = *fmt++;
                }
                else if (c >= 0 && isdigit(c)) {
                    prec = c - '0';
                    while (--fmtcnt >= 0) {
                        c = Py_CHARMASK(*fmt++);
                        if (!isdigit(c))
                            break;
                        if (prec > (INT_MAX - ((int)c - '0')) / 10) {
                            PyErr_SetString(
                                PyExc_ValueError,
                                "prec too big");
                            goto error;
                        }
                        prec = prec*10 + (c - '0');
                    }
                }
            } /* prec */
            if (fmtcnt >= 0) {
                if (c == 'h' || c == 'l' || c == 'L') {
                    if (--fmtcnt >= 0)
                        c = *fmt++;
                }
            }
            if (fmtcnt < 0) {
                PyErr_SetString(PyExc_ValueError,
                                "incomplete format");
                goto error;
            }
            if (c != '%') {
                v = getnextarg(args, arglen, &argidx);
                if (v == NULL)
                    goto error;
            }
            sign = 0;
            fill = ' ';
            switch (c) {
            case '%':
                pbuf = "%";
                len = 1;
                break;
            case 's':
// #ifdef Py_USING_UNICODE
//                 if (PyUnicode_Check(v)) {
//                     fmt = fmt_start;
//                     argidx = argidx_start;
//                     goto unicode;
//                 }
// #endif
                // temp = _PyObject_Str(v);
                temp = PyObject_Bytes(v);
// #ifdef Py_USING_UNICODE
//                 if (temp != NULL && PyUnicode_Check(temp)) {
//                     Py_DECREF(temp);
//                     fmt = fmt_start;
//                     argidx = argidx_start;
//                     goto unicode;
//                 }
// #endif
                /* Fall through */
            case 'r':
                if (c == 'r')
                    temp = PyObject_Repr(v);
                if (temp == NULL)
                    goto error;
                // if (!PyString_Check(temp)) {
                if (!PyBytes_Check(temp)) {
                    PyErr_SetString(PyExc_TypeError,
                      "%s argument has non-string str()");
                    Py_DECREF(temp);
                    goto error;
                }
                // pbuf = PyString_AS_STRING(temp);
                pbuf = PyBytes_AS_STRING(temp);
                // len = PyString_GET_SIZE(temp);
                len = PyBytes_GET_SIZE(temp);
                if (prec >= 0 && len > prec)
                    len = prec;
                break;
            case 'i':
            case 'd':
            case 'u':
            case 'o':
            case 'x':
            case 'X':
                if (c == 'i')
                    c = 'd';
                isnumok = 0;
                if (PyNumber_Check(v)) {
                    PyObject *iobj = NULL;

                    // if (PyInt_Check(v) || (PyLong_Check(v))) {
                    if (PyLong_Check(v)) {
                        iobj = v;
                        Py_INCREF(iobj);
                    }
                    else {
                        // iobj = PyNumber_Int(v);
                        // if (iobj == NULL) {
                        //     PyErr_Clear();
                        //     iobj = PyNumber_Long(v);
                        // }
                        iobj = PyNumber_Long(v);
                    }
                    if (iobj != NULL) {
                        // if (PyInt_Check(iobj)) {
                        if (PyLong_Check(iobj)) {
                            isnumok = 1;
                            pbuf = formatbuf;
                            len = formatint(pbuf,
                                            sizeof(formatbuf),
                                            flags, prec, c, iobj);
                            Py_DECREF(iobj);
                            if (len < 0)
                                goto error;
                            sign = 1;
                        }
                        else if (PyLong_Check(iobj)) {
                            int ilen;

                            isnumok = 1;
                            temp = _PyString_FormatLong(iobj, flags,
                                prec, c, &pbuf, &ilen);
                            Py_DECREF(iobj);
                            len = ilen;
                            if (!temp)
                                goto error;
                            sign = 1;
                        }
                        else {
                            Py_DECREF(iobj);
                        }
                    }
                }
                if (!isnumok) {
                    PyErr_Format(PyExc_TypeError,
                        "%%%c format: a number is required, "
                        "not %.200s", c, Py_TYPE(v)->tp_name);
                    goto error;
                }
                if (flags & F_ZERO)
                    fill = '0';
                break;
            case 'e':
            case 'E':
            case 'f':
            case 'F':
            case 'g':
            case 'G':
                temp = formatfloat(v, flags, prec, c);
                if (temp == NULL)
                    goto error;
                // pbuf = PyString_AS_STRING(temp);
                pbuf = PyBytes_AS_STRING(temp);
                // len = PyString_GET_SIZE(temp);
                len = PyBytes_GET_SIZE(temp);
                sign = 1;
                if (flags & F_ZERO)
                    fill = '0';
                break;
            case 'c':
// #ifdef Py_USING_UNICODE
//                 if (PyUnicode_Check(v)) {
//                     fmt = fmt_start;
//                     argidx = argidx_start;
//                     goto unicode;
//                 }
// #endif
                pbuf = formatbuf;
                len = formatchar(pbuf, sizeof(formatbuf), v);
                if (len < 0)
                    goto error;
                break;
            default:
                PyErr_Format(PyExc_ValueError,
                  "unsupported format character '%c' (0x%x) "
                  "at index %zd",
                  c, c,
                  // (Py_ssize_t)(fmt - 1 -
                  //              PyString_AsString(format)));
                  (Py_ssize_t)(fmt - 1 -
                               PyBytes_AsString(format)));
                goto error;
            }
            if (sign) {
                if (*pbuf == '-' || *pbuf == '+') {
                    sign = *pbuf++;
                    len--;
                }
                else if (flags & F_SIGN)
                    sign = '+';
                else if (flags & F_BLANK)
                    sign = ' ';
                else
                    sign = 0;
            }
            if (width < len)
                width = len;
            if (rescnt - (sign != 0) < width) {
                reslen -= rescnt;
                rescnt = width + fmtcnt + 100;
                reslen += rescnt;
                if (reslen < 0) {
                    Py_DECREF(result);
                    Py_XDECREF(temp);
                    return PyErr_NoMemory();
                }
                // if (_PyString_Resize(&result, reslen)) {
                if (_PyBytes_Resize(&result, reslen)) {
                    Py_XDECREF(temp);
                    return NULL;
                }
                // res = PyString_AS_STRING(result)
                //     + reslen - rescnt;
                res = PyBytes_AS_STRING(result)
                    + reslen - rescnt;
            }
            if (sign) {
                if (fill != ' ')
                    *res++ = sign;
                rescnt--;
                if (width > len)
                    width--;
            }
            if ((flags & F_ALT) && (c == 'x' || c == 'X')) {
                assert(pbuf[0] == '0');
                assert(pbuf[1] == c);
                if (fill != ' ') {
                    *res++ = *pbuf++;
                    *res++ = *pbuf++;
                }
                rescnt -= 2;
                width -= 2;
                if (width < 0)
                    width = 0;
                len -= 2;
            }
            if (width > len && !(flags & F_LJUST)) {
                do {
                    --rescnt;
                    *res++ = fill;
                } while (--width > len);
            }
            if (fill == ' ') {
                if (sign)
                    *res++ = sign;
                if ((flags & F_ALT) &&
                    (c == 'x' || c == 'X')) {
                    assert(pbuf[0] == '0');
                    assert(pbuf[1] == c);
                    *res++ = *pbuf++;
                    *res++ = *pbuf++;
                }
            }
            Py_MEMCPY(res, pbuf, len);
            res += len;
            rescnt -= len;
            while (--width >= len) {
                --rescnt;
                *res++ = ' ';
            }
            if (dict && (argidx < arglen) && c != '%') {
                PyErr_SetString(PyExc_TypeError,
                           "not all arguments converted during string formatting");
                Py_XDECREF(temp);
                goto error;
            }
            Py_XDECREF(temp);
        } /* '%' */
    } /* until end */
    if (argidx < arglen && !dict) {
        PyErr_SetString(PyExc_TypeError,
                        "not all arguments converted during string formatting");
        goto error;
    }
    if (args_owned) {
        Py_DECREF(args);
    }
    // if (_PyString_Resize(&result, reslen - rescnt))
    if (_PyBytes_Resize(&result, reslen - rescnt))
        return NULL;
    return result;

// #ifdef Py_USING_UNICODE
//  unicode:
//     if (args_owned) {
//         Py_DECREF(args);
//         args_owned = 0;
//     }
//     /* Fiddle args right (remove the first argidx arguments) */
//     if (PyTuple_Check(orig_args) && argidx > 0) {
//         PyObject *v;
//         Py_ssize_t n = PyTuple_GET_SIZE(orig_args) - argidx;
//         v = PyTuple_New(n);
//         if (v == NULL)
//             goto error;
//         while (--n >= 0) {
//             PyObject *w = PyTuple_GET_ITEM(orig_args, n + argidx);
//             Py_INCREF(w);
//             PyTuple_SET_ITEM(v, n, w);
//         }
//         args = v;
//     } else {
//         Py_INCREF(orig_args);
//         args = orig_args;
//     }
//     args_owned = 1;
//      /*Take what we have of the result and let the Unicode formatting
//        function format the rest of the input. */
//     rescnt = res - PyString_AS_STRING(result);
//     if (_PyString_Resize(&result, rescnt))
//         goto error;
//     fmtcnt = PyString_GET_SIZE(format) - \
//              (fmt - PyString_AS_STRING(format));
//     format = PyUnicode_Decode(fmt, fmtcnt, NULL, NULL);
//     if (format == NULL)
//         goto error;
//     v = PyUnicode_Format(format, args);
//     Py_DECREF(format);
//     if (v == NULL)
//         goto error;
//     /* Paste what we have (result) to what the Unicode formatting
//        function returned (v) and return the result (or error) */
//     w = PyUnicode_Concat(result, v);
//     Py_DECREF(result);
//     Py_DECREF(v);
//     Py_DECREF(args);
//     return w;
// #endif /* Py_USING_UNICODE */

 error:
    Py_DECREF(result);
    if (args_owned) {
        Py_DECREF(args);
    }
    return NULL;
}