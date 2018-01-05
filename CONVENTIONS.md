libnano styling and coding conventions
=====================================

#### Common argument names

    seq                             single DNA / RNA / amino acid sequence
    seq_length                      single DNA / RNA / amino acid sequence length
    seq1, seq2, ...                 multiple DNA / RNA / amino acid sequences
    seq1_length, seq2_length, ...   multiple DNA / RNA / amino acid sequence lengths
    temp_c                          temperature in degrees celsius
    temp_k                          temperature in degrees kelvin
    word_size                       sequence word size
    seq_int                         a seqint representation of a dna sequence

#### Common variable names

    i, j, k, l                      reserved as index variables for loops
    il, jl, kl, ll                  upper limit variables for loops

#### Common argument/variable prefixes / suffixes

    min_                            minimum value (e.g., min_temp_c)
    max_                            maximum value (e.g., max_temp_c)
    cur_                            current value, often used in loops
    _ptr                            pointer, often used in loops
    _arr                            c array or exposed c array
    _npy_arr                        numpy array (PyArrayObject)
    _frac                           fractional percent (e.g., .30 for 30%)
    _percent                        integer percent (e.g., 30 for 30%)
    
#### General coding conventions

    - We can get sequence lengths more or less for free from PyArg_ParseTuple
      so we will typically pass lengths around as integers to avoid having to
      call strlen or a similar function/macro to determine the length of a char
      array within a function

    - C functions that are not exposed to Python should be prefixed with two or
      more letters that unambigously reflect the file in which they are contained.
        e.g., all seqrepeat functions should be prefixed with "sr_"

        TODO: determine if C API-defined classes (e.g., RepeatCheck) should be
              part of this paradigm
