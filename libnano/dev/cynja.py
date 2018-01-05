'''
libnano.cynja

jinja2 templating and cython meet
'''

from Cython.Build import cythonize

import os.path
import io

CYNJA_EXT = '.pyxj'

from jinja2 import Environment, FileSystemLoader

jinja_env_py = {
    'block_start_string': '{%',
    'block_end_string': '%}',
    'variable_start_string': '{{',
    'variable_end_string': '}}'
}

jinja_env_c = {
    'block_start_string': '{%',
    'block_end_string': '%}',
    'variable_start_string': '<|',
    'variable_end_string': '|>'
}


def render(in_path, in_filename, out_path, out_filename, cynja_vars=None):
    if cynja_vars is None:
        cynja_vars = {}
    env = Environment(loader=FileSystemLoader(in_path), **jinja_env_py)
    template = env.get_template(in_filename)
    out_string = template.render(**cynja_vars)
    with io.open(os.path.join(out_path, out_filename), 'w') as fd:
        fd.write(out_string)


def cynjanize(extension_list, cynja_vars=None, **kwargs):
    """
    for use on cython only templates
    """
    for extension in extension_list:
        for i, source in enumerate(extension.sources):
            in_filename = os.path.basename(source)
            in_basename, ext = os.path.splitext(in_filename)
            if ext == CYNJA_EXT:
                path = os.path.dirname(source)
                rendered_file = in_basename + '.pyx'
                render(path, in_filename,
                        path, rendered_file,
                        cynja_vars=cynja_vars)
                extension.sources[i] = os.path.join(path, rendered_file)
    cythonize(extension_list, **kwargs)


def renderList(filelist, root_path=None):
    """ for use on mixtures of files or files whose destination differs from the
    source, or simply a filename change is desired

    provide a list of tuples of parameters
        (   source name (string): source filename (*.pyxj),
            source path (string): path to the source file,
            destination name (string or None): destination filename,
            destination path (string or None): path to the source file,
            cynja_vars (dict): dictionary of vars to be templated
        )
    """
    for src_name, src_path, dest_name, dest_path, cynja_vars in filelist:
        if root_path is not None:
            src_path = os.path.join(root_path, src_path)

        if dest_path is None:
            dest_path = src_path
        elif root_path is not None:
            dest_path = os.path.join(root_path, dest_path)

        in_basename, ext = os.path.splitext(src_name)
        if ext == CYNJA_EXT:
            if dest_name is None:
                rendered_file = in_basename + '.pyx'
            else:
                rendered_file = dest_name
            render(src_path, src_name,
                   dest_path, rendered_file,
                   cynja_vars=cynja_vars)
        else:
            raise ValueError("expected type .pyxj got".format(ext))
# end def

class CynjaModule(tuple):
    """ tuple subclass so we can get tuple unpacking with initialization
    logic.  use this if you want to locate templates and output modules in
    different locations
    """
    def __new__(cls,    source_name, source_path,
                        destination_name, destination_path=None,
                        cynja_vars=None):
        if destination_path is None:
            destination_path = source_path
        return super(CynjaModule, cls).__new__(cls, (   source_name,
                                                        source_path,
                                                        destination_name,
                                                        destination_path,
                                                        cynja_vars ) )
    # end def
    def __init__(self, source_name, source_path,
                        destination_name=None, destination_path=None,
                        cynja_vars=None):
        if cynja_vars is None or not isinstance(cynja_vars, dict):
            raise ValueError("cynja_vars must be a dictionary")
        self.src_name = source_name
        self.src_path = source_path
        self.dest_name = destination_name
        self.dest_path = source_path if destination_path is None else destination_path
        self.cynja_vars = cynja_vars
    # end def
# end class

if __name__ == '__main__':
    try:
        from setuptools import setup, Extension
    except ImportError:
        from distutils.core import setup, Extension

    test_list = [Extension('tests.test_cynja', sources=['tests/test_cynja.pyxj'])]
    cynjanize(test_list, cynja_vars={'key_t': 'int'})