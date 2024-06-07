from setuptools import setup, Extension
import os
import platform

openmm_dir = '@OPENMM_DIR@'
plugin_header_dir = '@RMSDPLUSFORCEPLUGIN_HEADER_DIR@'
plugin_library_dir = '@RMSDPLUSFORCEPLUGIN_LIBRARY_DIR@'

# setup extra compile and link arguments on Mac
extra_compile_args = ['-std=c++11']
extra_link_args = []

if platform.system() == 'Darwin':
    extra_compile_args += ['-stdlib=libc++', '-mmacosx-version-min=10.7']
    extra_link_args += ['-stdlib=libc++', '-mmacosx-version-min=10.7', '-Wl', '-rpath', openmm_dir+'/lib']

extension = Extension(name='_rmsdplusforceplugin',
                      sources=['RMSDPlusForcePluginWrapper.cpp'],
                      libraries=['OpenMM', 'RMSDPlusForcePlugin'],
                      include_dirs=[os.path.join(openmm_dir, 'include'), plugin_header_dir],
                      library_dirs=[os.path.join(openmm_dir, 'lib'), plugin_library_dir],
                      extra_compile_args=extra_compile_args,
                      extra_link_args=extra_link_args
                     )

setup(name='rmsdplusforceplugin',
      version='1.0',
      py_modules=['rmsdplusforceplugin'],
      ext_modules=[extension],
     )