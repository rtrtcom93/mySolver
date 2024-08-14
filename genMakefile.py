#!/usr/bin/env python

import os
import argparse

compiler_type = {
    'gnu': {'cpp':'g++', 'c':'gcc', 'fortran':'gfortran'},
    'clang': {'cpp':'clang++', 'c':'clang', 'fortran':'flang'},
    'intel': {'cpp':'icpc', 'c':'icc', 'fortran':'ifort'},
}

macro_variable = { 'cpp':'CXX', 'c':'CC', 'fortran':'FC' }
default_lang_spec = {'cpp':'c++11', 'c':'c99', 'fortran':'f95'}
default_lang_ext = {'cpp':'cpp', 'c':'c', 'fortran':'f90'}

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-l', '--lang', metavar='<language>', dest='lang', default='cpp', choices=['cpp', 'c', 'fortran'],
        help='programming language. available list is cpp, c, fortran. default is cpp')
    parser.add_argument('-c', '--compiler-type', metavar='<compiler-type>', dest='type', default='gnu',
        choices=['gnu', 'clang', 'intel'], help='compiler type. available list is gnu, clang, intel. default is gnu')
    parser.add_argument('-s', '--spec', metavar='<language spec>', dest='spec',
        help='language spec. default is c++11 for cpp, c99 for c and f90 for fortran')
    group = parser.add_mutually_exclusive_group()
    group.add_argument('--shared', action='store_true', dest='shared', help='shared library')
    group.add_argument('--static', action='store_true', dest='static', help='static library')
    parser.add_argument('-d', '--debug', action='store_true', help='turn on debug flag')
    parser.add_argument('-e', '--ext', metavar='<file extension>', dest='ext',
        help='file extension for source codes. default is cpp for cpp, c for c and f95 for fortran')
    parser.add_argument('-o', '--output', metavar='<output>', dest='out', default='main',
        help='output file name without extension')

    args = parser.parse_args()

if args.ext is not None:
    ext = args.ext
else:
    ext = default_lang_ext[args.lang]
    
sources = [f for f in os.listdir('.') if os.path.isfile(f) and f.endswith(ext)]
lines = []
lines.append('{0} = {1}'.format(macro_variable[args.lang], compiler_type[args.type][args.lang]))
if args.spec is not None:
    if args.type == 'intel' and args.lang == 'fortran':
        lines.append('{0}FLAGS = {1}'.format(macro_variable[args.lang], '-std{0}'.format(args.spec[1:])))
    else:
        lines.append('{0}FLAGS = {1}'.format(macro_variable[args.lang], '-std={0}'.format(args.spec)))
else:
    if args.type == 'intel' and args.lang == 'fortran':
        lines.append('{0}FLAGS = {1}'.format(macro_variable[args.lang],
            '-std{0}'.format(default_lang_spec[args.lang][1:])))
    else:
        lines.append('{0}FLAGS = {1}'.format(macro_variable[args.lang],
            '-std={0}'.format(default_lang_spec[args.lang])))

if args.shared:
    lines[-1] += (' -shared -fPIC')

if args.debug:
    lines[-1] += ' -g'
 
lines.append('LDFLAGS =')
lines.append('INCLUDE = -I.')
lines.append('LIB = -L.')
lines.append('SRC =')
for source in sources:
    lines[-1] += ' ' + source
lines.append('OBJ = $(SRC:%.{0}=%.o)'.format(ext))
prefix = ''
postfix = 'exe'
if args.shared:
    prefix = 'lib'
    postfix = 'so'
elif args.static:
    prefix = 'lib'
    postfix = 'a'

lines.append('OUT = {0}{1}.{2}\n'.format(prefix, args.out, postfix))
lines.append('%.o:%.{0}'.format(ext))
lines.append('\t$({0}) $({0}FLAGS) $(INCLUDE) -c $<\n'.format(macro_variable[args.lang]))
lines.append('all:build\n')
lines.append('build:$(OBJ)')
if not args.static:
    lines.append('\t$({0}) $({0}FLAGS) $(INCLUDE) -o $(OUT) $^ $(LDFLAGS) $(LIB)\n'.format(macro_variable[args.lang]))
else:
    lines.append('\tar rcs $(OUT) $^\n')

lines.append('clean:')
lines.append('\trm -f $(OBJ) $(OUT)\n')
lines.append('.PHONY:all build clean')

with open('Makefile', 'w') as fp:
    for line in lines:
        fp.write(line+'\n')