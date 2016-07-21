env = Environment(
    CPPPATH = ['src'],
    CCFLAGS='-std=c++11 -g -pg',
    LIBS=['pll','pthread'],
    LINKFLAGS=['-g', '-pg'])

# Doubles compilation time.
#env.Append(CCFLAGS='-O3 -msse2')

env.Program(
    target='test',
    source=[Glob('src/*.cpp')])
