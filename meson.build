project('treedepth', 'cpp')

if get_option('td-static')
	add_project_link_arguments('-static-libgcc', '-static-libstdc++',
		language: 'cpp')
endif

executable('td',
		'src/kernelization.cpp',
		'src/main.cpp',
		'src/mm.cpp',
		'src/naive-branching-solver.cpp',
		'src/simple-pid-solver.cpp',
		'src/utilities.cpp',
	cpp_args: '-std=c++17',
	install: true)
