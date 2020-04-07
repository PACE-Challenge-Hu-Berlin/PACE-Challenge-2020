project('treedepth', 'cpp')

conf_data = configuration_data()
conf_data.set('TD_SANITY_CHECKS', get_option('td-sanity-checks'))

if get_option('td-static')
	add_project_link_arguments('-static-libgcc', '-static-libstdc++',
		language: 'cpp')
endif

configure_file(configuration: conf_data,
	input: 'config.hpp.in',
	output: 'config.hpp')

executable('td',
		'src/kernelization.cpp',
		'src/main.cpp',
		'src/mm.cpp',
		'src/naive-branching-solver.cpp',
		'src/simple-pid-solver.cpp',
		'src/utilities.cpp',
	cpp_args: '-std=c++17',
	install: true)
