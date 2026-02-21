project('pineapfel', 'cpp',
  version : '0.1.0',
  default_options : ['cpp_std=c++20', 'warning_level=2', 'optimization=3'],
)

# Dependencies
pineappl_dep = dependency('pineappl_capi', required : true)
yaml_cpp_dep = dependency('yaml-cpp', required : true)

# APFEL++ via apfelxx-config
apfelxx_config = find_program('apfelxx-config', required : true)

apfelxx_cflags = run_command(apfelxx_config, '--cflags', check : true).stdout().strip()
apfelxx_ldflags = run_command(apfelxx_config, '--ldflags', check : true).stdout().strip()

message('APFEL++ cflags: ' + apfelxx_cflags)
message('APFEL++ ldflags: ' + apfelxx_ldflags)

# Parse APFEL++ include dirs from cflags (-I flags)
apfelxx_inc_dirs = []
if apfelxx_cflags != ''
  foreach flag : apfelxx_cflags.split()
    if flag.startswith('-I')
      apfelxx_inc_dirs += include_directories(flag.substring(2))
    endif
  endforeach
endif

# Fallback: derive include dir from apfelxx-config location
if apfelxx_inc_dirs.length() == 0
  apfelxx_config_path = apfelxx_config.full_path()
  apfelxx_prefix = run_command('dirname', run_command('dirname', apfelxx_config_path, check : true).stdout().strip(), check : true).stdout().strip()
  apfelxx_inc_dirs += include_directories(apfelxx_prefix / 'include')
  message('APFEL++ include dirs (fallback): ' + apfelxx_prefix / 'include')
endif

# Parse APFEL++ link flags from ldflags
apfelxx_link_args = []
if apfelxx_ldflags != ''
  apfelxx_link_args = apfelxx_ldflags.split()
endif

apfelxx_dep = declare_dependency(
  include_directories : apfelxx_inc_dirs,
  link_args : apfelxx_link_args,
)

pineapfel_lib = static_library('pineapfel',
  'pineapfel/src/cards.cpp',
  'pineapfel/src/evolution.cpp',
  'pineapfel/src/fill.cpp',
  'pineapfel/src/grid_gen.cpp',
  'pineapfel/src/sidis_init.cpp',
  include_directories : include_directories('pineapfel/include'),
  dependencies : [pineappl_dep, yaml_cpp_dep, apfelxx_dep],
  install : true,
)

pineapfel_dep = declare_dependency(
  include_directories : include_directories('pineapfel/include'),
  link_with : pineapfel_lib,
  dependencies : [pineappl_dep, yaml_cpp_dep, apfelxx_dep],
)

executable('pineapfel-evolve',
  'pineapfel_cli/main.cpp',
  dependencies : pineapfel_dep,
  install : true,
)

executable('pineapfel-build',
  'pineapfel_cli/build.cpp',
  dependencies : pineapfel_dep,
  install : true,
)

executable('create-grid-example',
  'examples/create_grid.cpp',
  dependencies : pineapfel_dep,
  install : false,
)

install_headers(
  'pineapfel/include/cards.h',
  'pineapfel/include/evolution.h',
  'pineapfel/include/fill.h',
  'pineapfel/include/grid_gen.h',
  'pineapfel/include/pineapfel.h',
  'pineapfel/include/sidis_api.h',
  install_dir : get_option('includedir'),
)
