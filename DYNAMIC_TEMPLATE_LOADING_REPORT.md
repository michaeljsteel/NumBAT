# NumBAT Dynamic Waveguide Template Loading

## Executive summary

NumBAT's current waveguide-template system works, but it uses file-based dynamic
imports in a way that is older than current Python conventions. The main
observable consequence is that dynamically loaded geometry classes can fail to
unpickle in Python 3.14 multiprocessing workers using the `forkserver` start
method.

The immediate compatibility fix is to serialize dynamically loaded geometry
instances by their implementation file and class name. That fix is now present
in `backend/usermesh.py`.

The more Pythonic long-term direction is to give built-in templates stable,
normal import paths, load each implementation module only once, and avoid
sending complete `Structure` objects through multiprocessing queues where
possible.

## Current implementation

Template definitions are listed in
`backend/msh/builtin_waveguides.json`. Each entry currently identifies:

- an implementation filename, such as `builtin_meshes.py`;
- a class name, such as `CircOnion1`;
- a user-facing shape name, such as `circ_onion1`.

The loader in `backend/meshing/templates.py` currently uses the implementation
filename and shape name to dynamically load modules:

```python
spec = importlib.util.spec_from_file_location(wgnm, pth_mod)
py_mod = spec.loader.load_module()
```

This has four important consequences:

1. `load_module()` is deprecated.
2. The same implementation file can be loaded multiple times under different
   module names.
3. Names such as `circ_onion1` are transient aliases rather than stable Python
   import paths.
4. Classes loaded this way are more difficult to use with multiprocessing,
   serialization, introspection, and packaging tools.

## Multiprocessing failure

The tutorial creates a structure, deep-copies it into each queued task, and
passes those tasks to worker processes. The structure contains a geometry
instance whose class was loaded from the transient module alias
`circ_onion1`.

Under Python 3.14 on Linux, the default multiprocessing behavior can use
`forkserver`. Queue payloads are pickled and unpickled in the worker. The
worker does not necessarily have the parent process's dynamically created
`circ_onion1` module, so unpickling fails with:

```text
ModuleNotFoundError: No module named 'circ_onion1'
```

## Compatibility fix

`UserGeometryBase` now provides a custom reduction method. When serialized, a
geometry instance records:

- the source implementation file;
- the class name;
- the instance dictionary.

The worker reloads the implementation file and reconstructs the instance from
that state. This avoids requiring the transient alias to be importable in the
worker and works with `spawn` and `forkserver` as well as `fork`.

This is a suitable compatibility layer for existing file-based user templates.
It also makes the existing tutorials work without requiring callers to force
the `fork` start method.

## Why the current approach is not fully modern Python

### Deprecated import API

`spec.loader.load_module()` is deprecated. The modern low-level import sequence
is:

```python
spec = importlib.util.spec_from_file_location(module_name, module_path)
module = importlib.util.module_from_spec(spec)
sys.modules[module_name] = module
spec.loader.exec_module(module)
```

Registering the module in `sys.modules` before executing it is important because
it gives the module normal import-system behavior and supports modules that
refer to themselves during initialization.

### Repeated loading of one file

Many built-in entries point to `builtin_meshes.py`, but the current loader uses
the shape name as the module name. Consequently, the same source file may be
executed repeatedly:

```text
rectangular
circular
circ_onion1
circ_onion2
...
```

This can produce multiple copies of the same class definitions and repeated
module-level initialization. A module cache keyed by resolved implementation
path would avoid this.

### Unstable module identity

The name `circ_onion1.CircOnion1` is not a normal package-qualified Python
symbol. A stable name such as
`msh.builtin_meshes.CircOnion1` is easier for Python to import, pickle, inspect,
and type-check.

## Recommended modernization

### Phase 1: modernize the existing file loader

Keep support for arbitrary user implementation files, but:

1. replace `load_module()` with `module_from_spec()` and `exec_module()`;
2. register modules in `sys.modules`;
3. cache modules by resolved file path;
4. use a stable generated module name based on the implementation path;
5. preserve the custom geometry reduction as a legacy compatibility layer.

This is relatively low risk and does not require changing the JSON format.

#### Implementation status

Implemented in `backend/meshing/template_loader.py`. The shared loader derives
a deterministic internal module name from each resolved implementation path,
registers the module in `sys.modules` before executing it, and caches only
successfully executed modules. It also removes a partially initialized module
from `sys.modules` if execution fails.

Both `backend/meshing/templates.py` and `backend/usermesh.py` use this shared
loader, so template registration and geometry reconstruction use the same
module identity and cache. The existing `UserGeometryBase.__reduce__`
compatibility mechanism remains in place for file-based user templates.

### Phase 2: use normal imports for built-in templates

Change built-in template metadata from filenames to import paths:

```json
{
  "module": "msh.builtin_meshes",
  "class": "CircOnion1",
  "inc_shape": "circ_onion1"
}
```

The loader can then use:

```python
from importlib import import_module

module = import_module(entry["module"])
template_class = getattr(module, entry["class"])
```

This gives built-in classes stable identities and makes them naturally
compatible with Python multiprocessing.

#### Implementation status

Implemented by replacing each built-in `wg_impl` entry with
`"module": "msh.builtin_meshes"` in `backend/msh/builtin_waveguides.json`.
The template loader now imports entries with `module` using
`importlib.import_module()` and retains `wg_impl` support for file-based user
templates. Built-in geometry classes consequently retain their natural
`msh.builtin_meshes` module identity.

### Phase 3: provide a package-based plugin mechanism

For third-party templates, the most maintainable approach is to use installed
Python packages and entry points rather than arbitrary source-file paths. A
plugin could declare an entry point such as:

```toml
[project.entry-points."numbat.waveguides"]
my_shape = "my_numbat_plugin:MyShape"
```

NumBAT could discover these using `importlib.metadata.entry_points()`.

Advantages include:

- normal dependency management;
- stable import paths;
- compatibility with child processes;
- clearer installation and distribution;
- less custom import machinery.

## Multiprocessing design recommendation

The current work queue sends deep-copied `Structure` objects to workers. Even
after fixing pickling, this may be expensive and can serialize more state than
the worker actually needs.

A longer-term design would send small, immutable task descriptions and create
worker-local state during worker initialization. For example:

```python
def worker_init(structure_parameters):
    global worker_structure
    worker_structure = build_structure(structure_parameters)

def calculate(task):
    index, wavelength = task
    return worker_structure.calc_EM_modes(...)
```

Potential benefits:

- smaller queue payloads;
- fewer serialization costs;
- no need to pickle mesh and geometry objects for every task;
- explicit ownership of non-thread-safe state;
- opportunity to reuse one structure per worker.

The trade-off is that structure construction must be measured. If mesh
construction is very expensive, constructing once per worker may still be
preferable to repeatedly transferring a complete structure.

## Recommendation

The recommended order is:

1. retain the current custom reduction as a compatibility fix;
2. modernize the file loader and add module caching;
3. migrate built-in templates to stable import paths;
4. add package-based plugin discovery if external templates are important;
5. consider redesigning multiprocessing tasks to avoid transferring complete
   structures.

The current fix is therefore appropriate for the immediate problem, but it
should be viewed as a bridge toward a normal importable plugin architecture,
not necessarily the final design.
